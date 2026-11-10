/*
 * assembler.c - Rebuild a virus genome from scrambled DNA fragments
 *               using a de Bruijn graph.
 *
 * Build (in WSL):   gcc -O2 -o assembler assembler.c
 *
 * Run:
 *   ./assembler data/sars_cov2.fasta                      new random fragments every run
 *   ./assembler data/sars_cov2.fasta -seed 42             repeat an earlier run exactly
 *   ./assembler data/sars_cov2.fasta -cov 30 -err 0.02    change the settings
 *   ./assembler data/sars_cov2.fasta -reads reads.fastq   assemble real or saved fragments
 *
 * The genome file is used for two things only: making the fragments (the
 * fake sequencing machine) and checking the answer at the end. The assembly
 * itself (stages 3-6) sees nothing but the fragments.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>

/* ===================== SETTINGS ===================== */
/* Defaults; each can be changed from the command line. */

static int      K         = 31;    /* chunk (k-mer) length, 2..31            */
static int      READ_LEN  = 100;   /* fragment length                        */
static double   COVERAGE  = 20;    /* times each letter is read, on average  */
static double   ERR_RATE  = 0.01;  /* chance that a letter is misread        */
static int      MIN_COUNT = 0;     /* chunks seen fewer times are typos;
                                      0 = choose automatically               */
static unsigned SEED      = 0;     /* 0 = pick from the clock (new run)      */
static int      STRANDS   = 2;     /* 2 = fragments come from both DNA
                                      strands; 1 = one strand only (needed
                                      for tiny k, like the toy example)     */

#define HIST_MAX 100000            /* counts above this share one bucket     */

static const char BASES[] = "ACGT";

/* ===================== SMALL HELPERS ===================== */

static void *xmalloc(size_t n)
{
    void *p = malloc(n);
    if (!p) { fprintf(stderr, "Out of memory\n"); exit(1); }
    return p;
}

static void *xcalloc(size_t n, size_t size)
{
    void *p = calloc(n, size);
    if (!p) { fprintf(stderr, "Out of memory\n"); exit(1); }
    return p;
}

static void *xrealloc(void *p, size_t n)
{
    p = realloc(p, n);
    if (!p) { fprintf(stderr, "Out of memory\n"); exit(1); }
    return p;
}

static double seconds_since(clock_t start)
{
    return (double)(clock() - start) / CLOCKS_PER_SEC;
}

/* Random number in [0, 1). */
static double rand01(void)
{
    return rand() / ((double)RAND_MAX + 1);
}

/* ===================== DNA LETTERS AND 2-BIT ENCODING ===================== */
/* A=00 C=01 G=10 T=11, so a chunk of up to 32 letters fits in one uint64_t. */

static int is_base(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

static int base_code(char c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    default:  return 3;   /* 'T' */
    }
}

/* The letter on the other DNA strand: A pairs with T, C pairs with G. */
static char complement(char c)
{
    switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default:  return 'N';
    }
}

/* Turn a number back into `len` letters. */
static void decode(uint64_t code, int len, char *out)
{
    for (int i = len - 1; i >= 0; i--) {
        out[i] = BASES[code & 3];
        code >>= 2;
    }
    out[len] = '\0';
}

/*
 * Reverse complement ("mirror"): the same stretch of DNA as read from the
 * other strand. Read it backwards and swap A<->T, C<->G. With our codes the
 * partner of letter b is 3 - b.
 */
static uint64_t revcomp_code(uint64_t code, int len)
{
    uint64_t rc = 0;
    for (int i = 0; i < len; i++) {
        rc = (rc << 2) | (3 - (code & 3));
        code >>= 2;
    }
    return rc;
}

static void revcomp_string(char *s, int len)
{
    for (int i = 0, j = len - 1; i <= j; i++, j--) {
        char a = s[i], b = s[j];
        s[i] = complement(b);
        s[j] = complement(a);
    }
}

/* ===================== HASH TABLE ===================== */
/*
 * Maps a 64-bit key (an encoded chunk) to a 32-bit value (a count, an ID or
 * a position). Open addressing with linear probing: if a slot is taken by a
 * different key, try the next slot. When the table gets half full, it
 * doubles in size (rehashing), so it stays fast for any amount of input.
 */

typedef struct {
    uint64_t      *keys;
    uint32_t      *vals;
    unsigned char *used;   /* 1 if the slot holds a key */
    uint64_t       cap;    /* number of slots, always a power of 2 */
    uint64_t       size;   /* number of keys stored */
} HashTable;

static HashTable ht_create(uint64_t max_keys)
{
    HashTable t;
    t.cap = 1;
    while (t.cap < 2 * max_keys) t.cap <<= 1;
    t.keys = xmalloc(t.cap * sizeof(uint64_t));
    t.vals = xmalloc(t.cap * sizeof(uint32_t));
    t.used = xcalloc(t.cap, 1);
    t.size = 0;
    return t;
}

static void ht_free(HashTable *t)
{
    free(t->keys);
    free(t->vals);
    free(t->used);
}

/* Slot that holds `key`, or the empty slot where it would go. */
static uint64_t ht_slot(const HashTable *t, uint64_t key)
{
    uint64_t h = key * 0x9E3779B97F4A7C15ULL;   /* scramble the bits */
    h ^= h >> 32;
    uint64_t i = h & (t->cap - 1);              /* same as h % cap   */
    while (t->used[i] && t->keys[i] != key)
        i = (i + 1) & (t->cap - 1);             /* linear probing    */
    return i;
}

/* Move every key into a table twice as big. */
static void ht_grow(HashTable *t)
{
    HashTable bigger = ht_create(t->cap);
    for (uint64_t i = 0; i < t->cap; i++) {
        if (!t->used[i]) continue;
        uint64_t j = ht_slot(&bigger, t->keys[i]);
        bigger.used[j] = 1;
        bigger.keys[j] = t->keys[i];
        bigger.vals[j] = t->vals[i];
        bigger.size++;
    }
    ht_free(t);
    *t = bigger;
}

/* Pointer to the value for `key`. A new key is added with value 0, and
 * *is_new tells the caller whether that happened. */
static uint32_t *ht_get_or_add(HashTable *t, uint64_t key, int *is_new)
{
    if (2 * (t->size + 1) > t->cap) ht_grow(t);   /* keep it at most half full */
    uint64_t i = ht_slot(t, key);
    *is_new = !t->used[i];
    if (*is_new) {
        t->used[i] = 1;
        t->keys[i] = key;
        t->vals[i] = 0;
        t->size++;
    }
    return &t->vals[i];
}

/* 1 and the value in *val if `key` is present, else 0. */
static int ht_find(const HashTable *t, uint64_t key, uint32_t *val)
{
    uint64_t i = ht_slot(t, key);
    if (!t->used[i]) return 0;
    *val = t->vals[i];
    return 1;
}

/* ===================== STAGE 1: LOAD THE GENOME ===================== */
/* FASTA: a '>' header line, then lines of letters. Reads the first record;
 * any letter other than A/C/G/T becomes 'N' (unknown). */

static char *read_fasta(const char *path, int *len_out)
{
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    int cap = 1024, len = 0, records = 0;
    char *seq = xmalloc(cap);
    char line[4096];

    while (fgets(line, sizeof line, f)) {
        if (line[0] == '>') {
            if (++records > 1) break;
            continue;
        }
        for (char *p = line; *p; p++) {
            char c = (char)toupper((unsigned char)*p);
            if (!isalpha((unsigned char)c)) continue;
            if (len + 1 >= cap) { cap *= 2; seq = xrealloc(seq, cap); }
            seq[len++] = is_base(c) ? c : 'N';
        }
    }
    fclose(f);
    seq[len] = '\0';
    *len_out = len;
    return seq;
}

/* ===================== STAGE 2: GET FRAGMENTS ===================== */
/* Pretends to be the sequencing machine: copies random pieces of the genome,
 * misreads some letters, and (like a real machine) reads about half of the
 * pieces from the other DNA strand. */

static char **simulate_reads(const char *genome, int g_len, int *n_out,
                             int *errors_out, int *flipped_out)
{
    int n = (int)(COVERAGE * g_len / READ_LEN);
    char **reads = xmalloc(n * sizeof(char *));
    int errors = 0, flipped = 0;

    for (int i = 0; i < n; i++) {
        int start = rand() % (g_len - READ_LEN + 1);
        char *r = xmalloc(READ_LEN + 1);
        memcpy(r, genome + start, READ_LEN);
        r[READ_LEN] = '\0';

        for (int j = 0; j < READ_LEN; j++) {
            if (rand01() < ERR_RATE) {
                /* replace with one of the 3 other letters */
                r[j] = BASES[(base_code(r[j]) + 1 + rand() % 3) % 4];
                errors++;
            }
        }
        if (STRANDS == 2 && rand() % 2) {   /* read from the other strand */
            revcomp_string(r, READ_LEN);
            flipped++;
        }
        reads[i] = r;
    }
    *n_out = n;
    *errors_out = errors;
    *flipped_out = flipped;
    return reads;
}

/*
 * Reads fragments from a file. Understands three layouts:
 *   FASTQ (real sequencer output): groups of 4 lines, letters on the 2nd
 *   FASTA: '>' label lines, skipped
 *   plain: one fragment per line
 * Any letter other than A/C/G/T becomes 'N' (unknown).
 */
static char **load_reads(const char *path, int *n_out)
{
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }

    int cap = 1024, n = 0, line_no = 0, fastq = -1;
    char **reads = xmalloc(cap * sizeof(char *));
    char line[4096];

    while (fgets(line, sizeof line, f)) {
        if (fastq == -1) fastq = (line[0] == '@');    /* decide from the first line */
        line_no++;
        if (fastq ? line_no % 4 != 2 : line[0] == '>') continue;

        char *r = xmalloc(strlen(line) + 1);
        int len = 0;
        for (char *p = line; *p; p++) {
            char c = (char)toupper((unsigned char)*p);
            if (isalpha((unsigned char)c)) r[len++] = is_base(c) ? c : 'N';
        }
        r[len] = '\0';
        if (len == 0) { free(r); continue; }
        if (n == cap) { cap *= 2; reads = xrealloc(reads, cap * sizeof(char *)); }
        reads[n++] = r;
    }
    fclose(f);
    *n_out = n;
    return reads;
}

static void save_reads(const char *path, char **reads, int n)
{
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    for (int i = 0; i < n; i++) fprintf(f, "%s\n", reads[i]);
    fclose(f);
}

/* ===================== STAGE 3: COUNT CHUNKS ===================== */
/*
 * Slides a window of K letters over every fragment, one letter at a time,
 * and counts how often each chunk appears.
 *
 * A chunk and its mirror are the same piece of DNA seen from the two
 * strands, so they are counted together under one name: whichever of the
 * two numbers is smaller (the "canonical" chunk). Windows containing an N
 * are skipped.
 */

static HashTable count_kmers(char **reads, int n_reads)
{
    HashTable counts = ht_create(1024);          /* grows as needed */
    uint64_t mask = (1ULL << (2 * K)) - 1;      /* keeps the last K letters */
    int shift = 2 * (K - 1);                    /* position of the first letter */
    int is_new;

    for (int i = 0; i < n_reads; i++) {
        uint64_t fwd = 0, rev = 0;
        int run = 0;                             /* valid letters in a row */
        for (int j = 0; reads[i][j]; j++) {
            char c = reads[i][j];
            if (!is_base(c)) { run = 0; continue; }   /* N: start again after it */
            int b = base_code(c);
            /* rolling update of the chunk and of its mirror */
            fwd = ((fwd << 2) | b) & mask;
            rev = (rev >> 2) | ((uint64_t)(3 - b) << shift);
            if (++run >= K) {
                uint64_t canon = (STRANDS == 2 && rev < fwd) ? rev : fwd;
                (*ht_get_or_add(&counts, canon, &is_new))++;
            }
        }
    }
    return counts;
}

/* ===================== STAGE 4: REMOVE TYPO CHUNKS ===================== */

/*
 * Typical count of a real chunk: the median count over all sightings of
 * chunks seen at least twice. Chunks seen once are left out because they
 * are almost all typos; real chunks then make up most of the sightings.
 */
static int typical_count(const long *hist)
{
    long total = 0, sum = 0;
    for (int c = 2; c <= HIST_MAX; c++) total += (long)c * hist[c];
    for (int c = 2; c <= HIST_MAX; c++) {
        sum += (long)c * hist[c];
        if (2 * sum >= total) return c;
    }
    return 1;
}

/*
 * Typo chunks make a big pile at count 1 that quickly falls away; real
 * chunks make a second hump at higher counts. The threshold is the bottom
 * of the dip between them: the first count where the histogram rises again.
 * The dip must be a low count (below twice the typical count); with very
 * little data there is no clear dip, and anything seen twice is kept.
 * With very deep data the same typo can turn up many times, so the
 * threshold is also at least a tenth of the typical count.
 */
static int auto_threshold(const long *hist)
{
    int typical = typical_count(hist);
    int dip = 2;
    for (int c = 1; c < 2 * typical && c < HIST_MAX; c++)
        if (hist[c + 1] > hist[c]) { dip = c; break; }
    int tenth = typical / 10;
    return dip > tenth ? dip : tenth;
}

static uint64_t *filter_kmers(const HashTable *counts, int threshold, int *n_out)
{
    uint64_t *solid = xmalloc((counts->size + 1) * sizeof(uint64_t));
    int n = 0;
    for (uint64_t i = 0; i < counts->cap; i++)
        if (counts->used[i] && counts->vals[i] >= (uint32_t)threshold)
            solid[n++] = counts->keys[i];
    *n_out = n;
    return solid;
}

/* ===================== STAGE 5: BUILD THE GRAPH ===================== */
/*
 * de Bruijn graph: every chunk is an EDGE from its first K-1 letters to its
 * last K-1 letters. The (K-1)-letter pieces are the NODES, numbered 0, 1, 2...
 * using a second hash table.
 *
 * Each kept chunk is added twice, as read and as its mirror, so the graph
 * holds both strands. Edge 2i is a chunk and edge 2i+1 its mirror, so the
 * mirror of edge e is always e ^ 1 (flip the last bit). With -strands 1,
 * each chunk is added once.
 *
 * Each node's outgoing edges form a linked list stored in arrays:
 * first_edge[node] is the first edge, next_edge[edge] the one after it,
 * and -1 means "no more".
 */

typedef struct {
    int            n_nodes, n_edges;
    uint64_t      *edge_kmer;    /* the chunk this edge stands for   */
    int           *edge_from;    /* node the edge leaves             */
    int           *edge_to;      /* node the edge enters             */
    int           *next_edge;    /* next outgoing edge of edge_from  */
    int           *first_edge;   /* first outgoing edge of a node    */
    int           *in_deg, *out_deg;
    unsigned char *used;         /* 1 once the edge or its mirror is walked */
} Graph;

static int node_id(HashTable *ids, uint64_t code, int *n_nodes)
{
    int is_new;
    uint32_t *id = ht_get_or_add(ids, code, &is_new);
    if (is_new) *id = (uint32_t)(*n_nodes)++;
    return (int)*id;
}

static Graph build_graph(const uint64_t *solid, int n_solid)
{
    Graph g;
    int n_edges   = STRANDS * n_solid; /* every chunk in both directions */
    int max_nodes = 2 * n_edges + 1;   /* each edge adds at most 2 nodes */
    g.n_nodes    = 0;
    g.n_edges    = n_edges;
    g.edge_kmer  = xmalloc(n_edges * sizeof(uint64_t));
    g.edge_from  = xmalloc(n_edges * sizeof(int));
    g.edge_to    = xmalloc(n_edges * sizeof(int));
    g.next_edge  = xmalloc(n_edges * sizeof(int));
    g.used       = xcalloc(n_edges, 1);
    g.first_edge = xmalloc(max_nodes * sizeof(int));
    g.in_deg     = xcalloc(max_nodes, sizeof(int));
    g.out_deg    = xcalloc(max_nodes, sizeof(int));
    for (int v = 0; v < max_nodes; v++) g.first_edge[v] = -1;   /* all lists empty */

    HashTable ids = ht_create(max_nodes);
    uint64_t suffix_mask = (1ULL << (2 * (K - 1))) - 1;

    for (int e = 0; e < n_edges; e++) {
        uint64_t kmer = solid[e / STRANDS];
        if (e % STRANDS == 1) kmer = revcomp_code(kmer, K);       /* the mirror */
        int from = node_id(&ids, kmer >> 2, &g.n_nodes);          /* first K-1 letters */
        int to   = node_id(&ids, kmer & suffix_mask, &g.n_nodes); /* last K-1 letters  */
        g.edge_kmer[e] = kmer;
        g.edge_from[e] = from;
        g.edge_to[e]   = to;
        g.next_edge[e] = g.first_edge[from];   /* add to the front of the list */
        g.first_edge[from] = e;
        g.out_deg[from]++;
        g.in_deg[to]++;
    }
    ht_free(&ids);
    return g;
}

/* ===================== STAGE 6: WALK THE GRAPH ===================== */
/*
 * Hierholzer's algorithm: walks the graph using every edge exactly once
 * (an Eulerian path). Walks start at nodes with more edges out than in,
 * which is where the genome begins. Walking an edge also marks its mirror
 * as used, so each stretch of DNA is spelled once, on one strand.
 *
 * Writes edge numbers into `order` in walking order, with -1 before each
 * new walk. Returns how many entries were written.
 */

/* First unused edge leaving v, or -1. Used edges at the front of the list
 * are dropped so they are never looked at again. */
static int next_unused(Graph *g, int v)
{
    int e = g->first_edge[v];
    while (e != -1 && g->used[e]) e = g->next_edge[e];
    g->first_edge[v] = e;
    return e;
}

static int eulerian_walks(Graph *g, int *order)
{
    int *stack_node = xmalloc((g->n_edges + 1) * sizeof(int));
    int *stack_edge = xmalloc((g->n_edges + 1) * sizeof(int));  /* edge used to arrive */
    int *walk       = xmalloc((g->n_edges + 1) * sizeof(int));
    int  n_order    = 0;

    /* pass 0: only proper start nodes; pass 1: whatever edges are left */
    for (int pass = 0; pass < 2; pass++) {
        for (int s = 0; s < g->n_nodes; s++) {
            if (next_unused(g, s) == -1) continue;
            if (pass == 0 && g->out_deg[s] <= g->in_deg[s]) continue;

            int top = 0, n_walk = 0;
            stack_node[0] = s;
            stack_edge[0] = -1;

            while (top >= 0) {
                int v = stack_node[top];
                int e = next_unused(g, v);
                if (e != -1) {
                    /* unused edge: take it (and retire its mirror), move forward */
                    g->used[e] = 1;
                    if (STRANDS == 2) g->used[e ^ 1] = 1;
                    g->first_edge[v] = g->next_edge[e];
                    top++;
                    stack_node[top] = g->edge_to[e];
                    stack_edge[top] = e;
                } else {
                    /* stuck: step back, recording the edge we came by */
                    if (stack_edge[top] != -1) walk[n_walk++] = stack_edge[top];
                    top--;
                }
            }

            /* edges come off the stack in reverse order */
            order[n_order++] = -1;
            for (int i = n_walk - 1; i >= 0; i--) order[n_order++] = walk[i];
        }
    }
    free(stack_node);
    free(stack_edge);
    free(walk);
    return n_order;
}

/* ===================== STAGE 6b: EDGES -> DNA (CONTIGS) ===================== */
/*
 * Splits the walks wherever two edges in a row don't connect, and spells each
 * piece: the first chunk in full, then the last letter of every next chunk.
 * Each piece is a contig: a stretch of genome we are sure about.
 */

typedef struct {
    char *seq;
    int   len;
} Contig;

static int build_contigs(const Graph *g, const int *order, int n_order, Contig *contigs)
{
    int n = 0;
    for (int i = 0; i < n_order; i++) {
        if (order[i] == -1) continue;

        /* find where this connected run ends */
        int j = i + 1;
        while (j < n_order && order[j] != -1 &&
               g->edge_from[order[j]] == g->edge_to[order[j - 1]])
            j++;

        int len = K + (j - i - 1);
        char *s = xmalloc(len + 1);
        decode(g->edge_kmer[order[i]], K, s);
        for (int t = i + 1; t < j; t++)
            s[K + (t - i - 1)] = BASES[g->edge_kmer[order[t]] & 3];
        s[len] = '\0';

        contigs[n].seq = s;
        contigs[n].len = len;
        n++;
        i = j - 1;
    }
    return n;
}

static int by_length_desc(const void *a, const void *b)
{
    return ((const Contig *)b)->len - ((const Contig *)a)->len;
}

/* N50: sort longest first, add lengths until you reach half the total. */
static int n50(const Contig *contigs, int n, long total)
{
    long sum = 0;
    for (int i = 0; i < n; i++) {
        sum += contigs[i].len;
        if (2 * sum >= total) return contigs[i].len;
    }
    return 0;
}

/* ===================== STAGE 7: CHECK THE ANSWER ===================== */

/* Edit distance: fewest single-letter changes (substitute, insert, delete)
 * turning a into b. An N in b (unknown letter) matches anything.
 * Dynamic programming, keeping only two rows. */
static int edit_distance(const char *a, int n, const char *b, int m)
{
    int *prev = xmalloc((m + 1) * sizeof(int));
    int *cur  = xmalloc((m + 1) * sizeof(int));

    for (int j = 0; j <= m; j++) prev[j] = j;
    for (int i = 1; i <= n; i++) {
        cur[0] = i;
        for (int j = 1; j <= m; j++) {
            int sub = prev[j - 1] + (a[i - 1] != b[j - 1] && b[j - 1] != 'N');
            int del = prev[j] + 1;
            int ins = cur[j - 1] + 1;
            int best = sub < del ? sub : del;
            cur[j] = best < ins ? best : ins;
        }
        int *tmp = prev; prev = cur; cur = tmp;
    }
    int result = prev[m];
    free(prev);
    free(cur);
    return result;
}

/* Where does contig s start in the genome? Uses its first chunk that also
 * appears in the genome. Returns 0 if no chunk does. */
static int place(const HashTable *where, const char *s, int len, int *start)
{
    uint64_t mask = (1ULL << (2 * K)) - 1, code = 0;
    for (int i = 0; i < len; i++) {
        code = ((code << 2) | base_code(s[i])) & mask;
        uint32_t pos;
        if (i >= K - 1 && ht_find(where, code, &pos)) {
            *start = (int)pos - (i - K + 1);
            return 1;
        }
    }
    return 0;
}

/*
 * For each contig: find where it belongs in the genome (trying both strands),
 * then compare it letter by letter against that stretch of genome with edit
 * distance. Also lists the genome positions where the contigs differ.
 */
static void evaluate(const char *genome, int g_len, const Contig *contigs, int n_contigs)
{
    /* index: chunk -> its position in the real genome */
    HashTable where = ht_create(g_len);
    uint64_t mask = (1ULL << (2 * K)) - 1, code = 0;
    int is_new, run = 0;
    for (int p = 0; p < g_len; p++) {
        if (!is_base(genome[p])) { run = 0; continue; }   /* skip chunks with N */
        code = ((code << 2) | base_code(genome[p])) & mask;
        if (++run >= K) {
            uint32_t *pos = ht_get_or_add(&where, code, &is_new);
            if (is_new) *pos = (uint32_t)(p - K + 1);
        }
    }

    unsigned char *covered = xcalloc(g_len, 1);
    char *diff_letter = xcalloc(g_len, 1);   /* our letter where it differs, else 0 */
    long placed_len = 0, total_edits = 0;
    int unplaced = 0, flipped = 0;

    for (int c = 0; c < n_contigs; c++) {
        int len = contigs[c].len, start;
        char *s = contigs[c].seq, *mirror = NULL;

        if (!place(&where, s, len, &start)) {
            /* maybe it was spelled from the other strand */
            mirror = xmalloc(len + 1);
            memcpy(mirror, s, len + 1);
            revcomp_string(mirror, len);
            s = mirror;
            if (!place(&where, s, len, &start)) { unplaced++; free(mirror); continue; }
            flipped++;
        }

        int skip = start < 0 ? -start : 0;   /* part hanging off the genome start */
        start += skip;
        int ref_len = len - skip < g_len - start ? len - skip : g_len - start;
        int edits = edit_distance(s + skip, len - skip, genome + start, ref_len);
        total_edits += edits;
        placed_len  += len;

        int mismatches = 0;
        for (int p = 0; p < ref_len; p++)
            mismatches += s[skip + p] != genome[start + p] && genome[start + p] != 'N';
        for (int p = 0; p < ref_len; p++) {
            covered[start + p] = 1;
            /* list single-letter differences, unless an insertion or
               deletion has shifted the letters out of line */
            if (mismatches <= edits && s[skip + p] != genome[start + p] &&
                genome[start + p] != 'N')
                diff_letter[start + p] = s[skip + p];
        }
        free(mirror);
    }

    /* N letters in the genome are unknown, so they are left out of the score */
    long n_covered = 0, n_known = 0;
    for (int p = 0; p < g_len; p++) {
        if (genome[p] == 'N') continue;
        n_known++;
        n_covered += covered[p];
    }

    if (n_known < g_len)
        printf("  Unknown (N) letters in genome      : %ld  (not scored)\n", g_len - n_known);
    printf("  Genome recovered (genome fraction) : %.2f%%  (%ld of %ld letters)\n",
           100.0 * n_covered / n_known, n_covered, n_known);
    printf("  Letters matching (identity)        : %.4f%%  (%ld differences in %ld letters)\n",
           placed_len ? 100.0 * (placed_len - total_edits) / placed_len : 0.0,
           total_edits, placed_len);
    printf("  Contigs on the other strand        : %d  (flipped before comparing)\n", flipped);
    printf("  Contigs not found in genome        : %d\n", unplaced);

    int shown = 0, n_diff = 0;
    for (int p = 0; p < g_len; p++) n_diff += diff_letter[p] != 0;
    if (n_diff) {
        printf("  Single-letter differences          : %d at genome positions:\n   ", n_diff);
        for (int p = 0; p < g_len && shown < 30; p++) {
            if (!diff_letter[p]) continue;
            printf(" %d:%c>%c", p + 1, genome[p], diff_letter[p]);
            if (++shown % 8 == 0) printf("\n   ");
        }
        if (n_diff > shown) printf(" ... (%d more)", n_diff - shown);
        printf("\n");
    }

    /* show which parts of the genome are missing */
    int gaps = 0;
    for (int p = 0; p < g_len; p++) {
        if (covered[p] || genome[p] == 'N') continue;
        int q = p;
        while (q < g_len && !covered[q] && genome[q] != 'N') q++;
        if (gaps < 5)
            printf("  Missing stretch                    : positions %d-%d (%d letters)\n",
                   p + 1, q, q - p);
        gaps++;
        p = q;
    }
    if (gaps > 5) printf("  ... and %d more missing stretches\n", gaps - 5);

    free(covered);
    free(diff_letter);
    ht_free(&where);
}

/* ===================== MAIN ===================== */

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage: %s genome.fasta [options]\n"
        "  -k N        chunk length, 2-31            (default 31)\n"
        "  -len N      fragment length               (default 100)\n"
        "  -cov X      coverage                      (default 20)\n"
        "  -err X      typo rate, e.g. 0.01 = 1%%     (default 0.01)\n"
        "  -min N      minimum chunk count to keep   (default: chosen automatically)\n"
        "  -seed N     repeat a run exactly          (default: new every run)\n"
        "  -reads FILE assemble fragments from FILE (FASTQ, FASTA or one per line)\n"
        "  -strands N  2 = fragments from both DNA strands, 1 = one strand (default 2)\n", prog);
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc < 2) usage(argv[0]);
    const char *reads_path = NULL;

    for (int i = 2; i < argc; i++) {
        if (i + 1 >= argc) usage(argv[0]);
        if      (!strcmp(argv[i], "-k"))     K         = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-len"))   READ_LEN  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-cov"))   COVERAGE  = atof(argv[++i]);
        else if (!strcmp(argv[i], "-err"))   ERR_RATE  = atof(argv[++i]);
        else if (!strcmp(argv[i], "-min"))   MIN_COUNT = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-seed"))  SEED      = (unsigned)strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "-reads")) reads_path = argv[++i];
        else if (!strcmp(argv[i], "-strands")) STRANDS = atoi(argv[++i]);
        else usage(argv[0]);
    }
    if (K < 2 || K > 31) { fprintf(stderr, "k must be between 2 and 31\n"); return 1; }
    if (STRANDS != 1 && STRANDS != 2) { fprintf(stderr, "-strands must be 1 or 2\n"); return 1; }
    if (SEED == 0) SEED = (unsigned)time(NULL);
    srand(SEED);

    clock_t t_all = clock(), t;
    char preview[64];

    /* ---- Stage 1 ---- */
    printf("\n=== Stage 1: Load genome ===\n");
    t = clock();
    int g_len;
    char *genome = read_fasta(argv[1], &g_len);
    if (g_len == 0) { fprintf(stderr, "No DNA letters found in %s\n", argv[1]); return 1; }
    printf("  File          : %s\n", argv[1]);
    printf("  Length        : %d letters\n", g_len);
    printf("  Starts with   : %.60s...\n", genome);
    printf("  Time          : %.3f s\n", seconds_since(t));

    /* ---- Stage 2 ---- */
    int n_reads;
    char **reads;
    t = clock();
    if (reads_path) {
        printf("\n=== Stage 2: Load fragments ===\n");
        reads = load_reads(reads_path, &n_reads);
        printf("  File          : %s\n", reads_path);
    } else {
        printf("\n=== Stage 2: Make fragments (simulated sequencing machine) ===\n");
        if (READ_LEN <= K || READ_LEN > g_len) {
            fprintf(stderr, "Fragment length must be more than k and at most the genome length\n");
            return 1;
        }
        int errors, flipped;
        reads = simulate_reads(genome, g_len, &n_reads, &errors, &flipped);
        save_reads("fragments.txt", reads, n_reads);
        printf("  Seed          : %u   (rerun with -seed %u to get these same fragments)\n", SEED, SEED);
        printf("  Settings      : length %d, coverage %.1fx, typo rate %.2f%%\n",
               READ_LEN, COVERAGE, ERR_RATE * 100);
        printf("  Typos added   : %d\n", errors);
        printf("  Other strand  : %d fragments read from the other strand\n", flipped);
        printf("  Saved to      : fragments.txt\n");
    }
    long total_letters = 0;
    int min_len = n_reads ? (int)strlen(reads[0]) : 0, max_len = 0, with_n = 0;
    for (int i = 0; i < n_reads; i++) {
        int len = (int)strlen(reads[i]);
        total_letters += len;
        if (len < min_len) min_len = len;
        if (len > max_len) max_len = len;
        if (strchr(reads[i], 'N')) with_n++;
    }
    printf("  Fragments     : %d  (lengths %d-%d, %d contain N)\n", n_reads, min_len, max_len, with_n);
    printf("  Letters       : %ld  (about %.0fx coverage)\n", total_letters, (double)total_letters / g_len);
    for (int i = 0; i < 3 && i < n_reads; i++)
        printf("  Fragment %d    : %.60s%s\n", i + 1, reads[i], strlen(reads[i]) > 60 ? "..." : "");
    printf("  Time          : %.3f s\n", seconds_since(t));

    /* ---- Stage 3 ---- */
    printf("\n=== Stage 3: Cut fragments into %d-letter chunks and count them ===\n", K);
    t = clock();
    HashTable counts = count_kmers(reads, n_reads);
    long *hist = xcalloc(HIST_MAX + 1, sizeof(long));   /* hist[c] = chunks seen c times */
    for (uint64_t i = 0; i < counts.cap; i++)
        if (counts.used[i]) hist[counts.vals[i] > HIST_MAX ? HIST_MAX : counts.vals[i]]++;
    long more = 0;
    for (int c = 11; c <= HIST_MAX; c++) more += hist[c];
    printf("  Different chunks : %llu  (a chunk and its mirror count as one)\n",
           (unsigned long long)counts.size);
    printf("  How many chunks were seen N times:\n");
    for (int c = 1; c <= 10; c++) printf("    N = %2d  : %ld\n", c, hist[c]);
    printf("    N = 11+ : %ld\n", more);
    printf("  Time          : %.3f s\n", seconds_since(t));

    /* ---- Stage 4 ---- */
    int threshold = MIN_COUNT ? MIN_COUNT : auto_threshold(hist);
    printf("\n=== Stage 4: Remove chunks seen fewer than %d times (typos) ===\n", threshold);
    t = clock();
    int n_solid;
    uint64_t *solid = filter_kmers(&counts, threshold, &n_solid);
    printf("  Typical count of a real chunk : %d\n", typical_count(hist));
    printf("  Threshold     : %d  (%s)\n", threshold,
           MIN_COUNT ? "set with -min"
                     : "automatic: the dip in the histogram, at least a tenth of the typical count");
    printf("  Kept          : %d\n", n_solid);
    printf("  Removed       : %llu\n", (unsigned long long)(counts.size - n_solid));
    printf("  Time          : %.3f s\n", seconds_since(t));
    ht_free(&counts);
    free(hist);

    /* ---- Stage 5 ---- */
    printf("\n=== Stage 5: Build de Bruijn graph (both strands) ===\n");
    t = clock();
    Graph g = build_graph(solid, n_solid);
    int starts = 0, ends = 0, branches = 0;
    for (int v = 0; v < g.n_nodes; v++) {
        if (g.out_deg[v] > g.in_deg[v]) starts++;
        if (g.in_deg[v] > g.out_deg[v]) ends++;
        if (g.out_deg[v] > 1) branches++;
    }
    printf("  Nodes (%d-letter pieces) : %d\n", K - 1, g.n_nodes);
    printf("  Edges (chunks x 2)       : %d\n", g.n_edges);
    printf("  Start-like nodes (out > in) : %d\n", starts);
    printf("  End-like nodes   (in > out) : %d\n", ends);
    printf("  Branching nodes  (out > 1)  : %d\n", branches);
    printf("  Time          : %.3f s\n", seconds_since(t));

    /* ---- Stage 6 ---- */
    printf("\n=== Stage 6: Walk the graph (Hierholzer) and spell out contigs ===\n");
    t = clock();
    int *order = xmalloc((2 * g.n_edges + 1) * sizeof(int));
    int n_order = eulerian_walks(&g, order);
    Contig *contigs = xmalloc((g.n_edges + 1) * sizeof(Contig));
    int n_contigs = build_contigs(&g, order, n_order, contigs);
    qsort(contigs, n_contigs, sizeof(Contig), by_length_desc);

    long total_len = 0;
    for (int i = 0; i < n_contigs; i++) total_len += contigs[i].len;

    FILE *out = fopen("contigs.fasta", "w");
    if (out) {
        for (int i = 0; i < n_contigs; i++)
            fprintf(out, ">contig_%d length=%d\n%s\n", i + 1, contigs[i].len, contigs[i].seq);
        fclose(out);
    }

    printf("  Contigs       : %d  (saved to contigs.fasta)\n", n_contigs);
    printf("  Longest       : %d letters\n", n_contigs ? contigs[0].len : 0);
    printf("  Total length  : %ld letters\n", total_len);
    printf("  N50           : %d\n", n50(contigs, n_contigs, total_len));
    for (int i = 0; i < 5 && i < n_contigs; i++) {
        snprintf(preview, sizeof preview, "%.40s", contigs[i].seq);
        printf("  Contig %d      : %d letters, %s%s\n", i + 1, contigs[i].len, preview,
               contigs[i].len > 40 ? "..." : "");
    }
    printf("  Time          : %.3f s\n", seconds_since(t));

    /* ---- Stage 7 ---- */
    printf("\n=== Stage 7: Compare with the reference genome ===\n");
    t = clock();
    evaluate(genome, g_len, contigs, n_contigs);
    printf("  Time          : %.3f s\n", seconds_since(t));

    printf("\nTotal time: %.3f s\n\n", seconds_since(t_all));

    /* ---- clean up ---- */
    for (int i = 0; i < n_contigs; i++) free(contigs[i].seq);
    free(contigs);
    free(order);
    free(g.edge_kmer); free(g.edge_from); free(g.edge_to); free(g.next_edge);
    free(g.used); free(g.first_edge); free(g.in_deg); free(g.out_deg);
    free(solid);
    for (int i = 0; i < n_reads; i++) free(reads[i]);
    free(reads);
    free(genome);
    return 0;
}
