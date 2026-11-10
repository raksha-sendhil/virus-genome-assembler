# Rebuilding a Virus Genome from Scrambled DNA Fragments using de Bruijn Graphs

DSA Experiential Learning project, written in C.

We take millions of short DNA fragments of a virus (SARS-CoV-2) and rebuild
the genome using only those fragments. The fragments are either real output of
a sequencing machine (public 2020 patient samples) or simulated by cutting the
reference genome at random. Because the correct genome is known in both cases,
we can measure exactly how accurate the rebuilt genome is.

A full walkthrough of the code, from zero, is in [EXPLANATION.md](EXPLANATION.md).

## Quick start (WSL / Linux)

```bash
gcc -O2 -o assembler assembler.c

bash demo.sh                       # main demo: pick a real sample, rebuild it, grade it
bash get_sample.sh ERR4204296      # download any real sample + its official genome

./assembler data/sars_cov2.fasta   # simulated fragments, new random cut every run
```

The real fragment files (60-100 MB each) are not stored in this repository.
`demo.sh` and `get_sample.sh` download them from the European Nucleotide
Archive (ENA) the first time a sample is used.

---

## 1. Problem statement

DNA sequencing machines cannot read a whole genome from start to finish.
Instead, they break DNA into millions of short, overlapping fragments
("reads", typically 100-150 letters long) and read each fragment separately.
Where each fragment came from is lost, and some letters are misread.
Software must rebuild the original genome from this pile of fragments.
This is called **genome assembly**.

Assembly matters most when the organism is new and there is no existing
genome to compare against. In January 2020, the first SARS-CoV-2 genome was
assembled from patient samples (Wu et al., *Nature*, 2020, mainly using the
de Bruijn graph assembler MEGAHIT). The sequence was published around 10-11
January 2020, and mRNA vaccine designs were produced within days, based
entirely on that genome.

Our project implements, from scratch in C, the core algorithm used by modern
short-read assemblers (Velvet, SPAdes, MEGAHIT): the **de Bruijn graph**
approach.

---

## 2. Biology background (what the data represents)

- One swab from one infected person contains **millions of copies** of the
  virus. The machine breaks all of those copies into fragments at random
  places. Because each copy breaks at different points, the fragments
  overlap, and that is what lets us put them back together.
  (Like many photocopies of the same book, each shredded differently.)
- The genome is **one specific sequence** (~30,000 letters for SARS-CoV-2)
  that every copy of the virus shares. It is not a "recurring pattern" across
  people.
- Different patients can carry slightly different versions of the virus
  (**variants**). Scientists assemble each sample separately, then compare
  the genomes to find mutations. That comparison is a separate step after
  assembly and is outside our scope.
- A real swab also contains the patient's own DNA, which real pipelines filter
  out first. Our simulated fragments come only from the virus genome.

---

## 3. Data source

Free, public genomes from **NCBI** (https://www.ncbi.nlm.nih.gov/labs/virus/
or GenBank), downloaded as **FASTA** files: one header line, then the DNA
letters (A, C, G, T).

```
>NC_045512.2 Severe acute respiratory syndrome coronavirus 2, complete genome
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA...
```

| Virus | Accession | Length | Use |
|---|---|---|---|
| SARS-CoV-2 | NC_045512.2 | ~29.9 kb | Main genome |
| Phage phiX174 | NC_001422 | ~5.4 kb | Small, quick tests |
| HIV-1 | NC_001802 | ~9.2 kb | Medium test |

In C, the genome is simply a `char` array.

---

## 4. What the system does

1. Loads a real public genome from a FASTA file.
2. Simulates a sequencing machine: cuts the genome into random overlapping
   fragments with configurable length, coverage and error rate.
3. Reassembles the genome using **only** those fragments.
4. Compares the rebuilt genome with the original and reports accuracy.
5. Shows everything in a web dashboard (HTML/CSS/JS) served by the C program.

**Important:** the original genome is used only to make the fragments and to
check the answer at the end. The assembly algorithm never sees it. In real
life nobody has the genome; recovering it is the whole goal. Using it during
assembly would be cheating, like solving a jigsaw puzzle by looking at the
picture on the box.

---

## 5. Implementation pipeline

### Stage 1 - Genome loader
- Parses a FASTA file (header line + DNA sequence) into a `char` buffer.
- Skips the header and newlines; converts lowercase to uppercase.

### Stage 2 - Read simulator (the "sequencing machine")
- Picks random start positions and copies substrings of length L.
- Number of reads = (coverage × genome length) / L.
- Injects random substitution errors at a configurable rate (e.g. 1%).
- Seeded random numbers so experiments can be repeated.
- This stage is **not** part of the assembly algorithm. It stands in for the
  sequencing machine, because we don't have one. The read length is set by
  the machine in real life; we can't choose it.

### Stage 3 - k-mer counting (cutting reads into chunks)
- Slides a window of size k (k = 31) over every read, **one letter at a
  time**. A read of length L gives L - k + 1 k-mers (100 letters → 70 k-mers).
- The chunks overlap; we do **not** split a read into separate side-by-side
  pieces.
- Each k-mer is packed into a 64-bit integer, 2 bits per base
  (A=00, C=01, G=10, T=11), using a rolling shift-and-mask update.
  k ≤ 32 fits in one `uint64_t`.
- Counts how often each k-mer appears in a hash table.

### Stage 4 - Error filtering
- k-mers caused by sequencing errors appear rarely (count 1-2).
  Real k-mers appear many times because of coverage.
- k-mers below a count threshold are deleted.
- Without this, a typo creates a false branch in the graph.

### Stage 5 - de Bruijn graph construction
- Each surviving k-mer becomes a directed **edge** from its prefix (first
  k-1 letters) to its suffix (last k-1 letters).
- The (k-1)-mers are the **nodes**, mapped to integer IDs using a second hash
  table.
- Stored as an adjacency list, with in-degree and out-degree arrays.

### Stage 6 - Eulerian path (Hierholzer's algorithm)
- Rebuilding the genome = finding a path that uses every edge exactly once
  (an Eulerian path).
- The start node is the one where out-degree - in-degree = 1.
- Implemented iteratively with an explicit stack, in O(V + E). Recursion would
  risk a stack overflow on graphs with tens of thousands of edges.
- If repeats or leftover errors break the graph, the program outputs
  **contigs** instead: the longest unambiguous stretches of sequence.
  (Real assemblers also output contigs rather than one single path.)

### Stage 7 - Sequence reconstruction and evaluation
- The path is turned back into a DNA string: the first node, then the last
  letter of each following node.
- Accuracy = edit distance (Levenshtein) between the rebuilt and original
  genomes, computed with dynamic programming using two rows (O(n) memory).
- Reports: % identity, number of contigs, N50, time taken, memory used.
- **N50**: sort contigs longest first and add up their lengths; N50 is the
  length of the contig at which the total reaches half the assembly length.
  Higher is better.

### Stage 8 - Web dashboard
- The C program runs a small HTTP server using **Mongoose** (an embedded
  single-file C library) and exposes JSON endpoints:
  `/api/simulate`, `/api/assemble`, `/api/stats`.
- Frontend (HTML/CSS/JS): sliders for k, read length, coverage and error
  rate; a graph view for small genomes; original vs rebuilt sequence with
  differences highlighted; result charts.

---

## 6. Worked example (toy genome)

Say the genome is only `ACGTTGCA` (the program pretends not to know this).
Fragments are 5 letters and chunks are 3 letters (31 in the real program).

**Step 1 - Make fragments.** Random 5-letter pieces; some overlap, some repeat,
one has a typo:

```
ACGTT  ACGTT  CGTTG  GTTGC  TTGCA  TTGCA  GTAGC  ← typo (should be GTTGC)
```

From here on, the program sees only this pile.

**Step 2 - Cut into 3-letter chunks (sliding window) and count them.**

| Chunk | ACG | CGT | GTT | TTG | TGC | GCA | GTA | TAG | AGC |
|---|---|---|---|---|---|---|---|---|---|
| Count | 2 | 3 | 4 | 4 | 3 | 2 | 1 | 1 | 1 |

Chunks from the typo appear only once, so `GTA`, `TAG` and `AGC` are deleted.

**Step 3 - Link each chunk to the next one (graph).** Neighbouring chunks
share 2 letters. After `ACG`, the next chunk must start with `CG`; the hash
table finds `CGT`.

```
ACG → CGT → GTT → TTG → TGC → GCA
```

If `GTA` had not been deleted, `CGT` would have two possible next chunks
(`GTT` and `GTA`), creating a false branch.

**Step 4 - Walk the graph from start to end.** Start at the chunk nothing
points to (`ACG`) and add the last letter of each following chunk:

```
ACG + T + T + G + C + A = ACGTTGCA ✓
```

**Step 5 - Check the result.** Edit distance against the original gives the
accuracy percentage.

---

## 7. Key questions and answers

**Why cut the fragments into chunks if our job is to put things back together?**
Fragments overlap by any amount, at any position. Joining them directly means
comparing every pair of fragments at every possible overlap length: for
~3,000 reads that is ~9 million pairs × up to 99 overlap lengths, and real
runs have millions of reads. Cutting every fragment into fixed-size chunks
turns "how much do these overlap, and where?" into "have I seen this exact
chunk before?", which a hash table answers instantly (~210,000 operations for
our data). Two more benefits:
- A typo only spoils the ~31 chunks that contain it, not the whole read.
- Duplicate reads collapse into one hash table entry with a higher count.

**Why not cut the genome into chunks directly?**
Because in real life nobody has the genome. We only have it in this project to
make realistic fragments and check the answer.

**Why do neighbouring chunks overlap by k - 1 letters?**
Because of how we cut them, not because of biology. The window moves one letter
at a time, so each chunk drops one letter at the front and adds one at the
back; neighbours always share k - 1 letters. To find what follows a chunk, take
its last k - 1 letters and try the only 4 possibilities (+A, +C, +G, +T): at
most 4 hash table lookups.

**Why 3 in the example but 31 in the real program?**
There are only 4³ = 64 possible 3-letter chunks, so in a 30,000-letter genome
every chunk would appear in hundreds of places and the graph would become a
tangle. There are 4³¹ ≈ 4.6 × 10¹⁸ possible 31-letter chunks, so almost every
chunk appears in only one place. Chunk size rule: small enough to find
overlaps, large enough to be unique.

**Where is the de Bruijn graph?**
It is the graph from Step 3. Formally, the (k-1)-letter overlaps (`AC`, `CG`,
`GT`, ...) are the nodes and each chunk is an edge between them.

**Is this the standard approach?**
Yes, for assembling a new genome from short reads. It was proposed by Idury &
Waterman (1995) and made popular by Pevzner, Tang & Waterman (PNAS, 2001), and
it is the basis of Velvet (2008), SPAdes (2012) and MEGAHIT (2015). Notes:
- For long reads (Nanopore, PacBio), tools mostly use the older
  Overlap-Layout-Consensus method or variants of it.
- When a reference genome already exists, labs usually line reads up against
  it (reference-guided assembly). De Bruijn assembly matters most for new
  pathogens.
- Real assemblers also clean the graph (remove dead-end tips and error
  bubbles) and output contigs.

---

## 8. Data structures and algorithms used

| Concept | Where it is used |
|---|---|
| Hash table (custom hash on 64-bit keys, open addressing or chaining, load factor, resizing) | k-mer → count; (k-1)-mer → node ID |
| Directed multigraph (adjacency list, in/out-degree arrays) | de Bruijn graph |
| Stack | Iterative Hierholzer's algorithm |
| Eulerian path (Hierholzer's algorithm), graph traversal | Rebuilding the genome; extracting contigs |
| Dynamic programming (space-optimised edit distance) | Accuracy check |
| Strings and bit manipulation (sliding window, 2-bit encoding, rolling update) | k-mer extraction |
| Dynamic arrays (`realloc`) | Reads, adjacency lists, output |
| Sorting | Contigs by length for N50 |
| Randomisation (seeded) | Read simulation and error injection |

**Complexity**
- k-mer counting: O(N × L) for N reads of length L (average O(1) hash operations).
- Graph construction and Eulerian path: O(V + E).
- Accuracy check: O(n × m) time, O(m) memory.

---

## 9. C programming concepts used

- Structs, pointers, pointers to structs
- Dynamic memory: `malloc`, `calloc`, `realloc`, `free` (with leak checking)
- Bitwise operators `<<`, `>>`, `&`, `|` for k-mer encoding
- Fixed-width integers (`stdint.h`: `uint64_t`, `uint32_t`)
- File I/O: `fopen`, `fgets`, `fread` for FASTA parsing
- Random numbers: `rand`, `srand` with seeds
- Multi-file design with `.h` / `.c` files and a Makefile
- Command-line arguments (`argc` / `argv`) for a CLI mode alongside the server
- HTTP serving via Mongoose
- Timing with `clock()`

---

## 10. Project layout

Everything is in one C file. The terminal version is done; the web dashboard
comes later, using the hand-written socket server from last year's C lab
project (no Mongoose). See [EXPLANATION.md](EXPLANATION.md) for a full
walkthrough of the code.

```
dsa_el/
├── assembler.c                all stages 1-7 (build: gcc -O2 -o assembler assembler.c)
├── demo.sh                    main demo: choose a sample, rebuild, compare with the answer
├── get_sample.sh              downloads a real sample and its official genome from ENA
├── README.md                  project plan (this file)
├── EXPLANATION.md             how the code works, from zero
├── data/
│   ├── sars_cov2.fasta        NC_045512.2 from NCBI
│   ├── toy.fasta              8-letter toy genome (ACGTTGCA)
│   ├── toy_fragments.txt      the 7 toy fragments, one with a typo
│   ├── other/                 genomes of six other viruses (phiX174, hepatitis B,
│   │                          HIV-1, dengue 2, Ebola, mpox) for simulated tests
│   └── real/                  real sequencing runs from ENA (COG-UK, 2020)
│       ├── <run>.fastq                  the fragments (downloaded, not in git)
│       ├── <run>.consensus.fasta        the sample's official genome (the answer)
│       ├── <run>.info                   lab method and collection date
│       └── <run>.contigs.fasta          our assembler's output
├── results/                   saved outputs of every test run
├── fragments.txt              written by each simulated run: the fragments it made
└── contigs.fasta              written by each run: the rebuilt DNA
```

Later, for the dashboard: `index.html`, `style.css`, `app.js`.

---

## 11. Experiments (report component)

We measure how assembly quality changes with:
- **k-mer size:** small k tangles the graph because of repeats; large k breaks
  it apart because of errors and gaps.
- **Coverage:** too little coverage leaves gaps in the genome.
- **Error rate:** more errors create false branches in the graph.

Results are plotted to show where assembly succeeds, and where and why it fails.

---

## 12. Why this project is useful

- Solves a real, foundational problem in bioinformatics: without genome
  assembly, new pathogens could not be characterised and vaccines or treatments
  could not be designed against them.
- Shows how core DSA concepts (hashing, graphs, stacks, DP) sit directly behind
  real scientific breakthroughs.
- Runs on real public data (NCBI genomes), not toy data.
- Produces a measurable result (accuracy against the known original), so
  correctness can be shown objectively.

---

## 13. Scope and limitations

- Handles fragments from both DNA strands (canonical k-mers), real FASTQ
  input, variable read lengths and N letters. Tested on two real SARS-CoV-2
  runs (results in EXPLANATION.md, section 9).
- One typo threshold for the whole genome: with very uneven coverage
  (amplicon data) weak regions are lost (94.4% recovered vs 99.96% on
  capture data).
- No graph clean-up (tip and bubble removal), no quality-score use, no
  adapter or primer trimming.
- Only substitution errors are simulated (no insertions or deletions).
- Repeated regions longer than k cannot be resolved uniquely; the system
  reports contigs in that case, as real assemblers do.
- No filtering of host (human) DNA; fragments come only from the virus.
- Designed for small genomes (viruses, ~10⁴-10⁵ bases), not human-scale
  genomes (~3 × 10⁹ bases).

---

## 14. References

- Idury, R. M. & Waterman, M. S. (1995). A new algorithm for DNA sequence
  assembly. *Journal of Computational Biology*.
- Pevzner, P. A., Tang, H. & Waterman, M. S. (2001). An Eulerian path approach
  to DNA fragment assembly. *PNAS*.
- Zerbino, D. R. & Birney, E. (2008). Velvet: algorithms for de novo short read
  assembly using de Bruijn graphs. *Genome Research*.
- Bankevich, A. et al. (2012). SPAdes: a new genome assembly algorithm.
  *Journal of Computational Biology*.
- Li, D. et al. (2015). MEGAHIT: an ultra-fast single-node solution for large
  and complex metagenomics assembly. *Bioinformatics*.
- Wu, F. et al. (2020). A new coronavirus associated with a human respiratory
  disease in China. *Nature*.
