# How the Genome Assembler Works

A complete walkthrough of [assembler.c](assembler.c): what it does, the data
structures and algorithms behind it (explained from zero), how each part is
coded, what comes out at every stage, and how it performs on both simulated
and **real** SARS-CoV-2 sequencing data.

---

## Contents

0. [Running it](#0-running-it)
1. [The big picture](#1-the-big-picture)
2. [Inputs and outputs](#2-inputs-and-outputs)
3. [Concepts from zero](#3-concepts-from-zero)
4. [The code, stage by stage](#4-the-code-stage-by-stage)
5. [A complete run on the toy example](#5-a-complete-run-on-the-toy-example)
6. [A complete run on simulated SARS-CoV-2 fragments](#6-a-complete-run-on-simulated-sars-cov-2-fragments)
7. [How we judge the result (metrics)](#7-how-we-judge-the-result-metrics)
8. [Experiments on simulated data](#8-experiments-on-simulated-data)
9. [Testing on real sequencing data](#9-testing-on-real-sequencing-data)
10. [Speed (complexity)](#10-speed-complexity)
11. [Simplifications and next steps](#11-simplifications-and-next-steps)
12. [Glossary](#12-glossary)

---

## 0. Running it

Open WSL, go to the project folder, compile once, then run:

```bash
cd /mnt/c/Users/raksh/Downloads/dsa_el
gcc -O2 -o assembler assembler.c

./assembler data/sars_cov2.fasta                          # new random fragments every run
./assembler data/sars_cov2.fasta -seed 42                 # repeat run number 42 exactly
./assembler data/sars_cov2.fasta -cov 30 -k 25            # change settings
./assembler data/sars_cov2.fasta -reads fragments.txt     # assemble a fragments file

# real sequencing data, checked against that sample's own official genome
./assembler data/real/ERR4204399.consensus.fasta -reads data/real/ERR4204399.fastq

# the toy example (tiny k needs -strands 1, see section 5)
./assembler data/toy.fasta -reads data/toy_fragments.txt -k 3 -min 2 -strands 1
```

`gcc -O2` means "compile with optimisations", which makes it run several times
faster. `-o assembler` names the program file.

**Tip for real data:** WSL reads files on the Windows C: drive slowly. Copying
the big `.fastq` files into Linux first makes runs much faster (about 2 s
instead of 6 s):
`mkdir -p ~/dsa && cp data/real/*.fastq ~/dsa/`, then use `-reads ~/dsa/ERR4204399.fastq`.

### Options

| Option | Meaning | Default |
|---|---|---|
| `-k N` | Chunk length (2 to 31) | 31 |
| `-len N` | Fragment length (simulation) | 100 |
| `-cov X` | Coverage: how many times each letter is read, on average (simulation) | 20 |
| `-err X` | Typo rate, 0.01 = 1% of letters misread (simulation) | 0.01 |
| `-min N` | Chunks seen fewer than N times are treated as typos | chosen automatically |
| `-seed N` | Seed for the random numbers (see below) | taken from the clock |
| `-reads FILE` | Use fragments from a file (FASTQ, FASTA or one per line) instead of making new ones | none |
| `-strands N` | 2 = fragments come from both DNA strands, 1 = one strand only | 2 |

### One genome, endless inputs

The genome is the same every time, but **each run cuts it into a different
random set of fragments**. The random choices come from a *seed* number. With no
`-seed`, the program uses the current time, so every run is new. It prints the
seed it used, so any run can be repeated exactly with `-seed`. That's how you
test an assembler: many different random inputs from one known answer. Then,
as a final test, real fragments from a real sequencing machine (section 9).

---

## 1. The big picture

### What the program does

```
 genome FASTA (the answer)
        │
        ▼
 Stage 1  Load genome ─────────────────────────────────────────────┐
        │                                                          │
        ▼                                                          │
 Stage 2  Make fragments (fake sequencing machine)                 │
        │      or load real/saved ones from a file                 │
        │                                                          │
        │   ┌──────────── the assembler only sees this part ──┐    │
        ▼   │                                                 │    │
 Stage 3  Cut fragments into chunks, count them (hash table)  │    │
 Stage 4  Delete rare chunks (typos)                          │    │
 Stage 5  Link chunks into a graph (de Bruijn graph)          │    │
 Stage 6  Walk the graph (Hierholzer) → contigs               │    │
        │   └─────────────────────────────────────────────────┘    │
        │      → contigs.fasta                                     │
        ▼                                                          │
 Stage 7  Compare with the genome ◄────────────────────────────────┘
          (genome fraction, identity, N50, list of differences)
```

The genome is used only at the start (to make fragments) and at the end (to
check the answer). Stages 3–6 never look at it. In real life nobody has the
genome, so if the assembler peeked at it, that would be cheating.

### How the file is organised

[assembler.c](assembler.c) is one file, split into sections that follow the
pipeline from top to bottom:

| Section | What's in it |
|---|---|
| SETTINGS | Default values of k, read length, coverage and so on |
| SMALL HELPERS | `xmalloc` (memory that never fails silently), timing, random numbers |
| DNA LETTERS AND 2-BIT ENCODING | Letters ↔ numbers, and mirror images (reverse complements) |
| HASH TABLE | The lookup table used by stages 3, 5 and 7 |
| STAGE 1 | `read_fasta`: load the genome |
| STAGE 2 | `simulate_reads`, `load_reads`, `save_reads`: get fragments |
| STAGE 3 | `count_kmers`: cut and count chunks |
| STAGE 4 | `typical_count`, `auto_threshold`, `filter_kmers`: delete typo chunks |
| STAGE 5 | `build_graph`: the de Bruijn graph |
| STAGE 6 | `eulerian_walks`, `build_contigs`, `n50`: walk the graph and spell DNA |
| STAGE 7 | `edit_distance`, `place`, `evaluate`: check the answer |
| MAIN | Reads the options, runs stages 1–7 in order, prints everything |

---

## 2. Inputs and outputs

### Input 1: the genome (FASTA file)

Downloaded free from NCBI. SARS-CoV-2 is accession **NC_045512.2**. The file
looks like this:

```
>NC_045512.2 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAA
CGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAAC
...
```

- The first line starts with `>`. It's a label (name of the sequence).
- Every other line is DNA letters, 70 per line. The line breaks mean nothing;
  joined together they form one string of 29,903 letters.
- A genome may contain `N`, meaning "unknown letter here". The official
  genomes of real patient samples have some (section 9).

### Input 2 (optional): fragments

The program accepts three layouts and works out which one it's given:

**Plain**, one fragment per line (what our simulator saves):
```
TCATTCTATTGGATTTGATTACGTCTATAATCCGTTTATGATTGATGTTCAACAATGGGG...
TTTGCAGAGTGGTTTTAAAAAAATGGCATTCCCATCTGGTAAAGTTGAGGGTTGTATGGT...
```

**FASTQ**, the real output format of sequencing machines. Each fragment takes
**4 lines**:
```
@ERR4204399.1 NB552079:30:HNMCFBGXF:1:13107:8648:2619/1      ← name, starts with @
CTCGTGTCCTGTCAACGACAGTAATTAGTTATTAATTATACTGCGTGAGTGCACT...   ← the letters
+                                                            ← separator
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE...   ← quality, one symbol per letter
```
The quality line says how confident the machine is about each letter. We
don't use it; the counting in stage 4 already removes typos.

**FASTA**: `>` label lines are skipped and the other lines are fragments.

### Outputs

| Output | What it is |
|---|---|
| Terminal | Results of every stage (sections 5, 6 and 9) |
| `fragments.txt` | The fragments made in this run (only when simulating) |
| `contigs.fasta` | The rebuilt DNA pieces, longest first, in FASTA format |

---

## 3. Concepts from zero

This section explains every idea the program uses, without assuming any
background. Section 4 then shows how each one is coded.

### 3.1 DNA as text

DNA is a long chain of four chemical "letters": **A**, **C**, **G**, **T**. For
the computer, a genome is simply text: a string of those four characters. In C,
a string is an **array of `char`** (one byte per letter) ending with a special
`'\0'` character that marks the end.

```
genome:  A T T A A A G G T T ...
index:   0 1 2 3 4 5 6 7 8 9 ...
```

### 3.2 Fragments, chunks, coverage

- **Fragment (read):** a short piece copied from a random place in the genome.
  That's what a sequencing machine outputs: 100 letters in our simulator, and
  31–151 letters in the real data.
- **Coverage:** how many fragments, on average, contain each letter. With 20×
  coverage and 100-letter fragments for a 29,903-letter genome:
  number of fragments = 20 × 29,903 / 100 = **5,980**.
- **Chunk (k-mer):** a piece of exactly k letters cut from a fragment. We use
  k = 31. "k-mer" is the textbook name; "mer" means "part".

#### The sliding window

Chunks are taken **at every position, moving one letter at a time**, so they
overlap:

```
Fragment:  A C G T T G
           [A C G]          → ACG
             [C G T]        → CGT
               [G T T]      → GTT
                 [T T G]    → TTG
```

A fragment of length L gives L − k + 1 chunks: 100 − 31 + 1 = **70** chunks
per 100-letter fragment.

**Why neighbouring chunks share k − 1 letters:** each step drops one letter at
the front and adds one at the back. So `ACG` and `CGT` share `CG`. That shared
part is what lets us link chunks later. It isn't a fact about biology; it comes
from how we cut.

#### Why cut fragments into chunks at all?

Fragments overlap at random positions and by random amounts. Checking every
pair of fragments for every possible overlap would take hundreds of millions of
comparisons. Once everything is cut into same-size chunks, overlapping
fragments contain **identical chunks**, and finding identical things is what a
hash table (3.5) does instantly.

### 3.3 Two strands: mirror images

DNA is a **double helix**: two strands wound together. The letters on the two
strands always pair up: **A with T, C with G**. And the two strands run in
**opposite directions**.

```
strand 1, read left to right:   A C G T T G
                                | | | | | |
strand 2, read right to left:   T G C A A C
```

A sequencing machine grabs pieces from **either strand at random**, and always
reads in that strand's own direction. A piece from strand 2 above would come out
as `CAACGT`: strand 2 read right to left.

`CAACGT` is the **reverse complement** of `ACGTTG`: reverse it (`GTTGCA`), then
swap every letter for its partner (G↔C, T↔A), giving `CAACGT`. We'll call it
the **mirror** for short. It's the *same* piece of DNA, just seen from the
other strand.

**Why this matters:** in real data about half the fragments are mirrors. Their
chunks look completely different from the other half, even though they describe
the same DNA. If we ignored this, the graph would split into two unrelated
halves.

**The fix (canonical chunks):** a chunk and its mirror are counted as **one**
chunk. Turn both into numbers (3.4) and always use whichever number is
**smaller** as the shared name. Example with k = 3: `ACG` is 6 and its mirror
`CGT` is 27, so both are stored under 6. This is the standard trick real
assemblers use.

**One catch:** with very short chunks, a genome can contain both a chunk *and*
its mirror at different places, and they get merged by mistake. The toy genome
`ACGTTGCA` contains both `ACG` and `CGT`. That's why the toy example runs with
`-strands 1` (one strand only). With 31-letter chunks this basically never
happens in a virus.

### 3.4 Bits and the 2-bit encoding

#### Bits

Computers store numbers in **binary**, using only 0 and 1. Each 0 or 1 is a
**bit**. For example 6 in binary is `110` (4 + 2 + 0).

A `uint64_t` in C is a whole number made of exactly 64 bits ("unsigned
integer, 64 bits"). It comes from the header `<stdint.h>`.

#### Packing DNA into bits

There are only 4 letters, and 2 bits can store exactly 4 different values:

| Letter | Code | Bits | Partner (mirror) |
|---|---|---|---|
| A | 0 | `00` | T (3) |
| C | 1 | `01` | G (2) |
| G | 2 | `10` | C (1) |
| T | 3 | `11` | A (0) |

So a chunk of 31 letters needs 62 bits and fits inside **one** 64-bit number.
A handy accident: a letter's partner is always **3 − code**.

```
ACG  →  00 01 10  →  binary 000110  →  the number 6
```

**Why bother?** Comparing two numbers is one CPU operation, while comparing
two 31-letter strings means checking up to 31 characters. Numbers are also
easy to feed into a hash function, and they take 8 bytes instead of 32.

#### The bit operators used

| Operator | Name | What it does | Example |
|---|---|---|---|
| `<<` | shift left | Moves all bits left, filling with 0s | `0110 << 2` = `011000` |
| `>>` | shift right | Moves all bits right, dropping the rightmost | `011011 >> 2` = `0110` |
| `\|` | OR | A bit is 1 if it's 1 in either number | `011000 \| 11` = `011011` |
| `&` | AND | A bit is 1 only if it's 1 in both | `1101111 & 111111` = `101111` |
| `^` | XOR | A bit is 1 if it's 1 in exactly one | `110 ^ 001` = `111` |

#### The rolling update

To get the next chunk, we don't re-encode all 31 letters. We update the
previous number:

```
code = ((code << 2) | new_letter) & mask
```

1. `code << 2`: shift left by 2 bits, making room for one new letter.
2. `| new_letter`: put the new letter's 2 bits in the empty space.
3. `& mask`: cut off the oldest letter, which is now too far to the left.
   `mask` is a number whose lowest 2k bits are all 1 (for k = 3: `111111`).

Worked example with k = 3, reading `A C G T T`:

| Letter read | code before | `<< 2` | `\| letter` | `& 111111` | Chunk |
|---|---|---|---|---|---|
| A (00) | `0` | `0` | `00` | `00` | (not 3 letters yet) |
| C (01) | `00` | `0000` | `0001` | `0001` | (not 3 letters yet) |
| G (10) | `0001` | `000100` | `000110` | `000110` | **ACG** (6) |
| T (11) | `000110` | `00011000` | `00011011` | `011011` | **CGT** (27) |
| T (11) | `011011` | `01101100` | `01101111` | `101111` | **GTT** (47) |

In the last row, the `01` (the C) was pushed past the 6th bit and the mask cut
it off. That's the window sliding forward.

#### Rolling the mirror too

The mirror of the window changes in the opposite way: the new letter's
**partner** enters at the **front**, and the oldest letter's partner falls off
the **back**:

```
rev = (rev >> 2) | ((3 - new_letter) << 2(k-1))
```

Reading `A C G T` with k = 3:

| Letter | Partner | rev after | Window | Its mirror |
|---|---|---|---|---|
| A | T (11) | `11 00 00` | (not full yet) | |
| C | G (10) | `10 11 00` | (not full yet) | |
| G | C (01) | `01 10 11` = CGT | ACG | **CGT** ✓ |
| T | A (00) | `00 01 10` = ACG | CGT | **ACG** ✓ |

So both the chunk and its mirror are updated in one step per letter.

Three more tricks used later:
- **First k − 1 letters (prefix):** `code >> 2` drops the last letter.
  CGT = `011011` → `0110` = CG.
- **Last k − 1 letters (suffix):** `code & (mask with 2(k−1) ones)` drops the
  first letter. CGT = `011011` & `1111` → `1011` = GT.
- **Last letter alone:** `code & 3` (3 is `11` in binary). CGT → `11` = T.

### 3.5 The hash table

#### The problem it solves

For real data we have about **20 million** chunks to count. For each one we
need to ask "have I seen this before, and how many times?" Searching a list
from start to end each time would take trillions of steps.

A **hash table** answers "have I seen this?" in about **one step**, no matter
how many items it holds.

#### How it works

Picture a row of numbered boxes called **slots**:

```
slot:   0     1     2     3     4     5     6     7
      [   ] [   ] [   ] [   ] [   ] [   ] [   ] [   ]
```

Each slot can hold one **key** (the chunk, as a number) and one **value**
(here, how many times we've seen it).

A **hash function** turns a key into a slot number. It scrambles the key's bits
so that different keys spread out evenly over the slots. To store or look up
a key, compute its slot and go straight there. There's no searching.

```
hash(ACG) = 5   → look in slot 5
hash(CGT) = 2   → look in slot 2
```

#### Collisions and linear probing

Sometimes two different keys hash to the same slot. That's a **collision**.
Our fix is **linear probing**: if the slot is taken by a *different* key, try
the next slot, then the next, until you find either the key or an empty slot.

```
hash(GTT) = 5, but slot 5 holds ACG  →  try slot 6: empty  →  GTT goes in slot 6
Later, looking up GTT: slot 5 is ACG (not it), slot 6 is GTT (found)
```

This storage style is called **open addressing**: everything lives in one
array and keys move to nearby slots when needed.

#### Keeping it fast: growing when half full

Probing is fast only if most slots are empty, so you hit an empty slot quickly.
The fraction of slots in use is the **load factor**. We keep it at or below
**50%**.

We don't know in advance how many different chunks the input will have (about
140,000 for our simulation, 1.6 million for one real sample). So the table
starts small and **grows**: when it's about to be more than half full, a new
table **twice as big** is made and every key is moved into it. This is called
**rehashing**, because each key's slot is recomputed: the slot depends on the
table size. Each move costs time, but since the size doubles, moves happen
rarely (11 times to grow from 2,048 to 4 million slots). On average, adding a key
still takes about one step.

The size is always a **power of 2** (1,024, 2,048, ... 4,194,304). Then
"slot = hash % size" can be written as "hash & (size − 1)", which is faster.
Example: size 8 = `1000`, size − 1 = `0111`, and `& 0111` keeps the lowest 3
bits, which is a number from 0 to 7.

#### Where we use it

| Stage | Key | Value |
|---|---|---|
| 3 | Chunk (31 letters, canonical) | How many times it was seen |
| 5 | Node (30 letters) | The node's ID number (0, 1, 2, ...) |
| 7 | Chunk of the real genome | Where it is in the genome |

### 3.6 Graphs

A **graph** is a set of **nodes** (dots) joined by **edges** (lines). In a
**directed** graph each edge is an arrow with a direction: `A → B` lets you go
from A to B, but not back.

```
   AC ──► CG ──► GT ──► TT ──► TG ──► GC ──► CA
```

- **Out-degree** of a node: how many arrows leave it.
- **In-degree**: how many arrows enter it.
- A node with out-degree > in-degree is where a path **starts** (more ways out
  than in). In-degree > out-degree means a path **ends** there.
- A node with out-degree > 1 is a **branch**: there are several ways forward.

#### The de Bruijn graph

Our graph is built like this:
- Every **chunk** (k letters) becomes an **edge**.
- The edge goes from the chunk's **first k − 1 letters** to its **last
  k − 1 letters**. Those (k − 1)-letter pieces are the **nodes**.

```
chunk ACG = edge  AC ──► CG
chunk CGT = edge  CG ──► GT
```

Because neighbouring chunks share k − 1 letters, the end node of one chunk is
the start node of the next. So the chunks automatically link up in the same
order they appear in the genome. Walking along the arrows and writing down
letters spells the genome. The graph is named after the mathematician Nicolaas
de Bruijn.

**Both strands:** each kept chunk is added as an edge **and** its mirror is
added as another edge. The graph then contains the genome path *and* its
mirror-image path. When we walk one, we mark the other as used, so each
stretch of DNA is spelled out once.

#### Adjacency list (how a graph is stored)

For each node, we keep a **list of its outgoing edges**. That's called an
**adjacency list** ("adjacent" = next to).

We store each list as a **linked list inside arrays**:
- `first_edge[node]` = number of the node's first outgoing edge (−1 if none)
- `next_edge[edge]` = number of the next edge in the same list (−1 if last)

```
Node CG has two outgoing edges, 4 and 9:

first_edge[CG] = 9     next_edge[9] = 4     next_edge[4] = -1
                 ↓                   ↓                    ↓
              edge 9  ──────────►  edge 4  ──────────►  (end)
```

Adding an edge to a node puts it at the **front** of the list: the new edge
points to the old first edge, then becomes the first edge. It's fast and needs
no resizing.

### 3.7 Eulerian path and Hierholzer's algorithm

#### Eulerian path

An **Eulerian path** is a walk through a graph that uses **every edge exactly
once**. (Leonhard Euler studied this in 1736 with the "bridges of Königsberg"
puzzle.)

Each edge is one chunk of the genome, and the genome contains each chunk once.
So walking every edge once, in the right order, spells the genome.

#### Why "just follow the arrows" isn't enough

Suppose part of the graph looks like this (a loop, which happens when DNA
repeats):

```
A ──► B ──► D
      ▲ │
      │ ▼
      └─C
```

Edges: A→B, B→C, C→B, B→D. The correct walk is **A → B → C → B → D**, which uses
all 4 edges. But if, at B, you happen to pick B→D first, you reach D and get
stuck, with B→C and C→B never used.

#### Hierholzer's algorithm

Hierholzer's algorithm fixes this by keeping a **stack**.

A **stack** is a pile: you can only add to the top (**push**) or remove from
the top (**pop**), like a stack of plates. The last thing pushed is the first
popped.

The algorithm:
1. Push the start node.
2. Look at the node on top of the stack.
   - **If it has an unused edge:** use the edge and push the node it leads to.
   - **If it has no unused edges:** pop it, and add the edge you arrived by to
     the answer.
3. Repeat until the stack is empty.
4. The answer comes out **backwards**, so reverse it.

Running it on the loop example (B's list happens to offer B→D first):

| Step | Stack (top on right) | Action | Answer so far |
|---|---|---|---|
| 1 | A | A has A→B: take it, push B | |
| 2 | A B | B has B→D: take it, push D | |
| 3 | A B D | D has no edges: pop, record B→D | B→D |
| 4 | A B | B still has B→C: take it, push C | B→D |
| 5 | A B C | C has C→B: take it, push B | B→D |
| 6 | A B C B | B has no edges left: pop, record C→B | B→D, C→B |
| 7 | A B C | C has none: pop, record B→C | B→D, C→B, B→C |
| 8 | A B | B has none: pop, record A→B | B→D, C→B, B→C, A→B |
| 9 | A | A has none: pop | |

Reversed: **A→B, B→C, C→B, B→D**. All 4 edges, in a valid order. When the walk
got stuck at D, the stack let it go back to B and fit the forgotten loop in
before D.

#### Why an explicit stack instead of recursion

Hierholzer's algorithm is often written with **recursion** (a function calling
itself). For our ~60,000 edges, that would mean up to ~60,000 nested calls,
which can overflow the program's call stack and crash it. Our own array-based
stack has no such limit.

#### When there isn't one perfect path: contigs

In practice the graph is broken into pieces. A missing chunk (seen too few
times) cuts the chain, and leftover typos add side branches. Then no single
path covers everything. Hierholzer still uses every edge once, but some
consecutive edges in its output won't actually connect. We cut the output at
those points. Each connected piece is a **contig** (from "contiguous"): a
stretch of genome we're confident about. Real assemblers also output contigs.

### 3.8 Dynamic programming and edit distance

#### Edit distance

The **edit distance** between two strings is the smallest number of
single-letter changes that turn one into the other. A change is one of:
- **substitute** a letter (`ACGT` → `AGGT`)
- **delete** a letter (`ACGT` → `AGT`)
- **insert** a letter (`AGT` → `ACGT`)

We use it to count how many letters in our rebuilt DNA differ from the real
genome. An `N` (unknown letter) in the genome matches anything, because nobody
knows what the true letter is.

#### Dynamic programming

**Dynamic programming (DP)** solves a big problem by solving small versions
first, storing their answers in a table, and building on them.

For edit distance: let `D[i][j]` = edit distance between the first `i` letters
of string a and the first `j` letters of string b. Then:

```
D[i][j] = smallest of:
    D[i-1][j-1] + (0 if a[i] == b[j], else 1)    substitute (or match)
    D[i-1][j]   + 1                              delete a[i]
    D[i][j-1]   + 1                              insert b[j]
```

The first row and column are easy: turning `i` letters into nothing takes `i`
deletions.

Example: a = `ACGT`, b = `AGT`:

|   | "" | A | G | T |
|---|---|---|---|---|
| **""** | 0 | 1 | 2 | 3 |
| **A** | 1 | 0 | 1 | 2 |
| **C** | 2 | 1 | 1 | 2 |
| **G** | 3 | 2 | 1 | 2 |
| **T** | 4 | 3 | 2 | **1** |

The bottom-right cell is the answer: **1** (delete the C).

#### Saving memory: two rows

Each row only needs the row above it. So instead of the whole table
(30,000 × 30,000 = 900 million cells, about 3.6 GB), we keep just **two rows**
(about 240 KB) and swap them after each row.

### 3.9 Sorting and N50

**Sorting** puts items in order. We sort contigs longest first using C's
built-in `qsort` (quicksort), which needs a small function that says which of
two items comes first.

**N50** describes how long the contigs typically are. Sort contigs longest
first and add up their lengths. N50 is the length of the contig at which the
running total reaches **half** of the total length.

Example: contigs 10, 8, 5, 3, 2 (total 28, half is 14).
10 → total 10 (not yet). 10 + 8 → 18 (≥ 14). **N50 = 8.**

It means "half of the assembled DNA is in contigs at least this long". Higher
is better, and a perfect assembly has N50 = genome length.

### 3.10 Random numbers and seeds

`rand()` gives "random" numbers, but they come from a formula that starts from a
**seed**. The same seed always produces the same sequence of numbers. So:
- a different seed gives different fragments (a new test input), and
- the same seed gives the exact same fragments (a repeatable experiment).

`srand(seed)` sets the seed. `rand() % n` gives a whole number from 0 to n−1.

---

## 4. The code, stage by stage

For each stage: **why** it exists, its **input → output**, and **how** it's
coded.

### 4.0 Settings and helpers

```c
static int      K         = 31;
static int      READ_LEN  = 100;
static double   COVERAGE  = 20;
static double   ERR_RATE  = 0.01;
static int      MIN_COUNT = 0;     /* 0 = choose automatically */
static unsigned SEED      = 0;
static int      STRANDS   = 2;
```

These are **global variables**: every function can read them. `static` means
"visible only inside this file". `main` overwrites them if you pass options.

```c
static void *xmalloc(size_t n)
{
    void *p = malloc(n);
    if (!p) { fprintf(stderr, "Out of memory\n"); exit(1); }
    return p;
}
```

`malloc(n)` asks the operating system for n bytes of memory and returns a
**pointer** (its address), or `NULL` if there isn't enough. `xmalloc` wraps it
so the program stops with a clear message instead of crashing later.
`xcalloc` does the same but also fills the memory with zeros, and `xrealloc`
grows an existing block.

The letter helpers:
- `is_base(c)` tells whether c is A/C/G/T.
- `base_code(c)` turns a letter into 0–3.
- `complement(c)` gives the partner letter.
- `decode(code, len, out)` turns a number back into letters.
- `revcomp_code(code, len)` and `revcomp_string(s, len)` compute mirrors.

```c
static uint64_t revcomp_code(uint64_t code, int len)
{
    uint64_t rc = 0;
    for (int i = 0; i < len; i++) {
        rc = (rc << 2) | (3 - (code & 3));   /* partner of the last letter */
        code >>= 2;                          /* move to the previous letter */
    }
    return rc;
}
```

This takes letters off the **end** of the chunk one by one (`code & 3`), and
appends each one's partner (`3 − letter`) to the result. Reading backwards
while swapping partners is exactly the reverse complement.

### 4.1 Stage 1: load the genome (`read_fasta`)

**Why:** we need the real genome to make fragments and to check the answer.

**Input → output:** file name → one `char` array (29,903 letters for
SARS-CoV-2).

```c
while (fgets(line, sizeof line, f)) {
    if (line[0] == '>') {
        if (++records > 1) break;   /* only the first sequence */
        continue;                   /* skip the label line */
    }
    for (char *p = line; *p; p++) {
        char c = (char)toupper((unsigned char)*p);
        if (!isalpha((unsigned char)c)) continue;           /* skip newlines etc. */
        if (len + 1 >= cap) { cap *= 2; seq = xrealloc(seq, cap); }
        seq[len++] = is_base(c) ? c : 'N';                  /* unknown → N */
    }
}
```

- `fgets` reads one line of the file at a time.
- Label lines (starting with `>`) are skipped.
- `toupper` turns lowercase into uppercase. Newlines and spaces are skipped.
  Any letter that isn't A/C/G/T becomes `N`. `N`s are **kept**, not removed,
  so every letter stays at its correct position.
- **Dynamic array:** we don't know the genome's length in advance. We start
  with room for 1,024 letters, and whenever it's full we **double** the space
  with `realloc`. Doubling means we only grow about 15 times for 30,000
  letters, instead of once per letter.

### 4.2 Stage 2: get fragments

#### Making them (`simulate_reads`)

**Why:** we don't have a sequencing machine, so we imitate one: random pieces,
occasional typos, and about half read from the other strand.

**Input → output:** genome → 5,980 strings of 100 letters (an array of
pointers, `char **`), also saved to `fragments.txt`.

```c
int n = (int)(COVERAGE * g_len / READ_LEN);           /* 20 × 29903 / 100 = 5980 */
for (int i = 0; i < n; i++) {
    int start = rand() % (g_len - READ_LEN + 1);      /* random start position */
    char *r = xmalloc(READ_LEN + 1);
    memcpy(r, genome + start, READ_LEN);              /* copy 100 letters */
    r[READ_LEN] = '\0';

    for (int j = 0; j < READ_LEN; j++) {
        if (rand01() < ERR_RATE) {                    /* 1% chance per letter */
            r[j] = BASES[(base_code(r[j]) + 1 + rand() % 3) % 4];
            errors++;
        }
    }
    if (STRANDS == 2 && rand() % 2) {                 /* coin flip: other strand */
        revcomp_string(r, READ_LEN);
        flipped++;
    }
    reads[i] = r;
}
```

- `rand() % (g_len - READ_LEN + 1)` picks a start from 0 to 29,803, so the
  fragment never runs past the end of the genome.
- `genome + start` is **pointer arithmetic**: the address of letter number
  `start`. `memcpy` copies 100 bytes from there.
- **Typo:** `base_code(r[j])` is the true letter's code (0–3). Adding 1, 2 or 3
  and taking `% 4` always gives a *different* letter. For example, for G (2),
  adding 1 gives 3 (T), adding 2 gives 0 (A), and adding 3 gives 1 (C).
- **Other strand:** `rand() % 2` is 0 or 1, a coin flip. On 1, the fragment is
  replaced by its mirror, just like a real machine reading the other strand.

#### Loading them (`load_reads`)

**Input → output:** a file → array of fragments, any lengths.

```c
while (fgets(line, sizeof line, f)) {
    if (fastq == -1) fastq = (line[0] == '@');    /* decide from the first line */
    line_no++;
    if (fastq ? line_no % 4 != 2 : line[0] == '>') continue;
    ...
        if (isalpha((unsigned char)c)) r[len++] = is_base(c) ? c : 'N';
    ...
}
```

- If the first line starts with `@`, the file is FASTQ. Then only line 2 of
  every group of 4 holds letters, so we keep lines where `line_no % 4 == 2`.
  `%` is the remainder: lines 2, 6, 10, ... give remainder 2.
- Otherwise, lines starting with `>` are skipped and every other line is a
  fragment.
- `cond ? a : b` means "if cond then a else b".
- Fragments can have any length; each one is stored with its own `malloc`.

### 4.3 Stage 3: cut into chunks and count (`count_kmers`)

**Why:** chunks turn messy random overlaps into identical pieces, and counting
them is how we find typos in stage 4.

**Input → output:** fragments → hash table of (canonical chunk → count).

```c
HashTable counts = ht_create(1024);          /* starts small, grows as needed */
uint64_t mask = (1ULL << (2 * K)) - 1;      /* 62 ones */
int shift = 2 * (K - 1);                    /* where the first letter sits */

for (int i = 0; i < n_reads; i++) {
    uint64_t fwd = 0, rev = 0;
    int run = 0;                             /* valid letters in a row */
    for (int j = 0; reads[i][j]; j++) {
        char c = reads[i][j];
        if (!is_base(c)) { run = 0; continue; }   /* N: start again after it */
        int b = base_code(c);
        fwd = ((fwd << 2) | b) & mask;                       /* the chunk */
        rev = (rev >> 2) | ((uint64_t)(3 - b) << shift);     /* its mirror */
        if (++run >= K) {
            uint64_t canon = (STRANDS == 2 && rev < fwd) ? rev : fwd;
            (*ht_get_or_add(&counts, canon, &is_new))++;     /* count it */
        }
    }
}
```

- `1ULL << 62` is 1 followed by 62 zeros. Subtracting 1 gives 62 ones, which
  is the mask. `ULL` means "this constant is a 64-bit number"; without it, the
  shift would overflow a 32-bit number.
- `fwd` and `rev` are the rolling chunk and rolling mirror from 3.4.
- **`run`** counts valid letters since the start or since the last `N`. A
  chunk is only counted once there are K valid letters in a row, so no chunk
  ever contains an `N` (we don't know what letter that is).
- **`canon`**: the smaller of chunk and mirror, so both strands count
  together (3.3).
- `ht_get_or_add` returns a **pointer to the count** in the table (adding the
  chunk with count 0 if it's new). `(*pointer)++` adds 1 to the count right
  there in the table.

#### The hash table code

```c
typedef struct {
    uint64_t      *keys;   /* the chunks             */
    uint32_t      *vals;   /* their counts / IDs     */
    unsigned char *used;   /* 1 if the slot is taken */
    uint64_t       cap;    /* number of slots        */
    uint64_t       size;   /* number of keys stored  */
} HashTable;
```

A `struct` groups related variables under one name. The table is three
parallel arrays: slot `i` holds `keys[i]`, `vals[i]` and `used[i]`.

```c
static uint64_t ht_slot(const HashTable *t, uint64_t key)
{
    uint64_t h = key * 0x9E3779B97F4A7C15ULL;   /* scramble the bits */
    h ^= h >> 32;
    uint64_t i = h & (t->cap - 1);              /* same as h % cap */
    while (t->used[i] && t->keys[i] != key)
        i = (i + 1) & (t->cap - 1);             /* linear probing */
    return i;
}
```

- **Hash function:** multiplying by a big odd constant mixes all the bits
  together (this constant is a common choice, derived from the golden ratio).
  `h ^= h >> 32` mixes the high bits into the low bits, because the next line
  only uses the low bits.
- **Probing:** while the slot is taken by a different key, move to the next
  one. `& (cap - 1)` wraps around from the last slot to slot 0.
- It returns either the slot holding `key` or the empty slot where it belongs.

```c
static void ht_grow(HashTable *t)
{
    HashTable bigger = ht_create(t->cap);        /* twice as many slots */
    for (uint64_t i = 0; i < t->cap; i++) {
        if (!t->used[i]) continue;
        uint64_t j = ht_slot(&bigger, t->keys[i]);   /* new slot in the big table */
        bigger.used[j] = 1;
        bigger.keys[j] = t->keys[i];
        bigger.vals[j] = t->vals[i];
        bigger.size++;
    }
    ht_free(t);                                  /* give back the old arrays */
    *t = bigger;
}

static uint32_t *ht_get_or_add(HashTable *t, uint64_t key, int *is_new)
{
    if (2 * (t->size + 1) > t->cap) ht_grow(t);  /* keep it at most half full */
    uint64_t i = ht_slot(t, key);
    ...
}
```

`ht_create(n)` makes the smallest power-of-2 table with at least 2n slots, so
`ht_create(t->cap)` gives exactly twice the current size. Before adding a key,
`ht_get_or_add` checks whether the table would become more than half full and
grows it first. For the real capture sample the table ends at 4,194,304 slots
(about 55 MB).

### 4.4 Stage 4: delete typo chunks

**Why:** a real chunk appears in many fragments because of coverage. A chunk
containing a typo almost always appears far fewer times, because the same random
typo at the same place rarely happens again. If typo chunks stayed, they'd add
false branches to the graph.

**Input → output:** hash table of counts → plain array of the chunks seen at
least *threshold* times.

#### Choosing the threshold automatically

The counts form two groups, which you can see in the histogram (section 6):
a huge pile of typo chunks at low counts, and a hump of real chunks at high
counts. The threshold should sit between them. The program finds it in two
steps:

```c
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
```

**1. Typical count of a real chunk.** `hist[c]` is how many different chunks
were seen exactly c times. A chunk seen c times accounts for c **sightings**.
Walk up the counts, adding sightings, until half of all sightings are reached.
That count is the **median count per sighting**. Real chunks are seen so often
that they make up most sightings, so this lands on a typical real chunk's
count. Chunks seen once are left out first, because with many typos they could
otherwise outnumber everything.

```c
static int auto_threshold(const long *hist)
{
    int typical = typical_count(hist);
    int dip = 2;
    for (int c = 1; c < 2 * typical && c < HIST_MAX; c++)
        if (hist[c + 1] > hist[c]) { dip = c; break; }
    int tenth = typical / 10;
    return dip > tenth ? dip : tenth;
}
```

**2. The threshold** is the larger of two numbers:
- **The dip:** the first count where the histogram goes back *up*. That's
  where the typo pile ends and the real hump begins. We only look below twice
  the typical count; further out, the numbers are small and bump up and down
  by chance.
- **A tenth of the typical count:** in very deep real data (~800×), the same
  typo or machine glitch can turn up dozens of times. Anything seen less than a
  tenth as often as a typical real chunk is treated as junk.

| Data | Typical count | Dip | A tenth | Threshold |
|---|---|---|---|---|
| Simulated, 20× | 11 | 3 | 1 | **3** |
| Simulated, 5× | low | none found (2) | 0 | **2** |
| Real capture sample, 832× | 637 | (noise) | 63 | **63** |
| Real amplicon sample, 790× | 794 | 9 | 79 | **79** |

You can always override it with `-min N`.

#### Filtering

```c
for (uint64_t i = 0; i < counts->cap; i++)
    if (counts->used[i] && counts->vals[i] >= (uint32_t)threshold)
        solid[n++] = counts->keys[i];
```

This walks through every slot and keeps the frequent chunks. The kept chunks
are often called "solid" k-mers.

### 4.5 Stage 5: build the graph (`build_graph`)

**Why:** the graph records which chunk can follow which, so that walking it
spells the genome.

**Input → output:** array of chunks → `Graph` (nodes, edges, adjacency lists,
degrees).

```c
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
```

Edge number `e` is described by `edge_kmer[e]`, `edge_from[e]` and
`edge_to[e]`. Nodes are numbered 0, 1, 2, ... so they can be used as array
indexes.

```c
for (int e = 0; e < n_edges; e++) {                  /* n_edges = 2 × chunks */
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
```

- **Both strands:** with `STRANDS = 2`, edges come in pairs: edge 0 is chunk 0,
  edge 1 is chunk 0's mirror, edge 2 is chunk 1, edge 3 is chunk 1's mirror,
  and so on. So `e / 2` is the chunk number and `e % 2` says whether this edge
  is the mirror. It also means **the mirror of edge e is always `e ^ 1`**:
  XOR with 1 flips the last bit, turning 6 into 7 and 7 into 6.
- `kmer >> 2` is the prefix and `kmer & suffix_mask` is the suffix (see 3.4).
- `node_id` uses a **second hash table** (30-letter piece → ID). The first time
  a piece is seen it gets the next free number. After that, the same piece
  always gets the same number. That's how the end of one chunk and the start of
  the next become the *same* node, which links them.
- The two `first_edge`/`next_edge` lines add edge `e` to the front of node
  `from`'s list (3.6).
- Degrees are counted as edges are added. Stage 6 uses them to find where the
  genome starts.

### 4.6 Stage 6: walk the graph (`eulerian_walks`)

**Why:** a walk that uses every edge once lists the chunks in genome order.

**Input → output:** graph → array `order` of edge numbers in walking order,
with −1 marking the start of each separate walk.

```c
static int next_unused(Graph *g, int v)
{
    int e = g->first_edge[v];
    while (e != -1 && g->used[e]) e = g->next_edge[e];   /* skip used edges */
    g->first_edge[v] = e;
    return e;
}
```

`next_unused` finds a node's first edge that hasn't been used yet. Any used
edges at the front of the list are dropped from it, so they're never looked at
again.

```c
for (int pass = 0; pass < 2; pass++) {
    for (int s = 0; s < g->n_nodes; s++) {
        if (next_unused(g, s) == -1) continue;                   /* nothing left here */
        if (pass == 0 && g->out_deg[s] <= g->in_deg[s]) continue;

        int top = 0, n_walk = 0;
        stack_node[0] = s;
        stack_edge[0] = -1;

        while (top >= 0) {
            int v = stack_node[top];
            int e = next_unused(g, v);
            if (e != -1) {                          /* unused edge: take it */
                g->used[e] = 1;
                if (STRANDS == 2) g->used[e ^ 1] = 1;   /* its mirror is now done too */
                g->first_edge[v] = g->next_edge[e];
                top++;
                stack_node[top] = g->edge_to[e];
                stack_edge[top] = e;
            } else {                                /* stuck: step back */
                if (stack_edge[top] != -1) walk[n_walk++] = stack_edge[top];
                top--;
            }
        }
        order[n_order++] = -1;
        for (int i = n_walk - 1; i >= 0; i--) order[n_order++] = walk[i];
    }
}
```

- **Choosing where to start:** pass 0 starts walks only at start-like nodes
  (out > in), which is where the genome, or a piece of it, begins. Pass 1 picks
  up any edges still unused, for example a piece that forms a loop.
- **The stack** is two arrays: `stack_node[top]` is the node on top, and
  `stack_edge[top]` is the edge we used to get there. `top` is the index of the
  top item; `top++` pushes and `top--` pops.
- **Using up an edge** marks it used, and also marks its **mirror** (`e ^ 1`)
  used. Once we've walked a stretch of DNA on one strand, the same stretch on
  the other strand is skipped, so it's not spelled out twice.
- **Recording:** when we pop, we save the edge we arrived by. These come out in
  reverse order, so they're copied into `order` backwards.

### 4.7 Stage 6b: turn edges into DNA (`build_contigs`)

**Why:** the walk is a list of edge numbers; we need letters.

**Input → output:** `order` → array of contigs (strings), sorted longest first
and saved to `contigs.fasta`.

```c
int j = i + 1;
while (j < n_order && order[j] != -1 &&
       g->edge_from[order[j]] == g->edge_to[order[j - 1]])
    j++;                                       /* extend while edges connect */

int len = K + (j - i - 1);
decode(g->edge_kmer[order[i]], K, s);          /* first chunk: all 31 letters */
for (int t = i + 1; t < j; t++)
    s[K + (t - i - 1)] = BASES[g->edge_kmer[order[t]] & 3];   /* then 1 letter each */
```

- A run of edges continues as long as each edge starts where the previous one
  ended. A walk boundary (−1) or a non-connecting pair ends the contig.
- **Spelling:** the first chunk gives 31 letters. Each following chunk adds only
  its **last** letter (`& 3`), because its other 30 letters overlap the
  previous chunk. A contig of m edges has 31 + (m − 1) letters.
- A contig may come out on either strand (as the genome or as its mirror).
  Stage 7 handles that.

Then:

```c
qsort(contigs, n_contigs, sizeof(Contig), by_length_desc);
```

`qsort` sorts any array. You give it the array, how many items, the size of
one item, and a **comparison function**. `by_length_desc` returns a positive
number when `b` is longer than `a`, meaning "b goes first", which gives
longest-first order.

`n50` then walks the sorted list exactly as in 3.9.

### 4.8 Stage 7: check the answer (`evaluate`)

**Why:** to measure, honestly and with numbers, how well the assembly worked.
This is the only stage after stage 2 that uses the genome.

**Input → output:** contigs + genome → printed metrics (section 7).

Steps:
1. **Index the genome:** a third hash table, (chunk of the genome → its
   position). Chunks containing `N` are skipped.
2. **Place each contig** (`place`): find the first chunk of the contig that
   exists in the genome. If that chunk is at genome position `pos` and was at
   offset `i` inside the contig, the contig starts at genome position
   `pos − i`. If no chunk matches, the contig might be on the other strand, so
   its mirror is tried. If that fails too, it's counted as "not found" (junk).
3. **Compare:** `edit_distance(contig, genome + start, length)` counts the
   letters that differ. An `N` in the genome matches anything.
4. **List differences:** wherever the contig and genome disagree letter by
   letter, the position and both letters are recorded, for example `241:C>T`
   (genome has C at position 241, our contig has T). In real data these are
   often real mutations (section 9).
5. **Mark coverage:** set `covered[p] = 1` for every genome position a contig
   sits on. Positions still 0 at the end (and not `N`) are missing stretches.

```c
static int edit_distance(const char *a, int n, const char *b, int m)
{
    for (int j = 0; j <= m; j++) prev[j] = j;           /* first row: 0,1,2,... */
    for (int i = 1; i <= n; i++) {
        cur[0] = i;                                     /* first column */
        for (int j = 1; j <= m; j++) {
            int sub = prev[j - 1] + (a[i - 1] != b[j - 1] && b[j - 1] != 'N');
            int del = prev[j] + 1;
            int ins = cur[j - 1] + 1;
            int best = sub < del ? sub : del;
            cur[j] = best < ins ? best : ins;
        }
        int *tmp = prev; prev = cur; cur = tmp;         /* swap the two rows */
    }
    return prev[m];
}
```

This is the formula from 3.8, with `prev` as the row above and `cur` as the row
being filled. The substitution cost is 1 only if the letters differ *and* the
genome letter isn't `N`. Swapping the two pointers reuses the memory without
copying.

### 4.9 `main`

1. **Read the options:** loops over `argv` (the words typed after the program
   name). For each option like `-k`, it reads the next word and converts it
   with `atoi` (text → int) or `atof` (text → decimal).
2. **Set the seed:** `srand(SEED)`, with the clock if no seed was given.
3. **Run stages 1–7 in order**, printing each stage's results and time
   (`clock()` measures CPU time).
4. **Free all memory:** every `malloc` has a matching `free`, so nothing leaks.

---

## 5. A complete run on the toy example

This is the example from the earlier chat: the genome `ACGTTGCA` and seven
5-letter fragments, one with a typo (`GTAGC` should be `GTTGC`).

```bash
./assembler data/toy.fasta -reads data/toy_fragments.txt -k 3 -min 2 -strands 1
```

`-strands 1` is needed because the toy genome contains both `ACG` and its
mirror `CGT` (3.3). With mirror handling on, they would be merged into one
chunk. `-min 2` sets the threshold by hand; with only 7 fragments there's no
meaningful histogram.

`data/toy_fragments.txt`:
```
ACGTT
ACGTT
CGTTG
GTTGC
TTGCA
TTGCA
GTAGC
```

**Stage 1**, input `toy.fasta` → output `ACGTTGCA` (8 letters).

**Stage 2**, input the file → output 7 fragments.

**Stage 3**: each 5-letter fragment gives 3 chunks (5 − 3 + 1):

```
ACGTT → ACG CGT GTT        (×2, it appears twice)
CGTTG → CGT GTT TTG
GTTGC → GTT TTG TGC
TTGCA → TTG TGC GCA        (×2)
GTAGC → GTA TAG AGC        ← all three contain the typo
```

| Chunk | ACG | CGT | GTT | TTG | TGC | GCA | GTA | TAG | AGC |
|---|---|---|---|---|---|---|---|---|---|
| Count | 2 | 3 | 4 | 4 | 3 | 2 | 1 | 1 | 1 |

Program output: `Different chunks : 9`, with 3 seen once, 2 seen twice, 2 seen
3 times and 2 seen 4 times. That matches the table.

**Stage 4** (threshold 2): `GTA`, `TAG`, `AGC` are removed. Program output:
`Kept 6, Removed 3`.

**Stage 5**: 6 edges, 7 nodes:

```
AC ─ACG─► CG ─CGT─► GT ─GTT─► TT ─TTG─► TG ─TGC─► GC ─GCA─► CA
```

Program output: `Nodes 7, Edges 6, Start-like 1 (AC), End-like 1 (CA),
Branching 0`. If `GTA` had been kept, node GT would have two ways out (GTT and
GTA): a false branch.

**Stage 6**: start at AC (out 1 > in 0). The walk goes all the way to CA,
gets stuck, and pops everything:

```
popped:   GCA, TGC, TTG, GTT, CGT, ACG
reversed: ACG, CGT, GTT, TTG, TGC, GCA
spelled:  ACG + T + T + G + C + A = ACGTTGCA
```

Program output: `Contigs 1, Longest 8, Contig 1: ACGTTGCA`.

**Stage 7**: `Genome recovered 100.00%`, `Letters matching 100.0000%`.

Full output: [results/toy.txt](results/toy.txt).

---

## 6. A complete run on simulated SARS-CoV-2 fragments

```bash
./assembler data/sars_cov2.fasta -seed 42
```

Defaults: k = 31, fragment length 100, coverage 20×, 1% typos, both strands,
automatic threshold. Below is the real output, with what each number means.
Full output: [results/simulated_seed42.txt](results/simulated_seed42.txt).

### Stage 1: load genome
```
  Length        : 29903 letters
  Starts with   : ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT...
```
The full SARS-CoV-2 reference genome.

### Stage 2: make fragments
```
  Typos added   : 5963
  Other strand  : 2992 fragments read from the other strand
  Fragments     : 5980  (lengths 100-100, 0 contain N)
  Letters       : 598000  (about 20x coverage)
  Fragment 1    : ACAATACAGATCATGGTTGCTTTGTAGGTTACCTGTAAAACCCCATTGTTGAACATCAAT...
```
5,980 = 20 × 29,903 / 100. About 1% of the 598,000 letters are typos (5,963),
and about half the fragments (2,992) are mirrors.

### Stage 3: cut and count
```
  Different chunks : 140570  (a chunk and its mirror count as one)
  How many chunks were seen N times:
    N =  1  : 108790      ← almost all typo chunks
    N =  2  : 1955
    N =  3  : 276         ← the dip: our threshold
    N =  4  : 566
    N =  5  : 1087
    N =  6  : 1846
    N =  7  : 2480
    N =  8  : 3122
    N =  9  : 3599        ← the hump of real chunks
    N = 10  : 3651
    N = 11+ : 13198
```
The genome has only 29,873 real chunks (29,903 − 31 + 1), but we saw 140,570
different ones. The extra ~110,000 are typo chunks: each typo spoils up to 31
chunks, and each spoiled chunk is unique. Expected count of a real chunk:
20 × 70/100 × 0.99³¹ ≈ 10 (coverage × fraction of positions where a chunk can
start × chance the chunk has no typo).

### Stage 4: filter
```
  Typical count of a real chunk : 11
  Threshold     : 3  (automatic: ...)
  Kept          : 29825
  Removed       : 110745
```
The program found the dip at 3 by itself. 29,825 kept is very close to the
29,873 real chunks.

### Stage 5: graph
```
  Nodes (30-letter pieces) : 59686
  Edges (chunks x 2)       : 59650
  Start-like nodes (out > in) : 40
  End-like nodes   (in > out) : 40
  Branching nodes  (out > 1)  : 2
```
Edges are twice the chunks because each chunk is added with its mirror. A
perfect graph would have 2 starts and 2 ends (one per strand). Each lost chunk
cuts the chain, adding a start and an end on both strands. 40 starts means
about 19 cuts. The 2 branches come from typo chunks that happened to survive.

### Stage 6: walk and spell
```
  Contigs       : 21  (saved to contigs.fasta)
  Longest       : 6672 letters
  Total length  : 30455 letters
  N50           : 2619
```
The genome came back in 21 pieces. The total (30,455) is a bit more than the
genome because neighbouring contigs overlap by up to 29 letters at each cut.
One contig starts with `TTTTTTTTTT...`: it's the genome's end (which finishes
with 33 A's) spelled from the other strand.

### Stage 7: check
```
  Genome recovered (genome fraction) : 99.94%  (29884 of 29903 letters)
  Letters matching (identity)        : 99.9967%  (1 differences in 30380 letters)
  Contigs on the other strand        : 11  (flipped before comparing)
  Contigs not found in genome        : 2
  Single-letter differences          : 1 at genome positions:
    1452:T>G
  Missing stretch                    : positions 1-17 (17 letters)
  Missing stretch                    : positions 29902-29903 (2 letters)
```
- **99.94%** of the genome is present in our contigs.
- Of the 30,380 letters we placed, **only 1** is wrong (position 1,452).
- 11 of the 21 contigs came out as mirrors and were flipped before comparing.
- 2 small contigs are made of leftover typos.
- **Only the very ends are missing.** Fragments start at random positions but
  can't start before letter 1, so the first few letters are covered by very
  few fragments. The first chunk can only come from fragments starting exactly
  at position 1, which is about 0.2 fragments on average, well below the
  threshold. Real sequencing has the same problem at genome ends.

The whole run takes **under half a second**.

---

## 7. How we judge the result (metrics)

A single "accuracy %" can't describe an assembly, because it can fail in
different ways. It can miss parts, get letters wrong, come back in many small
pieces, or include junk. Each metric measures one of those. These are the same
measures used by QUAST, the standard tool for evaluating genome assemblies.

| Metric | Question it answers | Perfect value |
|---|---|---|
| **Genome fraction** | How much of the genome did we recover? | 100% |
| **Identity** | Of the letters we output, how many match the genome? | 100% |
| **Number of contigs** | How many pieces did the genome come back in? | 1 |
| **Longest contig** | How big is the biggest piece? | Genome length |
| **N50** | How long are the pieces typically? | Genome length |
| **Contigs not found** | How much junk did we output? | 0 |
| **Single-letter differences** | Exactly where do we disagree with the genome? | None |
| **Missing stretches** | Exactly which parts are missing? | None |

Why both genome fraction and identity: you could score 100% identity by
outputting one short perfect contig (but genome fraction would be tiny). You
could get high genome fraction with error-filled contigs (but identity would be
low). A good assembly needs both.

How they're computed:
- **Genome fraction** = known genome positions covered by a placed contig ÷
  known genome positions (`N` positions are left out)
- **Identity** = (placed letters − total edit distance) ÷ placed letters

---

## 8. Experiments on simulated data

All runs use seed 42, so they're repeatable. Defaults: k = 31, length 100,
coverage 20, typos 1%, both strands, automatic threshold. Each row changes one
or two settings. Data: [results/simulated_experiments.txt](results/simulated_experiments.txt).

| Run | Threshold | Contigs | Longest | N50 | Genome fraction | Identity | Junk |
|---|---|---|---|---|---|---|---|
| Defaults | 3 | 21 | 6,672 | 2,619 | 99.94% | 99.997% | 2 |
| No typos (`-err 0`) | 2 | **1** | 29,885 | 29,885 | 99.94% | 100% | 0 |
| Coverage 30 | 4 | 3 | 15,766 | 15,766 | 99.92% | 100% | 0 |
| Coverage 40 | 4 | **1** | 29,883 | 29,883 | 99.93% | 100% | 0 |
| Coverage 5 | 2 | 353 | 511 | 119 | 91.44% | 99.991% | 4 |
| k = 21 | 3 | 8 | 16,486 | 16,486 | 99.94% | 99.993% | 1 |
| k = 15 | 4 | 13 | 11,749 | 5,481 | **71.59%** | **76.91%** | 0 |
| Threshold 2 (`-min 2`) | 2 | 150 | 6,222 | 1,873 | 99.98% | 99.932% | 128 |
| Typos 3% | 3 | 426 | 785 | 152 | 99.42% | 99.995% | 12 |
| Typos 3%, coverage 40 | 4 | 45 | 3,070 | 1,260 | 99.75% | 100% | 1 |

Five different seeds with default settings gave 9, 12, 12, 14 and 19 contigs,
with genome fraction 99.83–99.97% and identity 100% every time.

What they show:

- **No typos, or enough coverage (40×):** the whole genome comes back as **one
  single contig**. The algorithm works perfectly when the data is good enough.
- **Low coverage (5×):** real chunks are seen only ~2–3 times, so many drop out.
  The genome falls apart into 353 pieces and 9% of it is lost. **Enough coverage
  is essential.**
- **Small k (15):** **this is why we use k = 31.** With 15-letter chunks, many
  14-letter pieces appear in more than one place in the genome (and some appear
  as each other's mirrors). Their nodes merge, the graph tangles, and the walk
  joins parts of the genome that don't belong together. Identity drops to 77%.
  With k = 21 or 31, almost every piece is unique.
- **Threshold too low (2):** typo chunks seen twice survive, which adds false
  branches. The result is 150 pieces and 128 junk contigs. The automatic
  threshold avoids this.
- **More typos (3%):** only 0.97³¹ ≈ 39% of chunks are typo-free, so real chunks
  are seen fewer times and many drop out. Extra coverage makes up for it.
- **Randomness:** results vary from seed to seed. In the report, run each
  setting with several seeds and report the average.

---

## 9. Testing on real sequencing data

Simulated data is the standard first test, because you know the true answer.
But our simulator only creates problems our assembler already handles. The
real test is data from an actual sequencing machine.

### 9.1 The data

Real sequencing runs are stored for free in NCBI's **SRA** (Sequence Read
Archive), mirrored by Europe's **ENA**. We picked two runs from **March–April
2020** in the UK, from the COG-UK project:

| | Sample 1 | Sample 2 |
|---|---|---|
| Run | ERR4204399 | ERR4239992 |
| Collected | 4 April 2020 | 14 March 2020 |
| Method | **Targeted capture**: probes fish virus DNA out of the sample | **Amplicon**: the genome is copied in ~100 overlapping pieces by PCR |
| Machine | Illumina NextSeq 550 | Illumina MiSeq |
| Fragments | 167,610 | 193,018 |
| Lengths | 54–151 letters | 31–151 letters |
| Coverage | ~832× | ~790× |

**Why early 2020:** the virus had only a handful of mutations compared to the
reference genome then, so differences are easy to interpret.

**The true answer:** for each sample, COG-UK published an official **consensus
genome**: the sample's actual sequence, worked out by the UK's standard
pipeline. That's what we compare against. Comparing against the 2020 reference
would count the virus's real mutations as "errors".

### 9.2 How it was downloaded

```bash
mkdir -p data/real && cd data/real
B=https://ftp.sra.ebi.ac.uk/vol1
# sample 1: fragments (two files, one per end of each DNA piece) and official genome
curl -O $B/fastq/ERR420/009/ERR4204399/ERR4204399_1.fastq.gz
curl -O $B/fastq/ERR420/009/ERR4204399/ERR4204399_2.fastq.gz
curl -o ERR4204399.consensus.fasta.gz $B/analysis/ERZ188/ERZ18889613/PHEC-1B679.200417_NS552079_0030_AHNMCFBGXF.climb.ena-a.fasta.gz
# sample 2
curl -O $B/fastq/ERR423/002/ERR4239992/ERR4239992_1.fastq.gz
curl -O $B/fastq/ERR423/002/ERR4239992/ERR4239992_2.fastq.gz
curl -o ERR4239992.consensus.fasta.gz $B/analysis/ERZ182/ERZ18277375/NIRE-1021CE.270420_M70780_0023.climb.ena-a.fasta.gz
# unzip; put both ends of each run into one file
for r in ERR4204399 ERR4239992; do
  gzip -dc ${r}_1.fastq.gz ${r}_2.fastq.gz > $r.fastq
done
gunzip *.consensus.fasta.gz
```

The machine reads each DNA piece from **both ends**, giving two files (`_1`
and `_2`). The two ends come from opposite strands, which our mirror handling
takes care of, so we just put both files together.

### 9.3 What's different about real data

| Problem | In real data | How the program handles it |
|---|---|---|
| Both strands | About half the fragments are mirrors. The very first read pair in sample 1 shows it: read 2 starts `AGGTAACAAACCAACC` (the genome as written), read 1 is the same region from the other strand | Canonical chunks + mirror edges (3.3) |
| Format | FASTQ, 4 lines per fragment | `load_reads` detects it |
| Variable lengths | 31–151 letters | Each fragment's own length is used |
| Unknown letters | Some fragments contain `N`; the official genomes have `N` where nobody could tell | Chunks with `N` skipped; `N` in the genome not scored |
| Very deep coverage | ~800×, so typos repeat many times | Threshold at least a tenth of the typical count |
| Size unknown in advance | 1.6 million different chunks | Hash table grows |
| Uneven coverage | Some regions read far more than others, especially amplicon data | **Not fully handled** (see 9.6) |
| Real biology | Mutations, subgenomic RNA | Shown in the difference list (9.5) |

### 9.4 Results

Command (sample 1; sample 2 is the same with the other run name):

```bash
./assembler data/real/ERR4204399.consensus.fasta -reads data/real/ERR4204399.fastq
```

Compared against each sample's **own official genome**:

| | Simulated (20×) | Real sample 1 (capture) | Real sample 2 (amplicon) |
|---|---|---|---|
| Threshold (automatic) | 3 | 63 | 79 |
| Contigs | 21 | 21 | 39 |
| Longest contig | 6,672 | 7,546 | 4,012 |
| N50 | 2,619 | 5,407 | 2,153 |
| **Genome fraction** | 99.94% | **99.96%** | **94.40%** |
| **Identity** | 99.997% | **99.86%** | **99.99%** |
| Differences | 1 in 30,380 | 42 in 30,330 | 3 in 28,405 |
| Junk contigs | 2 | 8 | 2 |
| Run time | 0.5 s | 2 s | 1 s |

Full outputs: [sample 1](results/ERR4204399_vs_consensus.txt),
[sample 2](results/ERR4239992_vs_consensus.txt).

**Sample 1 (capture), stage by stage:**
- **Stage 2:** 167,610 fragments, 24.9 million letters. The first three
  fragments shown are identical: the lab process copies some pieces several
  times ("duplicates").
- **Stage 3:** 1,656,876 different chunks, but only ~30,000 are real. The
  histogram falls steadily from 1.2 million chunks seen once, and 35,401 chunks
  are seen 11+ times. The real chunks are out at about 637.
- **Stage 4:** typical count 637, threshold 63. 30,016 chunks kept and
  1,626,860 removed (98% of all different chunks were junk).
- **Stage 5:** 60,024 nodes, 60,032 edges, 38 start-like nodes and 23
  branches.
- **Stage 6:** 21 contigs, longest 7,546.
- **Stage 7:** 99.96% recovered; only 13 known letters are missing
  (positions 13–14, 6,845–6,853 and 29,869–29,870).

So on sample 1, the result on real data is **about as good as on simulated
data**.

### 9.5 What the differences reveal

#### Real mutations

Compared against the **2020 Wuhan reference** instead
([sample 1](results/ERR4204399_vs_wuhan.txt),
[sample 2](results/ERR4239992_vs_wuhan.txt)), the difference list contains
real mutations of these UK viruses:

| Sample | Positions found | What they are |
|---|---|---|
| 1 | 241 C>T, 3037 C>T, 14408 C>T, 23403 A>G | The four mutations of lineage **B.1**, the variant that spread across Europe in early 2020. 23403 A>G is the well-known **D614G** spike mutation |
| 1 | 28881–28883 GGG>AAC | Three letters changed together: marks lineage **B.1.1** |
| 1 | 18568 C>T, 25831 C>T | Further mutations of this particular virus |
| 2 | 3037 C>T, 14408 C>T, 23403 A>G | B.1 again (241 is in a region missing from our contigs) |
| 2 | 1059 C>T, 25563 G>T | The two mutations that define Nextstrain clade **20C** |

**Every one of these disappears** when we compare against the samples' official
genomes instead. So our assembler and the official UK pipeline **agree** on all
of them. Built from nothing but scrambled fragments, the program found the
mutations that mark these virus lineages.

#### Subgenomic RNA

In sample 1, 12 of the remaining 17 single-letter differences sit together at
positions 78–96. Checking the contigs there
([results/ERR4204399_leader_junctions.txt](results/ERR4204399_leader_junctions.txt))
shows what happened. After the genome's first ~75 letters, two short contigs
don't continue along the genome. Instead they jump to:
- position **27,894**, the start of gene **ORF8**, and
- position **25,391**, the start of gene **ORF3a**.

That's real biology. Inside infected cells, coronaviruses make shortened copies
of their genes with the first ~75 letters (the "leader") stuck on the front.
These are called **subgenomic RNAs**. They're in the sample, so their chunks
are in the data, and our assembler built them. Compared letter by letter to the
genome, they look like differences, but they're real molecules, not assembly
errors.

The few remaining differences (small clusters around 25,334, 25,826, 26,262
and 27,021) weren't investigated further. Several sit next to gene starts, so
they may be similar junctions.

### 9.6 The weak spot: uneven coverage

Sample 2 recovered only **94.4%** of its genome. The missing stretches are
whole blocks of ~250 letters (positions 59–324, 2,563–2,831, 8,634–8,879, ...).

Amplicon sequencing copies the genome in ~100 overlapping PCR pieces, and some
pieces copy much worse than others. A few regions are read thousands of times;
weak pieces only dozens of times. Our program uses **one threshold for the
whole genome** (79 here), so real chunks in weak regions fall below it and get
deleted along with the typos.

Changing the threshold by hand shows the trade-off
([results/real_threshold_sweep.txt](results/real_threshold_sweep.txt)):

| Threshold | Sample 1: fraction / identity / contigs | Sample 2: fraction / identity / contigs |
|---|---|---|
| 10 | 100.00% / 99.13% / 498 | **99.76%** / 99.48% / 72 |
| 30 | 99.99% / 99.86% / 89 | 93.39% / 97.14% / 24 |
| 60 | 99.96% / 99.86% / 28 | 96.24% / 99.94% / 33 |
| 100 | 99.69% / 99.95% / 18 | 92.89% / 100% / 45 |
| 200 | 98.71% / 99.98% / 64 | 85.32% / 100% / 68 |

- **Low threshold:** keeps the weak regions (99.76% of sample 2) but also keeps
  more junk, giving more pieces and lower identity.
- **High threshold:** only the cleanest data, so identity reaches 100%, but
  weak regions are lost.
- In between (20–30 for sample 2), some contigs get **joined wrongly**
  (identity 97%). Leftover junk chunks connect regions that don't belong
  together.

The automatic threshold picks a safe middle: very few errors, at the cost of
some missing regions in uneven data. Real assemblers handle this with a
**local** threshold, comparing each chunk with the coverage of its
neighbourhood instead of with the whole genome.

Some of sample 2's gaps are also missing in the official genome (its `N`
regions at 19,273 and 22,321). No amplicon there worked at all, so no method
could recover them.

### 9.7 Verdict

- On **capture data** (fairly even coverage), the assembler recovers 99.96% of
  the genome at 99.86% identity, about as well as on simulated data. Most
  remaining differences are real subgenomic RNA.
- On **amplicon data** (very uneven coverage), it's extremely accurate (99.99%)
  but recovers less (94.4%) with the automatic threshold.
- On both, it independently finds the same mutations as the official UK
  pipeline.

---

## 10. Speed (complexity)

**Big-O notation** describes how running time grows with input size. O(n)
means "doubling the input doubles the time". O(n²) means "doubling the input
makes it 4 times slower".

| Stage | Time | Why |
|---|---|---|
| 1. Load genome | O(G) | Reads each letter once (G = genome length) |
| 2. Make/load fragments | O(N × L) | N fragments of L letters |
| 3. Count chunks | O(N × L) | One rolling update and one hash lookup (≈ O(1)) per position; growing the table adds only a constant factor on average |
| 4. Filter | O(table size) | One pass over the table |
| 5. Build graph | O(E) | Two hash lookups per edge (E = number of edges) |
| 6. Walk graph | O(V + E) | Every edge is pushed and popped once (V = number of nodes) |
| 6b. Sort contigs | O(C log C) | C = number of contigs |
| 7. Edit distance | O(n × m) per contig | Every cell of the table; memory only O(m) with two rows |

Measured times: simulated run 0.5 s in total; real capture sample (24.9 million
letters, about 20 million chunks) about 2 s, of which counting chunks is about 1 s.
Peak memory on real data is about 110 MB.

**Why the hash table matters:** without it, counting 20 million chunks would
mean comparing each against all the others, roughly 10¹⁴ steps instead of
about 2 × 10⁷.

---

## 11. Simplifications and next steps

What the program now handles, and what's still simplified:

| Feature | Status |
|---|---|
| Fragments from both strands | ✅ Handled (canonical chunks, mirror edges) |
| FASTQ input, variable lengths, `N` letters | ✅ Handled |
| Hash table growth | ✅ Handled (rehashing when half full) |
| Typo threshold | ✅ Automatic (dip + a tenth of typical count) |
| **Uneven coverage** | ⚠️ One threshold for the whole genome; weak regions of amplicon data are lost. Next step: compare each chunk with its local neighbourhood |
| **Graph clean-up** | ⚠️ Not done. Leftover junk adds short side branches ("tips") and small loops ("bubbles") that break contigs. Next step: remove tips shorter than k before walking |
| Quality scores | ⚠️ Ignored. Could skip letters the machine marked as unsure |
| Adapters and primers | ⚠️ Not trimmed. Lab-added sequences at fragment ends become junk chunks |
| Subgenomic RNA | ℹ️ Assembled like any other molecule. Real virus pipelines filter by mapping to a reference |
| Simulated typos | ℹ️ Substitutions only. Real short-read machines mostly make substitutions, so this is realistic for Illumina |
| One `.c` file, terminal output | ℹ️ By choice. Dashboard to come via the socket server from the C lab project |

---

## 12. Glossary

| Term | Meaning |
|---|---|
| **Genome** | The complete DNA sequence of an organism |
| **FASTA** | Text file format: a `>` label line, then DNA letters |
| **FASTQ** | Sequencer output format: 4 lines per fragment (name, letters, `+`, quality) |
| **Fragment / read** | A short piece of DNA read by the sequencing machine |
| **Coverage** | Average number of fragments covering each letter |
| **k-mer / chunk** | A piece of exactly k letters |
| **Sliding window** | Taking a chunk at every position, moving one letter at a time |
| **Strand** | One of the two chains of the DNA double helix |
| **Reverse complement / mirror** | The same DNA read from the other strand: reversed, with A↔T and C↔G swapped |
| **Canonical chunk** | The smaller (as a number) of a chunk and its mirror, used as their shared name |
| **N** | An unknown letter |
| **Bit** | A single 0 or 1 |
| **2-bit encoding** | Storing each DNA letter in 2 bits (A=00, C=01, G=10, T=11) |
| **Mask** | A number used with `&` to keep only certain bits |
| **XOR (`^`)** | Bit operation: 1 where exactly one side is 1; `e ^ 1` flips the last bit |
| **Hash table** | A lookup table that finds a key in about one step |
| **Hash function** | Turns a key into a slot number |
| **Collision** | Two keys mapping to the same slot |
| **Linear probing** | On a collision, try the next slot |
| **Load factor** | Fraction of slots in use |
| **Rehashing** | Moving every key into a bigger table when the current one gets too full |
| **Histogram** | A count of how many items have each value |
| **Median** | The middle value: half the items are below, half above |
| **Graph** | Nodes joined by edges |
| **Directed edge** | A one-way connection (an arrow) |
| **In/out-degree** | Number of edges entering/leaving a node |
| **Adjacency list** | Each node's list of outgoing edges |
| **Linked list** | Items where each one points to the next |
| **de Bruijn graph** | Graph where chunks are edges between their (k−1)-letter ends |
| **Eulerian path** | A walk using every edge exactly once |
| **Hierholzer's algorithm** | Finds an Eulerian path using a stack |
| **Stack** | A pile: add and remove only at the top |
| **Recursion** | A function that calls itself |
| **Contig** | A continuous stretch of rebuilt DNA |
| **N50** | Contig length at which half the total assembled length is reached |
| **Edit distance** | Fewest single-letter changes turning one string into another |
| **Dynamic programming** | Solving a problem by building up answers to smaller versions |
| **Genome fraction** | Percentage of the genome present in the contigs |
| **Identity** | Percentage of output letters that match the genome |
| **SRA / ENA** | Free public archives of real sequencing data (US / Europe) |
| **Consensus genome** | A sample's sequence as determined by a standard analysis pipeline |
| **Targeted capture** | Lab method: probes fish virus DNA out of a sample |
| **Amplicon sequencing** | Lab method: the genome is copied in many overlapping pieces by PCR |
| **Mutation** | A letter that differs from the reference genome |
| **Lineage / clade** | A family of related virus variants, defined by shared mutations |
| **Subgenomic RNA** | Shortened gene copies that coronaviruses make, with a leader sequence on the front |
| **Seed** | Starting number for random numbers; same seed gives same results |
| **Pointer** | A variable holding a memory address |
| **malloc / free** | Ask for memory / give it back |
| **Big-O** | How running time grows as input grows |
