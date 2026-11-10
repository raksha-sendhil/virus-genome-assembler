#!/bin/bash
# Main demo: pick a real sample, rebuild its genome from the real fragments,
# then compare the result with that patient's official genome.
#
# Usage (in WSL):  bash demo.sh

cd "$(dirname "$0")" || exit 1
CACHE=~/dsa      # WSL reads the C: drive slowly, so runs use a copy in Linux

# build the program if it is missing or older than the source
if [ ! -x assembler ] || [ assembler.c -nt assembler ]; then
    gcc -O2 -o assembler assembler.c || exit 1
fi

pause() { read -r -p "Press Enter to continue..."; echo; }

# ---- choose the input ----
echo
echo "=== Rebuilding a virus genome from real DNA fragments ==="
echo
echo "Choose a sample (real SARS-CoV-2 sequencing runs, UK, 2020):"
samples=()
for f in data/real/*.consensus.fasta; do
    [ -e "$f" ] && samples+=("$(basename "$f" .consensus.fasta)")
done
i=1
for s in "${samples[@]}"; do
    info=$(cat "data/real/$s.info" 2>/dev/null)
    status=ready
    [ -s "data/real/$s.fastq" ] || status="downloads first"
    printf "  %d) %-11s %-9s %-10s (%s)\n" $i "$s" ${info:-? ?} "$status"
    i=$((i + 1))
done
printf "  %d) Download a new sample by ID\n" $i
read -r -p "> " choice

if [ "$choice" = "$i" ]; then
    read -r -p "Sample ID (for example ERR4204296): " run
elif [[ "$choice" =~ ^[0-9]+$ ]] && [ "$choice" -ge 1 ] && [ "$choice" -lt $i ]; then
    run=${samples[$((choice - 1))]}
else
    echo "Not a valid choice."; exit 1
fi

if [ ! -s "data/real/$run.fastq" ] || [ ! -s "data/real/$run.consensus.fasta" ]; then
    echo
    bash get_sample.sh "$run" || exit 1
fi
if [ ! -s "data/real/$run.consensus.fasta" ]; then
    echo "This sample has no official genome to compare with."; exit 1
fi

# ---- step 1: the input ----
echo
echo "=== Step 1: The input ==="
n=$(( $(wc -l < "data/real/$run.fastq") / 4 ))
echo "Sample $run: $n real fragments from a sequencing machine. The first five:"
awk 'NR % 4 == 2' "data/real/$run.fastq" | head -5 | cut -c1-80 | sed 's/^/  /'
echo
echo "The answer (this patient's official genome) is only used at the end, to grade the result."
pause

# ---- step 2: run the pipeline ----
echo "=== Step 2: Rebuild the genome ==="
mkdir -p "$CACHE" results
cp -u "data/real/$run.fastq" "$CACHE/"
out="results/${run}_vs_consensus.txt"
./assembler "data/real/$run.consensus.fasta" -reads "$CACHE/$run.fastq" | tee "$out"
cp contigs.fasta "data/real/$run.contigs.fasta"
pause

# ---- step 3: the result compared with the answer ----
echo "=== Step 3: Result compared with the answer ==="
field() { sed -n "s/.*$1 *: //p" "$out" | head -1; }
echo "  Sample           : $run"
echo "  Genome recovered : $(field '(genome fraction)')"
echo "  Letters correct  : $(field '(identity)')"
echo "  Pieces (contigs) : $(field 'Contigs' | sed 's/ *(saved.*//'), longest $(field 'Longest')"
echo "  Rebuilt DNA      : data/real/$run.contigs.fasta"
echo "  Full output      : $out"
echo
