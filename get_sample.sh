#!/bin/bash
# Downloads one real sequencing run from ENA plus its official genome.
#
# Usage (in WSL):  bash get_sample.sh ERR4204399
#
# Writes:
#   data/real/<run>.fastq             the fragments (both read files joined)
#   data/real/<run>.consensus.fasta   the sample's official genome (the answer)
#   data/real/<run>.info              lab method and collection date

run=$1
if [ -z "$run" ]; then echo "Usage: bash get_sample.sh RUN_ACCESSION"; exit 1; fi

api=https://www.ebi.ac.uk/ena/portal/api
cd "$(dirname "$0")" || exit 1
mkdir -p data/real
cd data/real || exit 1

# 1. Where are the fragment files, which patient sample is this, and how was it made?
info=$(curl -sS "$api/filereport?accession=$run&result=read_run&fields=fastq_ftp,sample_accession,library_strategy,collection_date" | tail -1)
files=$(echo "$info" | cut -f2)
sample=$(echo "$info" | cut -f3)
if [ -z "$files" ]; then echo "No read files found for $run"; exit 1; fi
method=$(echo "$info" | cut -f4)
case "$method" in
    Targeted-Capture) method=capture ;;
    AMPLICON)         method=amplicon ;;
esac
echo "$method $(echo "$info" | cut -f5)" > "$run.info"

# 2. Download the fragments and join both read files into one
if [ -s "$run.fastq" ]; then
    echo "Fragments of $run already downloaded."
else
    echo "Downloading fragments of $run ..."
    for f in $(echo "$files" | tr ';' ' '); do
        case "$f" in
            *_1.fastq.gz|*_2.fastq.gz) curl -sS "https://$f" | gzip -dc >> "$run.fastq" ;;
        esac
    done
    if [ ! -s "$run.fastq" ]; then          # single-file run
        curl -sS "https://$(echo "$files" | cut -d';' -f1)" | gzip -dc > "$run.fastq"
    fi
fi

# 3. Find the official genome for the same sample (COG-UK assembly if there is one)
if [ ! -s "$run.consensus.fasta" ]; then
    echo "Looking up the official genome of sample $sample ..."
    list=$(curl -sS -G "$api/search" --data-urlencode "result=analysis" \
             --data-urlencode "query=sample_accession=\"$sample\"" \
             --data-urlencode "fields=analysis_type,submitted_ftp")
    answer=$(echo "$list" | awk -F'\t' '$2=="SEQUENCE_ASSEMBLY"{print $3; exit}' | tr ';' '\n' | grep 'fasta.gz$')
    [ -z "$answer" ] && answer=$(echo "$list" | awk -F'\t' '$2=="COVID19_CONSENSUS"{print $3; exit}' | tr ';' '\n' | grep 'fasta.gz$')
    if [ -z "$answer" ]; then
        echo "No official genome found for this sample; compare with data/sars_cov2.fasta instead."
    else
        curl -sS "https://${answer//#/%23}" | gzip -dc > "$run.consensus.fasta"   # '#' must be written %23 in a URL
    fi
fi

echo
echo "Fragments : data/real/$run.fastq ($(( $(wc -l < "$run.fastq") / 4 )) fragments)"
[ -s "$run.consensus.fasta" ] && echo "Answer    : data/real/$run.consensus.fasta ($(head -1 "$run.consensus.fasta"))"
echo "Method    : $(cat "$run.info")"
