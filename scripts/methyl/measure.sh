#!/usr/bin/env bash
# Measure methylation-storage variants on a BAM slice: size after `lrsqueeze pack` (default mode)
# and CpG pileup agreement with the original (modkit 0.6.4).
# usage: measure.sh IN.bam REF.fa OUTDIR
set -euo pipefail
IN=$1 REF=$2 OUT=$3
LRSQUEEZE=${LRSQUEEZE:-"$(cd "$(dirname "$0")/../.." && pwd)/target/release/lrsqueeze"}
# python3 (with pysam), samtools and modkit, from PATH unless named.
PY=${PYTHON:-python3}
SAM=${SAMTOOLS:-samtools}
MODKIT=${MODKIT:-modkit}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"; cd "$OUT"

for v in original none no6ma cpg cpg_no6ma; do
    if [ "$v" = original ]; then src=$IN; else
        [ -s "$v.bam" ] || "$PY" -I "$HERE/prototype.py" "$v" "$IN" "$v.bam" "$REF"
        src=$v.bam
    fi
    [ -s "$v.cram" ] || "$LRSQUEEZE" pack "$src" -T "$REF" -o "$v.cram" -@ 8 --platform ont 2> "$v.pack.log"
    "$SAM" cram-size -v "$v.cram" > "$v.streams"
    echo -e "$v\t$(stat -c %s "$v.cram")" >> sizes.tsv
done

# CpG pileups (combined strands, fixed threshold) from the archives as stored.
for v in original cpg cpg_no6ma; do
    [ -s "$v.cpg.bed.gz" ] || "$MODKIT" pileup "$v.cram" - --ref "$REF" --cpg --combine-strands \
        --modified-bases 5mC 5hmC --filter-threshold C:0.8 -t 8 --log-filepath "$v.cpg.log" | gzip -1 > "$v.cpg.bed.gz"
done
echo MEASURE DONE
