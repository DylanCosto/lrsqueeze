#!/usr/bin/env bash
# SV calls (lrpipe sv, a Sniffles2 2.8.1 port) on original BAMs and lrsqueeze archives, for
# concordance. Usage: sv.sh   (runs every pair below; outputs in $W/sv/)
set -euo pipefail
. "$(dirname "$0")/env.sh"
need ONT_BAM ONT_REF PB_BAM PB_REF
THREADS=${THREADS:-4}
cd "$W"; mkdir -p sv logs

sv() { # NAME INPUT REF
    [ -s "sv/$1.vcf.gz" ] && return 0
    echo "[$(date -u +%FT%TZ)] start sv_$1" >> logs/events.log
    /usr/bin/time -v -o "logs/sv_$1.time" "$LRPIPE" sv -i "$2" --reference "$3" -t "$THREADS" \
        -v "sv/$1.vcf.gz" > "logs/sv_$1.log" 2>&1
    echo "[$(date -u +%FT%TZ)] done sv_$1" >> logs/events.log
}

sv ont.original "$ONT_BAM" "$ONT_REF"
# ont.default is 8 bins without homopolymer protection; ont.q8hp is the shipped default.
sv ont.default archives/ont.default.cram "$ONT_REF"
sv ont.q8hp archives/ont.q8hp.cram "$ONT_REF"
sv ont.max archives/ont.max.shipped.cram "$ONT_REF"
sv pb.original "$PB_BAM" "$PB_REF"
sv pb.default archives/pb.default.cram "$PB_REF"
echo "[$(date -u +%FT%TZ)] SV DONE" >> logs/events.log
