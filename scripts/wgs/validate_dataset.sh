#!/usr/bin/env bash
# Validate lrsqueeze on one aligned dataset: encode default and --max archives, call small variants
# (lrpipe snv, CPU exact) and SVs (lrpipe sv) on the original and both archives, and score small
# variants with RTG vcfeval. Resumable; results append to $W/eval/accuracy.tsv and sizes.tsv.
#
# usage: validate_dataset.sh NAME BAM REF SDF MODEL PLATFORM SAMPLE TRUTH BED CONTIGS
#   SDF: a directory name under $W/sdf/;  CONTIGS: comma-separated (calls and scoring use only these)
set -euo pipefail
NAME=$1 BAM=$2 REF=$3 SDF=$4 MODEL=$5 PLATFORM=$6 SAMPLE=$7 TRUTH=$8 BED=$9 CONTIGS=${10}
. "$(dirname "$0")/env.sh"
cd "$W"; mkdir -p logs archives calls eval sv

log() { echo "[$(date -u +%FT%TZ)] $*" >> logs/events.log; }
timed() { local label=$1; shift; log "start $label"; /usr/bin/time -v -o "logs/$label.time" "$@" > "logs/$label.log" 2>&1; log "done $label ($(grep Elapsed "logs/$label.time" | awk '{print $NF}'))"; }
once() { [ -e "logs/$1.done" ] && return 1; return 0; }

# Scoring regions: the truth BED restricted to CONTIGS.
REGIONS=eval/$NAME.regions.bed
[ -s "$REGIONS" ] || awk -v c="$CONTIGS" 'BEGIN{n=split(c,a,","); for(i=1;i<=n;i++) k[a[i]]=1} k[$1]' "$BED" > "$REGIONS"

for mode in default max; do
    flag=""; [ $mode = max ] && flag=--max
    if once "encode_$NAME.$mode"; then
        timed "encode_$NAME.$mode" "$LRSQUEEZE" pack "$BAM" -T "$REF" -o "archives/$NAME.$mode.cram" -@ 16 $flag
        touch "logs/encode_$NAME.$mode.done"
    fi
done
echo -e "$NAME.original\t$(stat -c %s "$BAM")\n$NAME.default\t$(stat -c %s "archives/$NAME.default.cram")\n$NAME.max\t$(stat -c %s "archives/$NAME.max.cram")" >> eval/sizes.tsv

for cond in original default max; do
    input=$BAM; [ $cond = original ] || input=archives/$NAME.$cond.cram
    if once "call_$NAME.$cond"; then
        rm -rf "calls/$NAME.$cond"
        timed "call_$NAME.$cond" "$LRPIPE" snv run -b "$input" -f "$REF" -m "$MODEL" -o "calls/$NAME.$cond" \
            -t 14 -p "$PLATFORM" --ctg-name "$CONTIGS" --sample-name "$SAMPLE" --remove-intermediate-dir
        touch "logs/call_$NAME.$cond.done"
    fi
    if once "sv_$NAME.$cond"; then
        args=(); IFS=, read -ra cs <<< "$CONTIGS"; for c in "${cs[@]}"; do args+=(-c "$c"); done
        timed "sv_$NAME.$cond" "$LRPIPE" sv -i "$input" --reference "$REF" -t 8 "${args[@]}" -v "sv/$NAME.$cond.vcf.gz"
        touch "logs/sv_$NAME.$cond.done"
    fi
    if once "eval_$NAME.$cond"; then
        for kind in all:snps,indels snps:snps indels:indels; do
            k=${kind%%:*} types=${kind#*:}
            t="eval/truth.$NAME.$k.vcf.gz"
            [ -s "$t" ] || { bcftools view -v "$types" -Oz -o "$t" "$TRUTH"; tabix -f -p vcf "$t"; }
            c="eval/$NAME.$cond.$k.vcf.gz"
            bcftools view -v "$types" -f PASS -Oz -o "$c" "calls/$NAME.$cond/merge_output.vcf.gz"; tabix -f -p vcf "$c"
            rm -rf "eval/$NAME.$cond.regions.$k"
            timed "vcfeval_${NAME}.${cond}_$k" rtg vcfeval -b "$t" -c "$c" -e "$REGIONS" -t "sdf/$SDF" \
                -o "eval/$NAME.$cond.regions.$k" --threads 14
            echo -e "$NAME.$cond\tregions\t$k\t$(awk '$1=="None"{print $2,$3,$4,$5}' OFS='\t' "eval/$NAME.$cond.regions.$k/summary.txt")" >> eval/accuracy.tsv
        done
        touch "logs/eval_$NAME.$cond.done"
    fi
done
log "DATASET DONE $NAME"
