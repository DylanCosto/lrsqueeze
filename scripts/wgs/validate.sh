#!/usr/bin/env bash
# Whole-genome validation of lrsqueeze: original BAM vs lrsqueeze archives, called with lrpipe snv
# (Clair3 2.0.3 port, CPU, exact) on autosomes and scored with RTG vcfeval against GIAB v4.2.1.
#
# Stages are resumable (each writes a .done marker). Every command's wall time and peak memory go
# to logs/*.time. Usage: validate.sh [STAGE...]   (default: all stages, in order)
set -euo pipefail

# --- configuration (paths: scripts/wgs/local.env.example) ----------------------------------------
. "$(dirname "$0")/env.sh"
need MODELS ONT_BAM ONT_REF ONT_TRUTH ONT_PRIMARY_BED ONT_AUTOSOMES_BED PB_BAM PB_REF PB_TRUTH PB_BED
THREADS=${THREADS:-14}
AUTOSOMES=$(printf 'chr%s,' {1..22}); AUTOSOMES=${AUTOSOMES%,}
ONT_MODEL=$MODELS/r1041_e82_400bps_sup_v500         # HG003 PAY87954, Dorado SUP v5.0
PB_MODEL=$MODELS/hifi_revio
# -------------------------------------------------------------------------------------------------

mkdir -p "$W"/{logs,archives,calls,eval,sdf}
cd "$W"

timed() { # LABEL CMD...
    local label=$1; shift
    echo "[$(date -u +%FT%TZ)] start $label: $*" >> logs/events.log
    /usr/bin/time -v -o "logs/$label.time" "$@" > "logs/$label.log" 2>&1
    echo "[$(date -u +%FT%TZ)] done $label ($(grep 'Elapsed' "logs/$label.time" | awk '{print $NF}'))" >> logs/events.log
}
stage() { [ -e "logs/$1.done" ] && return 1; echo "[$(date -u +%FT%TZ)] stage $1" >> logs/events.log; return 0; }
done_() { touch "logs/$1.done"; }

sdf() { # NAME FASTA
    stage "sdf_$1" || return 0
    rm -rf "sdf/$1"; timed "sdf_$1" rtg format -o "sdf/$1" "$2"; done_ "sdf_$1"
}

encode() { # NAME BAM REF ARGS...
    local name=$1 bam=$2 ref=$3; shift 3
    stage "encode_$name" || return 0
    timed "encode_$name" "$LRSQUEEZE" pack "$bam" -T "$ref" -o "archives/$name.cram" -@ 16 "$@"
    done_ "encode_$name"
}

call() { # NAME INPUT REF MODEL PLATFORM SAMPLE
    local name=$1 input=$2 ref=$3 model=$4 platform=$5 sample=$6
    stage "call_$name" || return 0
    rm -rf "calls/$name"
    timed "call_$name" "$LRPIPE" snv run -b "$input" -f "$ref" -m "$model" -o "calls/$name" \
        -t "$THREADS" -p "$platform" --ctg-name "$AUTOSOMES" --sample-name "$sample" \
        --remove-intermediate-dir
    done_ "call_$name"
}

evaluate() { # NAME TRUTH BED SDF TAG
    local name=$1 truth=$2 bed=$3 sdfdir=$4 tag=$5
    stage "eval_${name}_$tag" || return 0
    for kind in all:snps,indels snps:snps indels:indels; do
        local k=${kind%%:*} types=${kind#*:}
        local t="eval/truth.$sdfdir.$k.vcf.gz"
        [ -s "$t" ] || { bcftools view -v "$types" -Oz -o "$t" "$truth"; tabix -f -p vcf "$t"; }
        local c="eval/$name.$k.vcf.gz"
        # Made again each time, so calls made again are what is scored.
        bcftools view -v "$types" -f PASS -Oz -o "$c" "calls/$name/merge_output.vcf.gz"; tabix -f -p vcf "$c"
        rm -rf "eval/$name.$tag.$k"
        timed "vcfeval_${name}_${tag}_$k" rtg vcfeval -b "$t" -c "$c" -e "$bed" -t "sdf/$sdfdir" \
            -o "eval/$name.$tag.$k" --threads "$THREADS"
        echo -e "$name\t$tag\t$k\t$(awk '$1=="None"{print $2,$3,$4,$5}' OFS='\t' "eval/$name.$tag.$k/summary.txt")" >> eval/accuracy.tsv
    done
    done_ "eval_${name}_$tag"
}

size_of() { stat -c %s "$1"; }

all_stages() {
    sdf ont "$ONT_REF"
    sdf pb "$PB_REF"
    # ONT HG003: original, default, --max
    encode ont.default "$ONT_BAM" "$ONT_REF"
    call ont.original "$ONT_BAM" "$ONT_REF" "$ONT_MODEL" ont HG003
    call ont.default "archives/ont.default.cram" "$ONT_REF" "$ONT_MODEL" ont HG003
    for n in ont.original ont.default; do
        evaluate $n "$ONT_TRUTH" "$ONT_PRIMARY_BED" ont primary
        evaluate $n "$ONT_TRUTH" "$ONT_AUTOSOMES_BED" ont autosomes
    done
    # PacBio HG002 (no @RG PL in this BAM): original, default
    encode pb.default "$PB_BAM" "$PB_REF" --platform pacbio
    call pb.original "$PB_BAM" "$PB_REF" "$PB_MODEL" hifi HG002
    call pb.default "archives/pb.default.cram" "$PB_REF" "$PB_MODEL" hifi HG002
    for n in pb.original pb.default; do evaluate $n "$PB_TRUTH" "$PB_BED" pb autosomes; done
    # ONT --max last: slowest to encode
    encode ont.max "$ONT_BAM" "$ONT_REF" --max
    call ont.max "archives/ont.max.cram" "$ONT_REF" "$ONT_MODEL" ont HG003
    evaluate ont.max "$ONT_TRUTH" "$ONT_PRIMARY_BED" ont primary
    evaluate ont.max "$ONT_TRUTH" "$ONT_AUTOSOMES_BED" ont autosomes
    # Appended: the other scripts log here the sizes of archives they delete.
    [ -s eval/sizes.tsv ] || echo -e "file\tbytes" > eval/sizes.tsv
    for f in "$ONT_BAM" "$PB_BAM" archives/*.cram; do echo -e "$(basename "$f")\t$(size_of "$f")"; done >> eval/sizes.tsv
    echo "[$(date -u +%FT%TZ)] ALL DONE" >> logs/events.log
}

if [ $# -eq 0 ]; then all_stages; else for s in "$@"; do "$s"; done; fi
