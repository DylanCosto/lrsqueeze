#!/usr/bin/env bash
# Screen ONT homopolymer-protection variants on HAC chr1-6: encode each variant, call small
# variants genome-wide (calls restricted to a BED phase with fewer SNPs and are not comparable),
# score SNPs and indels with vcfeval in all scoring regions and in the long-homopolymer stratum,
# and delete the archive. `original` and `default` reuse validate_dataset.sh's calls.
# Resumable; results append to $W/eval/screen.tsv.
#
# usage: screen_homopolymers.sh VARIANT...   (VARIANT = original | default | lLEN_fFLANK, e.g. l5_f2)
set -euo pipefail
LRSQUEEZE=${LRSQUEEZE:?"set LRSQUEEZE to a pre-0.1.0 build with homopolymer-flank.patch applied"}
. "$(dirname "$0")/../../scripts/wgs/env.sh"
need ONT_HAC_BAM ONT_REF MODELS
BAM=$ONT_HAC_BAM REF=$ONT_REF MODEL=$MODELS/r1041_e82_400bps_hac_v500
cd "$W"; mkdir -p logs archives calls eval policies

log() { echo "[$(date -u +%FT%TZ)] $*" >> logs/events.log; }
timed() { local label=$1; shift; log "start $label"; /usr/bin/time -v -o "logs/$label.time" "$@" > "logs/$label.log" 2>&1; log "done $label ($(grep Elapsed "logs/$label.time" | awk '{print $NF}'))"; }
once() { [ -e "logs/$1.done" ] && return 1; return 0; }

for v in "$@"; do
    name=ont_hac.$v archive=archives/ont_hac.$v.cram
    case $v in
        original | default) ;;
        l*_f*)
            len=${v#l}; len=${len%%_*}; flank=${v#*_f}
            if once "encode_$name"; then
                "$LRSQUEEZE" policy --platform ont \
                    | sed "s/^preserve_homopolymers = .*/preserve_homopolymers = $len\nhomopolymer_flank = $flank/" \
                    > "policies/$v.toml"
                grep -q "homopolymer_flank = $flank" "policies/$v.toml"
                timed "encode_$name" "$LRSQUEEZE" pack "$BAM" -T "$REF" -c "policies/$v.toml" -o "$archive" -@ 16
                echo -e "$name\t$(stat -c %s "$archive")" >> eval/sizes.tsv
                touch "logs/encode_$name.done"
            fi
            if once "call_$name"; then
                rm -rf "calls/$name"
                timed "call_$name" "$LRPIPE" snv run -b "$archive" -f "$REF" -m "$MODEL" -o "calls/$name" -t 14 \
                    -p ont --ctg-name chr1,chr2,chr3,chr4,chr5,chr6 --sample-name HG002 --remove-intermediate-dir
                touch "logs/call_$name.done"
            fi
            rm -f "$archive" "$archive.crai" ;;
        *) echo "unknown variant $v" >&2; exit 1 ;;
    esac
    for k in snps indels; do
        c="eval/$name.$k.vcf.gz"
        [ -s "$c" ] || { bcftools view -v "$k" -f PASS -Oz -o "$c" "calls/$name/merge_output.vcf.gz"; tabix -f -p vcf "$c"; }
        for scope in regions:eval/ont_hac.regions.bed longhp:eval/ont_hac.longhp.bed; do
            s=${scope%%:*} bed=${scope#*:} out="eval/$name.${scope%%:*}.$k"
            if [ ! -s "$out/summary.txt" ]; then
                rm -rf "$out"
                timed "vcfeval_${name}_${s}_$k" rtg vcfeval -b "eval/truth.ont_hac.$k.vcf.gz" -c "$c" -e "$bed" \
                    -t sdf/ont -o "$out" --threads 14
            fi
            grep -q "^$name	$s	$k	" eval/screen.tsv 2>/dev/null ||
                echo -e "$name\t$s\t$k\t$(awk '$1=="None"{print $2,$3,$4,$5}' OFS='\t' "$out/summary.txt")" >> eval/screen.tsv
        done
    done
done
log "SCREEN DONE $*"
