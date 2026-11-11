#!/usr/bin/env bash
# Screen read-based homopolymer protection (for FASTQ, which has no alignment) on ONT HAC chr1-6
# (DATASET=ont_hac, the default) or ONT SUP HG003 autosomes (DATASET=ont), and FASTQ's PacBio
# --max (one quality per read) on PacBio Revio HG002 (DATASET=pb_revio) and
# PacBio HiFi HG002, which has no `rq` (DATASET=pb):
# encode each variant from the aligned BAM, call small variants genome-wide with lrpipe snv
# (Clair3), score SNPs and indels with vcfeval in all scoring regions and in the long-homopolymer
# stratum, and delete the archive. Aligners ignore base qualities, so a binned FASTQ aligns as the
# binned BAM does, and these calls are what the FASTQ would give.
#
# Variants: readN (8 ONT bins, native qualities in read runs of N bp or more), nohp (8 bins, no
# protection), max_meanerr and max_meanphred (PacBio --max with each read's mean error, or mean
# Phred score, in place of the preset's). `original` and `default` reuse the validation calls.
# Resumable; each run's rows replace earlier ones in $W/eval/screen.tsv.
#
# usage: screen_read_homopolymers.sh VARIANT...
set -euo pipefail
. "$(dirname "$0")/env.sh"
need MODELS
DATASET=${DATASET:-ont_hac}
REF=${ONT_REF:-} PLATFORM=ont SDF=ont PACK_ARGS=()
case $DATASET in
    ont_hac)
        need ONT_HAC_BAM ONT_REF; BAM=$ONT_HAC_BAM
        MODEL=$MODELS/r1041_e82_400bps_hac_v500 CONTIGS=chr1,chr2,chr3,chr4,chr5,chr6 SAMPLE=HG002
        SCOPES="regions:eval/ont_hac.regions.bed longhp:eval/ont_hac.longhp.bed" ;;
    ont)
        need ONT_BAM ONT_REF ONT_AUTOSOMES_BED ONT_LONGHP_BED; BAM=$ONT_BAM
        MODEL=$MODELS/r1041_e82_400bps_sup_v500 CONTIGS=$(printf 'chr%s,' {1..22}) SAMPLE=HG003
        CONTIGS=${CONTIGS%,}
        SCOPES="autosomes:$ONT_AUTOSOMES_BED longhp:$ONT_LONGHP_BED" ;;
    pb_revio)
        need PB_REVIO_BAM PB_REF; BAM=$PB_REVIO_BAM REF=$PB_REF
        MODEL=$MODELS/hifi_revio CONTIGS=$(printf 'chr%s,' {1..22}) SAMPLE=HG002
        CONTIGS=${CONTIGS%,} PLATFORM=hifi SDF=pb PACK_ARGS=(--max)
        SCOPES="regions:eval/pb_revio.regions.bed" ;;
    pb)
        # eval/pb.autosomes.bed links to GIAB v4.2.1's HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed.
        need PB_BAM PB_REF; BAM=$PB_BAM REF=$PB_REF
        MODEL=$MODELS/hifi_revio CONTIGS=$(printf 'chr%s,' {1..22}) SAMPLE=HG002
        CONTIGS=${CONTIGS%,} PLATFORM=hifi SDF=pb PACK_ARGS=(--max)
        SCOPES="autosomes:eval/pb.autosomes.bed" ;;
    *) echo "unknown DATASET $DATASET" >&2; exit 1 ;;
esac
cd "$W"; mkdir -p logs archives calls eval policies

log() { echo "[$(date -u +%FT%TZ)] $*" >> logs/events.log; }
timed() { local label=$1; shift; log "start $label"; /usr/bin/time -v -o "logs/$label.time" "$@" > "logs/$label.log" 2>&1; log "done $label ($(grep Elapsed "logs/$label.time" | awk '{print $NF}'))"; }
once() { [ -e "logs/$1.done" ] && return 1; return 0; }

for v in "$@"; do
    name=$DATASET.$v archive=archives/$DATASET.$v.cram
    case $v in
        original | default) ;;
        read* | nohp | max_meanerr | max_meanphred)
            if once "encode_$name"; then
                # Newer built-in policies name their homopolymer source; drop it.
                if [ "${v#max_}" != "$v" ]; then
                    source=$([ "$v" = max_meanerr ] && echo mean-error || echo mean-phred)
                    "$LRSQUEEZE" policy --platform pacbio --max | sed "s/^source = .*/source = \"$source\"/" > "policies/$v.toml"
                    grep -q "source = \"$source\"" "policies/$v.toml"
                elif [ "$v" = nohp ]; then
                    "$LRSQUEEZE" policy --platform ont | grep -v '^preserve_homopolymers\|^homopolymers_in' > "policies/$v.toml"
                else
                    "$LRSQUEEZE" policy --platform ont | grep -v '^homopolymers_in' |
                        sed "s/^preserve_homopolymers = .*/preserve_homopolymers = ${v#read}\nhomopolymers_in = \"read\"/" > "policies/$v.toml"
                    grep -q 'homopolymers_in = "read"' "policies/$v.toml"
                fi
                timed "encode_$name" "$LRSQUEEZE" pack "$BAM" -T "$REF" -c "policies/$v.toml" "${PACK_ARGS[@]}" -o "$archive" -@ 16
                echo -e "$name\t$(stat -c %s "$archive")" >> eval/sizes.tsv
                touch "logs/encode_$name.done"
            fi
            if once "call_$name"; then
                rm -rf "calls/$name"
                timed "call_$name" "$LRPIPE" snv run -b "$archive" -f "$REF" -m "$MODEL" -o "calls/$name" -t 14 \
                    -p "$PLATFORM" --ctg-name "$CONTIGS" --sample-name "$SAMPLE" --remove-intermediate-dir
                touch "logs/call_$name.done"
            fi
            rm -f "$archive" "$archive.crai" ;;
        *) echo "unknown variant $v" >&2; exit 1 ;;
    esac
    for k in snps indels; do
        # Each step is redone when its input is newer, so calls made again are scored again.
        calls=calls/$name/merge_output.vcf.gz c="eval/$name.$k.vcf.gz"
        [ "$c" -nt "$calls" ] || { bcftools view -v "$k" -f PASS -Oz -o "$c" "$calls"; tabix -f -p vcf "$c"; }
        for scope in $SCOPES; do
            s=${scope%%:*} bed=${scope#*:} out="eval/$name.$s.$k"
            if [ ! "$out/summary.txt" -nt "$c" ]; then
                rm -rf "$out"
                timed "vcfeval_${name}_${s}_$k" rtg vcfeval -b "eval/truth.$DATASET.$k.vcf.gz" -c "$c" -e "$bed" \
                    -t "sdf/$SDF" -o "$out" --threads 14
            fi
            # Replace this row if an earlier run wrote one.
            touch eval/screen.tsv
            grep -v "^$name	$s	$k	" eval/screen.tsv > eval/screen.tsv.new || true
            echo -e "$name\t$s\t$k\t$(awk '$1=="None"{print $2,$3,$4,$5}' OFS='\t' "$out/summary.txt")" >> eval/screen.tsv.new
            mv eval/screen.tsv.new eval/screen.tsv
        done
    done
done
log "READ HP SCREEN DONE $DATASET $*"
