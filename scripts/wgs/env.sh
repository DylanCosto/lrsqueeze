# Settings shared by the validation scripts, sourced by each. Paths come from the environment or
# from local.env beside this file (untracked: copy local.env.example and fill it in). Tools (rtg,
# bcftools, tabix, samtools, lrpipe) are found on PATH, searched after TOOLS_BIN when it is set.
ENV_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
if [ -f "$ENV_DIR/local.env" ]; then . "$ENV_DIR/local.env"; fi
: "${W:?set W to the work directory (see scripts/wgs/local.env.example)}"
LRSQUEEZE=${LRSQUEEZE:-"$ENV_DIR/../../target/release/lrsqueeze"}
LRPIPE=${LRPIPE:-lrpipe}
if [ -n "${TOOLS_BIN:-}" ]; then PATH=$TOOLS_BIN:$PATH; fi
export PATH OMP_NUM_THREADS=1 RTG_MEM=${RTG_MEM:-12G}

# Fails unless every named variable is set.
need() {
    local v
    for v in "$@"; do
        if [ -z "${!v:-}" ]; then echo "set $v (see scripts/wgs/local.env.example)" >&2; exit 1; fi
    done
}
