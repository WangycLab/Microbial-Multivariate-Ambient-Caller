#!/usr/bin/env bash
# build_gene_module.sh — build the gene -> KEGG module map (gene_module.tsv) from cds.fa
#
# cds.fa is written by mva_extract_features.py and carries gene-tagged headers
# (><seqid>#<gene_id>). The map is consumed by MVAcaller.py --kegg-map.
#
# Steps:
#   1. seqkit translate cds.fa > protein.faa
#   2. emapper.py -i protein.faa --itype proteins -o eggnog_out --cpu N --dbmem
#   3. grep -v '^#' eggnog_out.emapper.annotations | cut -f1,14 | sed 's/^[^#]*#//'
#      -> gene_module.tsv (column 14 is KEGG_Module; sed drops the <seqid># prefix)
set -euo pipefail

CDS=""; OUT=""; DATA_DIR=""; PREFIX="eggnog_out"
THREADS=64; DBMEM=1; FORCE=0; DRY_RUN=0

log()  { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die()  { printf '\nERROR: %s\n' "$*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "missing dependency: $1"; }

usage() {
cat <<'EOF'
Usage: build_gene_module.sh -i ref/cds.fa [-o ref/gene_module.tsv] [options]

Required
  -i, --cds FILE         cds.fa from mva_extract_features.py (gene-tagged headers)

Options
  -o, --out FILE         Output map [<cds dir>/gene_module.tsv]
  -d, --data-dir DIR     eggNOG database directory passed to emapper.py --data_dir
  -p, --prefix STR       emapper.py output prefix [eggnog_out]
  -t, --threads INT      emapper.py --cpu [64]
      --no-dbmem         Do not pass --dbmem to emapper.py (keep the DB on disk)
      --force            Rebuild even if the output already exists
      --dry-run          Print every command without executing it
  -h, --help             Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  if [[ "$1" == --*=* ]]; then set -- "${1%%=*}" "${1#*=}" "${@:2}"; fi
  case "$1" in
    -i|--cds)       CDS="$2"; shift 2 ;;
    -o|--out)       OUT="$2"; shift 2 ;;
    -d|--data-dir)  DATA_DIR="$2"; shift 2 ;;
    -p|--prefix)    PREFIX="$2"; shift 2 ;;
    -t|--threads)   THREADS="$2"; shift 2 ;;
    --no-dbmem)     DBMEM=0; shift ;;
    --force)        FORCE=1; shift ;;
    --dry-run)      DRY_RUN=1; shift ;;
    -h|--help)      usage; exit 0 ;;
    *)              die "unknown option: $1 (see --help)" ;;
  esac
done

[[ -n "$CDS" ]] || { usage; die "missing required argument: --cds"; }
[[ -r "$CDS" ]] || die "cannot read: $CDS"
[[ -n "$OUT" ]] || OUT="$(cd "$(dirname "$CDS")" && pwd)/gene_module.tsv"
if [[ -n "$DATA_DIR" && ! -d "$DATA_DIR" ]]; then die "not a directory: $DATA_DIR"; fi

REF_DIR="$(dirname "$OUT")"
mkdir -p "$REF_DIR"
PROT="$REF_DIR/protein.faa"
PREFIX_PATH="$REF_DIR/$PREFIX"
ANNOT="$PREFIX_PATH.emapper.annotations"

if [[ -s "$OUT" && $FORCE == 0 ]]; then
  log "reusing $OUT"
  exit 0
fi

if [[ $DRY_RUN == 0 ]]; then need seqkit; need emapper.py; fi

# ------------------------------------------------------------ 1. CDS -> protein
log "Step 1 | seqkit translate $CDS"
if [[ $DRY_RUN == 1 ]]; then
  log "seqkit translate $CDS > $PROT"
else
  seqkit translate "$CDS" > "$PROT"
  log "wrote $PROT ($(grep -c '^>' "$PROT") proteins)"
fi

# --------------------------------------------------------------- 2. eggNOG-mapper
log "Step 2 | emapper.py annotation"
EMAPPER_ARGS=(-i "$PROT" --itype proteins -o "$PREFIX_PATH" --cpu "$THREADS")
if [[ -n "$DATA_DIR" ]]; then EMAPPER_ARGS+=(--data_dir "$DATA_DIR"); fi
if [[ $DBMEM == 1 ]]; then EMAPPER_ARGS+=(--dbmem); fi
if [[ $DRY_RUN == 1 ]]; then
  log "emapper.py ${EMAPPER_ARGS[*]}"
else
  emapper.py "${EMAPPER_ARGS[@]}"
  [[ -s "$ANNOT" ]] || die "emapper.py produced no $ANNOT"
fi

# ------------------------------------------------- 3. KEGG module column -> map
log "Step 3 | gene -> KEGG module map"
if [[ $DRY_RUN == 1 ]]; then
  log "grep -v '^#' $ANNOT | cut -f1,14 | sed 's/^[^#]*#//' > $OUT"
else
  grep -v '^#' "$ANNOT" \
    | cut -f1,14 \
    | sed 's/^[^#]*#//' \
    | awk -F'\t' 'NF == 2 && $2 != "" && $2 != "-"' > "$OUT.tmp"
  mv "$OUT.tmp" "$OUT"
  log "wrote $OUT ($(wc -l < "$OUT") genes with a KEGG module)"
fi

log "Done."
