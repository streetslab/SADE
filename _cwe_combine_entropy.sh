#!/usr/bin/env bash
# Usage: ./_cwe_combine_entropy.sh [barcodes_file] [chromosomal_entropy_dir] [output_file]
set -euo pipefail

BARCODES="${1:-_barcodes.csv}"
Chromosomal_Entropy_Dir="${2:-.}"
OUT="${3:-selected_entropy_all.tsv}"


files=( "$Chromosomal_Entropy_Dir"/chr*_calculated_barcode_entropy_df.tsv )
[[ -e "${files[0]}" ]] || { echo "No matching files found" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"


awk -F'[\t,]' -v OFS='\t' '
    # Load barcodes of interest (keeps their original order)
    FNR == NR {
        gsub(/\r/, "")
        if ($0 != "" && !($0 in want)) { want[$0]; order[++n] = $0 }
        next
    }
    # New file: derive chromosome from the file name (text before first "_")
    FNR == 1 {
        f = FILENAME; sub(/.*\//, "", f); sub(/_.*/, "", f)
        chr = f; chrs[++nc] = chr; e = 0; p = 0; l = 0; m = 0
    }
    # Header line: locate the columns by name
    {
        for (i = 1; i <= NF; i++) {
            if ($i == "Entropy")        e = i
            if ($i == "P_closed_state") p = i
            if ($i == "mle_lambda") l = i
            if ($i == "Entropy_open_region") m = i
        }
    }
    e && p && l && m && ($1 in want) && !(chr SUBSEP $1 in ent) && $e != "Entropy" {
        ent[chr, $1] = $e
        pcs[chr, $1] = $p
        mle[chr, $1] = $l
        ent_open[chr, $1] = $m
    }
    END {
        # Header
        line = "barcode"
        for (c = 1; c <= nc; c++) line = line OFS chrs[c] "_Entropy" OFS chrs[c] "_P_closed_state" OFS chrs[c] "_mle_lambda" OFS chrs[c] "_Entropy_open_region"
        print line
        # One row per barcode
        for (r = 1; r <= n; r++) {
            b = order[r]; line = b
            for (c = 1; c <= nc; c++) {
                k = chrs[c] SUBSEP b
                line = line OFS ((k in ent) ? ent[k] : "NA") OFS ((k in pcs) ? pcs[k] : "NA") OFS ((k in mle) ? mle[k] : "NA") OFS ((k in ent_open) ? ent_open[k] : "NA")
            }
            print line
        }
    }
' "$BARCODES" "${files[@]}" > "$OUT"

echo "Wrote $(($(wc -l < "$OUT") - 1)) barcodes x ${#files[@]} chromosomes to $OUT"
