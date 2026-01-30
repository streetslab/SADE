#!/bin/bash
# Call peaks using MACS2/MACS3 directly from a fragments file.
# This avoids the expensive samtools sort + BAM processing required by Genrich.
#
# MACS calls peaks from Tn5 insertion sites extracted from the fragments file,
# using standard ATAC-seq parameters (--nomodel --shift -75 --extsize 150).

set -o errexit
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

Usage="Usage: $0
  -f <fragments_file. Required. TSV with columns: chr, start, end, barcode, count>
  -o <output_dir. Required>
  -g <genome. Required: 'hg38', 'mm10' etc.>
  "

while getopts "f:o:g:" opt; do
  case $opt in
    f) fragments_file="$OPTARG" ;;
    o) output_dir="$OPTARG"
       output_dir=${output_dir%/}
       ;;
    g) genome="$OPTARG" ;;
    \?)
      echo -e "Invalid option: -$OPTARG \n$Usage" >&2
      exit 1
      ;;
  esac
done

if [[ -z ${fragments_file} || -z ${output_dir} || -z ${genome} ]]; then
    echo "Missing required arguments."
    echo -e "$Usage" >&2
    exit 1
fi

mkdir -p "${output_dir}"

# Map genome name to MACS effective genome size parameter
case "${genome}" in
    hg38|hg19) macs_genome="hs" ;;
    mm10|mm9)  macs_genome="mm" ;;
    *)
        echo "WARNING: Unknown genome '${genome}' for MACS genome size. Using genome name directly." >&2
        macs_genome="${genome}"
        ;;
esac

# Select MACS command based on configured version
macs_cmd="${MACS_VERSION:-macs3}"
if ! command -v "${macs_cmd}" &>/dev/null; then
    echo "ERROR: ${macs_cmd} not found. Install with: pip install ${MACS_VERSION:-MACS3}" >&2
    exit 1
fi

# Convert fragments to Tn5 insertion sites (BED format)
# Each fragment has two Tn5 cut sites at the 5' ends of both reads:
#   - Forward cut: (start, start+1)
#   - Reverse cut: (end-1, end)
tn5_sites_file="${output_dir}/tn5_sites.bed"
echo "Extracting Tn5 insertion sites from fragments..."
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $2+1; print $1, $3-1, $3}' \
    "${fragments_file}" > "${tn5_sites_file}"

echo "Calling peaks with ${macs_cmd} (genome: ${macs_genome})..."

# ATAC-seq MACS parameters:
#   --nomodel: don't build the shifting model (we handle Tn5 sites explicitly)
#   --shift -75 --extsize 150: extend 150bp centered on each Tn5 insertion site
#   --keep-dup all: fragments are already deduplicated by upstream pipeline
#   --call-summits: identify sub-peaks within each peak region
${macs_cmd} callpeak \
    -t "${tn5_sites_file}" \
    -f BED \
    -g "${macs_genome}" \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    --keep-dup all \
    --call-summits \
    -n peaks \
    --outdir "${output_dir}" \
    2>&1 | tee "${output_dir}/macs_log.txt"

# Convert MACS narrowPeak output to the BED format expected by downstream scripts.
# narrowPeak is BED6+4; downstream expects peaks_w_blacklistregion.bed (BED format).
macs_peaks="${output_dir}/peaks_peaks.narrowPeak"
output_peaks="${output_dir}/peaks_w_blacklistregion.bed"

if [[ -f "${macs_peaks}" ]]; then
    # Use first 6 columns (chr, start, end, name, score, strand) for BED compatibility
    cut -f1-6 "${macs_peaks}" > "${output_peaks}"
    echo "Peak calling complete: $(wc -l < "${output_peaks}") peaks written to ${output_peaks}"
else
    echo "ERROR: MACS peak file not found: ${macs_peaks}" >&2
    exit 1
fi

# Clean up intermediate Tn5 sites file
rm -f "${tn5_sites_file}"

echo "Call peaks with ${macs_cmd} .... Done"
