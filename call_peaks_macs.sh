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

source ${PYTHON_ENV}

# Select MACS command based on configured version
macs_cmd="${MACS_VERSION:-macs3}"
if ! command -v "${macs_cmd}" &>/dev/null; then
    echo "ERROR: ${macs_cmd} not found. Install with: pip install ${MACS_VERSION:-MACS3}" >&2
    exit 1
fi


echo "Calling peaks with ${macs_cmd} (genome: ${macs_genome})..."

# Call peaks with MACS
${macs_cmd} callpeak \
    -t "${fragments_file}" \
    -f FRAG \
    -g "${macs_genome}" \
    -n peaks \
    --outdir "${output_dir}" \
    2>&1 | tee "${output_dir}/macs_log.txt"

# Convert MACS narrowPeak output to the BED format expected by downstream scripts.
# narrowPeak is BED6+4; downstream expects peaks_w_blacklistregion.bed (BED format).
macs_peaks="${output_dir}/peaks_peaks.narrowPeak"
output_peaks="${output_dir}/peaks_w_blacklistregion.bed"

if [[ -f "${macs_peaks}" ]]; then
    # Cap the integer score at 1000 for compatibility with downstream tools. Ref from MACS docs: https://macs3-project.github.io/MACS/docs/callpeak.html
    awk -v OFS="\t" '{$5=$5>1000?1000:$5} {print}' "${macs_peaks}" > "${output_peaks}"
    echo "Peak calling complete: $(wc -l < "${output_peaks}") peaks written to ${output_peaks}"
else
    echo "ERROR: MACS peak file not found: ${macs_peaks}" >&2
    exit 1
fi


echo "Call peaks with ${macs_cmd} .... Done"
