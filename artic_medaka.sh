#!/bin/bash
set -e

##########################################
# A script to automate the processing of SARS-CoV-2
# genome assembly directly from raw output FASTQ files
# using the artic bioinformatics protocol
# https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
# Version 0.0.1
# Written by George Githinji
# License GNU
#########################################

USAGE="USAGE $0 -f sample_manifest_file.csv -r run_name -i input_fastq_files -p path/to/primer/scheme -n normalisation -t threads -k bam"

echo "Activating the artic workflow...."

# Activate conda enviroment, assumes that miniconda is installed in the home directory
source "${HOME}/miniconda3/etc/profile.d/conda.sh"

# Activate the artic conda enviroment,
conda activate artic

SCHEME_DIR="${HOME}/artic-ncov2019/primer_schemes" # Artic scheme directory
MEDAKA_MODEL="r941_min_high_g360"                  # medaka model
NORMALISE_VALUE=200                                # default normalisation value
MIN_LENGTH=300                                     # default minimum amplicon length
MAX_LENGTH=750                                     # default maximum amplicon length
KEEP=""

intermediate_files=("*.vcf"  "*.gz" "*.tbi" "*.txt" "*.er" "*.depths" "*.muscle.*" "*.preconsensus.*" "*.bam" "*.bai" "*.hdf")
intermediate_files_less_bam=("*.vcf"  "*.gz" "*.tbi" "*.txt" "*.er" "*.depths" "*.muscle.*" "*.preconsensus.*" "*.hdf")

# A function that deletes all intermediate files
delete_all(){
    for ((i=0; i<${#intermediate_files[@]}; i++)); do
        find . -maxdepth 1 -name "${intermediate_files[$i]}" -type f -delete || continue
    done
}

# Delete all finds except the bam files
keep_bam(){
    for ((i=0; i<${#intermediate_files_less_bam[@]}; i++)); do
        find . -maxdepth 1 -name "${intermediate_files_less_bam[$i]}" -type f -delete || continue
    done
}

while getopts f:r:i:p:n:t:k:h opt
do
    case "$opt" in
        f) SAMPLE_NAMES="$OPTARG";;
        r) RUN_NAME="$OPTARG";;
        i) FASTQ_FILES="$OPTARG";;
        p) SCHEME_DIR="$OPTARG";;
        n) NORMALISE_VALUE="$OPTARG";;
        t) THREADS="$OPTARG";;
        k) KEEP="$OPTARG";;
        h) echo $USAGE
            exit 1;;
    esac
done

# enforce "run name", "sample manifest file path", "fastq directory path, and number of threads command option are provided by the user
if [[ $RUN_NAME == "" || $SAMPLE_NAMES == "" || $FASTQ_FILES == "" || $THREADS == "" ]]; then
    echo $USAGE
    exit 0
fi

# read the sample manifest file to process the samples
# TODO: optimize this loop to run in parallel
while IFS=, read -r barcode sample_name
do
	if [ -z "${barcode}" ] || [ -z "${sample_name}" ]; then

        echo "Missing barcode or matching sample_id, skipping this entry"

    else

        echo "Processing "$sample_name" "

        test -z "$barcode" && continue

	    artic guppyplex \
            --skip-quality-check \
            --min-length "$MIN_LENGTH" \
            --max-length "$MAX_LENGTH" \
            --directory "$FASTQ_FILES/$barcode" \
            --prefix "$RUN_NAME" && echo "running artic minion ..." || echo "Missing files, failed to process $barcode $sample_name"

        artic minion --medaka \
            --medaka-model "$MEDAKA_MODEL" \
            --normalise "$NORMALISE_VALUE" \
            --threads "$THREADS" \
            --scheme-directory "$SCHEME_DIR" \
            --read-file "${RUN_NAME}_${barcode}.fastq" \
            nCoV-2019/V4.1 "${sample_name}_${RUN_NAME}" || continue

        if [ "$KEEP" == "none" ]; then
            echo "Removing intermediate files"
            delete_all
        elif [ "$KEEP" == "bam" ]; then
            echo "Retained the BAM and corresponding bam in file"
            keep_bam
        elif [ "$KEEP" == "all" ]; then
            echo "Retaining intermediate files..."
        fi
    fi

done <$SAMPLE_NAMES

cat *.consensus.fasta >"$RUN_NAME.fasta"

echo "Done!"

conda deactivate

echo "Deactivating artic workflow....."
