#!/bin/bash


directories=(
    "/Users/vardaspavarde/Desktop/dataForTauras/rawData/day0"
    "/Users/vardaspavarde/Desktop/dataForTauras/rawData/barcode03"
    "/Users/vardaspavarde/Desktop/dataForTauras/rawData/barcode04"
    "/Users/vardaspavarde/Desktop/dataForTauras/rawData/barcode05"
    "/Users/vardaspavarde/Desktop/dataForTauras/rawData/barcode06"
)


reference="/Users/vardaspavarde/Desktop/dataForTauras/PA14_baseline_reference_genome_flye_medaka.fasta"


exec > >(tee log.txt) 2>&1

# Function to process directories
process_directories() {
    for dir in "${directories[@]}"; do
        
        combined_fastq="$dir/combined_reads.fastq"
        cat "$dir"/*.fastq > "$combined_fastq"

        
        sam_output="$dir/output.sam"
        minimap2 -a -x map-ont "$reference" "$combined_fastq" > "$sam_output"

    
        bam_output="$dir/output.bam"
        samtools view -b -o "$bam_output" "$sam_output"

    
        rm "$sam_output"

  
        dir_name=$(basename "$dir")

        # Sort BAM file with the directory name as the prefix
        sorted_bam_output="$dir/${dir_name}.bam"
        samtools sort "$bam_output" -o "$sorted_bam_output"

        
        samtools index "$sorted_bam_output"

        
        rm "$bam_output"

        
        run_nanocaller "$sorted_bam_output"

        echo "Processed and sorted data in directory: $dir with final BAM as ${dir_name}.bam"
    done
}


run_nanocaller() {
    local bam_file=$1

    
    local dir=$(dirname "$bam_file")

   
    cd "$dir" || return

    
    NanoCaller --bam "$(basename "$bam_file")" --ref "$reference" --cpu 10 --mode snps

    
    cd - || return
}

# Main function to run the script
main() {
    # Process directories
    process_directories
}

# Run the main function
main

