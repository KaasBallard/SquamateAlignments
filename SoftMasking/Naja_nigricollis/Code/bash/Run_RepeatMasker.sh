#!/bin/bash
: <<'ScriptDescription'
Date: 2025/03/17
This script is designed to run RepeatMasker, and is based off of Sid's script and Daren's blog post.
It performs hard masking and then soft masking.
RepeatMasker can be found here:
https://www.repeatmasker.org/
https://github.com/Dfam-consortium/RepeatMasker

The aformentioned blog post can be found here:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
ScriptDescription

# Reference genome
reference_genome="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/from_Sekar/_completeAssemblies/Naja_nigricollis_najNig1/Assembly/najNig2.ragtag.scaffold_naNa.REHEADER.MT.fasta"

# Genome basename
species_name=$(basename "$reference_genome" .fasta)

# Set the RepeatMasker directory
repeat_masker_dir="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Naja_nigricollis/Results/2_RepeatMasker"

# Make the logs directory if it doesn't exist
[ ! -d "$repeat_masker_dir/Logs" ] && mkdir -p "$repeat_masker_dir/Logs"

# Set the script to create a log file for errors
log_file="$repeat_masker_dir/Logs/RepeatMasker_run_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$log_file") 2>&1

# Change the directory to the RepeatMasker directory if it isn't there already
cd "$repeat_masker_dir" || exit 1

# Print the current working directory
echo "Current working directory: $(pwd)"

# Exit immediately if a command exits with a non-zero status
set -e

# Function to report error
error_exit() {
    echo "Error occurred in script at line: $1"
    exit 1
}

# Trap errors and report the line number
trap 'error_exit $LINENO' ERR

# Create a set of masking rounds
round1="01_SSR_Mask"
round2="02_Squamate_BovB_CR1"
round3="03_RepBase_Tetrapoda_Mask"
round4="04_19Snake_Known_Mask"
round5="05_19Snake_Unknown_Mask"
round6="06_Full_Mask"

# Create a set of masking rounds
rounds=("$round1" "$round2" "$round3" \
        "$round4" "$round5" "$round6")

# Make the output directories if they don't already exist
for round in "${rounds[@]}"; do
    mkdir -p "$round"
done

# Set the number of threads for RepeatMasker -pa to use
t=40

# TODO: You need to check if this actually works
# Create a function to rename masked files to something more useful
rename_masked_file() {
    local directory="$1"   # Directory where the file is located
    local new_suffix="$2"  # New suffix to replace ".masked.fasta"

    # Find the masked genome file
    local masked_file
    masked_file=$(find "$directory" -type f -name "*.masked.fasta" | head -n 1)

    # Check if a file was found
    if [[ -z "$masked_file" ]]; then
        echo "Error: No masked genome file found in $directory"
        exit 1
    fi

    # Construct the new filename
    local new_masked_file="${masked_file/.masked.fasta/.$new_suffix.fasta}"

    # Rename the file
    mv "$masked_file" "$new_masked_file"

    # Return the new filename (optional)
    echo "$new_masked_file"
}

# Round 1: SSR masking
# Run simple sequence repeat (SSR) masking
# Note that Daren likes to do this first as it sames time to mask all of the simple repeats, 
# and it keeps simple, complex, and interspersed rpeats (e.g. TEs) seperate
echo -e "\e[31mRunning step 1: SSR masking\e[0m"
RepeatMasker \
	-pa "$t" \
	-a \
	-e ncbi \
	-dir "$round1" \
	-noint \
	-xsmall \
	"$reference_genome" 2>&1 | tee "Logs/$round1.log"

# Rename some of the files to make it more clear what they are
rename 's/fasta/simple_mask/g' "$round1"/*
rename 's/.masked$/.masked.fasta/g' "$round1"/*


# Set the repeat elements library for round 2
bovb_cr1="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/known_repeat_elements/CR1_BovB_Squamates_TElib.fasta"

# Find the output file from Round 1
round1_genome=$(find "$round1" -type f -name "*.masked.fasta" | head -n 1)

# Round 2: BovB/CR1 masking
# run BovB/CR1 masking
echo -e "\e[31mRunning step 2: BovB/CR1 masking\e[0m"
RepeatMasker \
	-pa "$t" \
	-engine ncbi \
	-lib "$bovb_cr1" \
	-a \
	-dir "$round2" \
	"$round1_genome" \
	-nolow 2>&1 | tee "Logs/$round2.log"

# Rename the files to something more useful
round2_genome=$(rename_masked_file "$round2" "BovB_masked")

# Round 3: RepBase Tetrapoda masking based on repease 20181026
# Run RepBase Tetrapoda masking based on repease 20181026
echo -e "\e[31mRunning step 3: RepBase Tetrapoda masking\e[0m"
RepeatMasker \
	-pa "$t" \
	-engine ncbi \
	-species tetrapoda \
	-a \
	-dir "$round3" \
	"$round2_genome" \
	-nolow 2>&1 | tee "Logs/$round3.log"

# Rename the file to something more useful
round3_genome=$(rename_masked_file "$round3" "Tetrapoda_masked")


# Round 4: "19snake" known masking
# Run "19snake" known masking
# Note that this comes from the Repclassifier output
# TODO: This is incomplete, but I need the complete repclassifier out for it
echo -e "\e[31mRunning step 4: 18 snake (Schield 2022) + croAtr2 known repeats\e[0m"
RepeatMasker -pa "$t" \
	-engine ncbi \
	-lib ./z_repeat_libraries/19Snakes.Known.clust.fasta \
	-a \
	-dir 04_19snake_known_mask \
	"$round3_genome" \
	-nolow 2>&1 | tee "Logs/$round5.log"

# Rename the file to something more useful
round4_genome=$(rename_masked_file "$round4" "04_19Snake_Known_Mask")

# Round 5: 05_19 Snake unknown masking
# run "19snake" unknown masking
echo -e "\e[31mRunning step 5: 18 snake (Schield 2022) + croAtr2 unknown repeats\e[0m"
RepeatMasker \
	-pa "$t" \
	-engine ncbi \
	-lib ./z_repeat_libraries/19Snakes.Unknown.clust.fasta \
	-a \
	-dir "$round5" \
	"$round4_genome" \
	-nolow 2>&1 | tee "Logs/$round5.log"

# Rename the file to something more useful
round5_genome=$(rename_masked_file "$round5" "04_19Snake_Known_Mask")


# Summarize/combine full output
cp "$round5_genome" "$round6/$species_name.complex_hard-mask.masked.fasta"


# .out
echo -e "\e[31mConcatenating .out outputs...\e[0m"
cat <(cat 01_SSR_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.simple_mask.out) \
	<(cat 02_squamate_bovb_cr1/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.bovb_mask.out | tail -n +4) \
	<(cat 03_repbase_tetrapoda_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.tetrapoda_mask.out | tail -n +4) \
	<(cat 04_19snake_known_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.snake_known_mask.out | tail -n +4) \
	<(cat 05_19snake_unknown_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.snake_unknown_mask.out | tail -n +4) > 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.out

# .cat
echo -e "\e[31mConcatenating .cat.gz outputs...\e[0m"
cat 01_SSR_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.simple_mask.cat.gz \
	02_squamate_bovb_cr1/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.bovb_mask.cat.gz \
	03_repbase_tetrapoda_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.tetrapoda_mask.cat.gz \
	04_19snake_known_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.snake_known_mask.cat.gz \
	05_19snake_unknown_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.snake_unknown_mask.cat.gz > 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.cat.gz

# .align
echo -e "\e[31mConcatenating .align outputs...\e[0m"
cat 01_SSR_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.simple_mask.align \
	02_squamate_bovb_cr1/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.bovb_mask.align \
	03_repbase_tetrapoda_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.tetrapoda_mask.align \
	04_19snake_known_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.snake_known_mask.align \
	05_19snake_unknown_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.snake_unknown_mask.align > 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.align

# .tbl
echo -e "\e[31mProcessing repeats...\e[0m"
ProcessRepeats -species tetrapoda 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.cat.gz

# create repeat landscape files
echo -e "\e[31mCreating repeat landscape files...\e[0m"
faToTwoBit /home/administrator/ExtraSSD2/Sid/Crotalus_atrox_genomics/z_Final_Annotation_Genome_Oct2024/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.fasta /home/administrator/ExtraSSD2/Sid/Crotalus_atrox_genomics/z_Final_Annotation_Genome_Oct2024/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.2bit
calcDivergenceFromAlign.pl -s 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.landscape 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.align
createRepeatLandscape.pl -div 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.landscape \
	-twoBit /home/administrator/ExtraSSD2/Sid/Crotalus_atrox_genomics/z_Final_Annotation_Genome_Oct2024/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.2bit \
	> 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.landscape.html

# create GFFs
echo -e "\e[31mCreating GFF3...\e[0m"
rmOutToGFF3.pl 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.out > 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.gff3

echo -e "\e[31mReformatting GFF3 (Daren Card method)...\e[0m"
cat 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.gff3 | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.reformat.gff3

# filter out simple repeats
# echo -e "\e[31mFiltering out simple repeats...\e[0m"
# grep -v -e "Satellite" -e ")n" -e "-rich" 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.reformat.gff3 > 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.reformat.complex.gff3

# soft-mask genome
echo -e "\e[31mSoft masking the reference genome...\e[0m"

# Find original Ns (come from scaffolding of the original pre-masked genome) and then find all Ns in the hard-masked genome
# Remove original Ns from the all N list
# Add to this the regions that were soft-masked in step 1
# Combine these sets of coordinates to get all masked regions and soft mask them

# TODO: Fix this to make it faster
# FIXME: This is slow as shit
### seqkit takes extremely long and should be optimized

seqkit locate --bed -rPp "N+" /home/administrator/ExtraSSD2/Sid/Crotalus_atrox_genomics/z_Final_Annotation_Genome_Oct2024/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.fasta > original_N_coords.bed
bedtools merge -i original_N_coords.bed > consolidated_original_N_coords.bed

seqkit locate --bed -rPp "N+" 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.complex_hard-mask.masked.fasta > masked_N_coords.bed
bedtools merge -i masked_N_coords.bed > consolidated_masked_N_coords.bed

bedtools subtract -a consolidated_masked_N_coords.bed -b consolidated_original_N_coords.bed > consolidated_new_masked_N_coords.bed

seqkit locate --bed -rPp "[acgtryswkmbdhvn]+" 01_SSR_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.simple_mask.masked.fasta > soft_masked_coords.bed
bedtools merge -i soft_masked_coords.bed > consolidated_soft_masked_coords.bed

cat consolidated_new_masked_N_coords.bed consolidated_soft_masked_coords.bed | sort -k1,1 -k2,2n > all_masked_coords.bed

bedtools maskfasta -soft -fi /home/administrator/ExtraSSD2/Sid/Crotalus_atrox_genomics/z_Final_Annotation_Genome_Oct2024/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.fasta \
	-bed all_masked_coords.bed \
	-fo 06_full_mask/croAtr2.ragtag.scaffold_cVir.cAda.REHEADER.MT.full_mask.soft.fasta
	
rm *.bed

# compress outputs
echo -e "\e[31mCompressing outputs...\e[0m"
gzip */*.out
gzip */*.fasta
gzip */*.align
