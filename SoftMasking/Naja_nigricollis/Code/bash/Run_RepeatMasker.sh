#!/bin/bash
: <<'ScriptDescription'
Date: 2025/03/17
This script is designed to run RepeatMasker, and is based off of Sid's script and Daren's blog post.
It performs hard masking and then soft masking.
RepeatMasker can be found here:
https://www.repeatmasker.org/
https://github.com/Dfam-consortium/RepeatMasker

The aforementioned blog post can be found here:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
ScriptDescription

# Reference genome
reference_genome="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/Sekar_Genomes/Scaffold_Assemblies/Elapidae/Naja/Naja_nigricollis_najNig1/Assembly/najNig2.ragtag.scaffold_naNa.REHEADER.MT.fasta"

# Set a directory for the genomic files create by this script
reference_genome_extra_dir="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/Sekar_Genomes/Scaffold_Assemblies/Elapidae/Naja/Naja_nigricollis_najNig1/Assembly/FromRepeatMaskerProcess"

# Make the above directory if it does not exist
[ ! -d "$reference_genome_extra_dir" ] && mkdir -p "$reference_genome_extra_dir"

# Genome basename
species_name=$(basename "$reference_genome" .fasta)

# Species name abreviation
species_abbreviation="najNig2"

# Set the RepeatMasker directory
repeat_masker_dir="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Naja_nigricollis/Results/2_RepeatMasker"

# Make the logs directory if it doesn't exist
[ ! -d "$repeat_masker_dir/Logs" ] && mkdir -p "$repeat_masker_dir/Logs"

# Set the script to create a log file for errors
log_file="$repeat_masker_dir/Logs/RepeatMasker_run_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$log_file") 2>&1

# Activate the mamba environment before doing anything else
echo "Activating mamba environment: RepeatMaskAnnot"
source /home/administrator/mambaforge/bin/activate RepeatMaskAnnot

# Check if the mamba environment is active
if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
	echo "No mamba environment is active. Exiting."
	exit 1
else
	echo "Mamba environment active: $CONDA_DEFAULT_ENV"
fi

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

# Function to check if the each round of RepeatMasker was successful
check_round() {
	local round_dir=$1

	find "$round_dir" -maxdepth 1 \( -name "*.fasta" -o -name "*.align" -o -name "*.cat.gz" -o -name "*.out" -o -name "*.tbl" \) | grep -q .
}

# Function to check if the each round of RepeatMasker was successful using the base run before files where renamed
check_round2() {
	local round_dir=$1

	find "$round_dir" -maxdepth 1 \( -name "*.fasta.masked" -o -name "*.fasta.align" -o -name "*.fasta.cat.gz" -o -name "*.fasta.out" -o -name "*.fasta.tbl" \) | grep -q .
}

# Function to check for required commands
check_requirements() {
	local missing_tools=()
	for cmd in RepeatMasker faToTwoBit twoBitInfo calcDivergenceFromAlign.pl createRepeatLandscape.pl rmOutToGFF3.pl seqkit bedtools; do
		if ! command -v "$cmd" &> /dev/null; then
			missing_tools+=("$cmd")
		fi
	done

	if [ ${#missing_tools[@]} -ne 0 ]; then
		echo "Error: The following required tools are missing:"
		for tool in "${missing_tools[@]}"; do
			echo "  - $tool"
		done
		echo "Please install these tools before running this script."
		exit 1
	fi
}

# Call this function to check for required commands
check_requirements



# Trap errors and report the line number
trap 'error_exit $LINENO' ERR

# Create a set of masking rounds
round1="01_SSR_Mask"
round2="02_Squamate_BovB_CR1"
round3="03_RepBase_Tetrapoda_mask"
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


# Round 1: SSR masking
# Run simple sequence repeat (SSR) masking
# Note that Daren likes to do this first as it sames time to mask all of the simple repeats, 
# and it keeps simple, complex, and interspersed rpeats (e.g. TEs) seperate
echo -e "\e[31mRunning step 1: SSR masking\e[0m"
if check_round "$round1"; then
	echo -e "\e[31Files found for $round1. Skipping.\e[0m"
else
	echo -e "\e[31Files not found for $round1. Running RepeatMasker now.\e[0m"

	# Run RepeatMasker
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
	rename 's/.masked$/.fasta/g' "$round1"/*

	# Check if the files were created, exit otherwise
	if check_round2 "$round1"; then
		echo "Files found for $round1. Continuing to the next round."
	else
		echo "Error: Files for $round1 were not created."
		exit 1
	fi
fi



# Set the repeat elements library for round 2
bovb_cr1="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/known_repeat_elements/CR1_BovB_Squamates_TElib.fasta"

# Find the output file from Round 1
round1_genome=$(find "$round1" -type f -name "*.simple_mask.fasta" | head -n 1)

# Round 2: BovB/CR1 masking
# run BovB/CR1 masking
echo -e "\e[31mRunning step 2: BovB/CR1 masking\e[0m"
if check_round "$round2"; then
	echo -e "\e[31Files found for $round2. Skipping.\e[0m"
else
	echo -e "\e[31Files not found for $round2. Running RepeatMasker now.\e[0m"

	# Run RepeatMasker
	RepeatMasker \
		-pa "$t" \
		-engine ncbi \
		-lib "$bovb_cr1" \
		-a \
		-dir "$round2" \
		"$round1_genome" \
		-nolow 2>&1 | tee "Logs/$round2.log"

	# Rename some of the files to make it more clear what they are
	rename 's/simple_mask.fasta/BovB_mask/g' "$round2"/*
	rename 's/.masked$/.fasta/g' "$round2"/*

	# Check if the files were created, exit otherwise
	if check_round2 "$round2"; then
		echo "Files found for $round2. Continuing to the next round."
	else
		echo "Error: Files for $round2 were not created."
		exit 1
	fi
fi



# Find the output file from Round 1
round2_genome=$(find "$round2" -type f -name "*.BovB_mask.fasta" | head -n 1)

# Round 3: RepBase Tetrapoda masking based on repease 20181026
# Run RepBase Tetrapoda masking based on repease 20181026
echo -e "\e[31mRunning step 3: RepBase Tetrapoda masking\e[0m"
if check_round "$round3"; then
	echo -e "\e[31Files found for $round3. Skipping.\e[0m"
else
	echo -e "\e[31Files not found for $round3. Running RepeatMasker now.\e[0m"

	# Run RepeatMasker
	RepeatMasker \
		-pa "$t" \
		-engine ncbi \
		-species tetrapoda \
		-a \
		-dir "$round3" \
		"$round2_genome" \
		-nolow 2>&1 | tee "Logs/$round3.log"

	# Rename some of the files to make it more clear what they are
	rename 's/BovB_mask.fasta/Tetrapoda_mask/g' "$round3"/*
	rename 's/.masked$/.fasta/g' "$round3"/*


	# Check if the files were created, exit otherwise
	if check_round2 "$round3"; then
		echo "Files found for $round3. Continuing to the next round."
	else
		echo "Error: Files for $round3 were not created."
		exit 1
	fi
fi


# Round 4: "19snake" known masking
# Find the FASTA from round 3
round3_genome=$(find "$round3" -type f -name "*.Tetrapoda_mask.fasta" | head -n 1)

# Set path for the 10th round of repclassifier
known_library="../1_Repclassifier/round-10_Self/round-10_Self.known"
# Define the custom known file for round 6 that Todd made containing 18 snake repeat elements
custom_18snake_known_repeats="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/known_repeat_elements/18Snakes.Known.clust.fasta"
# Set the new file path
new_known_19snakes_library="../1_Repclassifier/round-10_Self/19Snakes_with_$species_abbreviation.Known.clust.fasta"

# Concatenate the two files togehter
cat "$custom_18snake_known_repeats" "$known_library" > "$new_known_19snakes_library"

# Run "19snake" known masking
# Note that this comes from the Repclassifier output
echo -e "\e[31mRunning step 4: 18 snake (Schield 2022) + $species_abbreviation known repeats\e[0m"
if check_round "$round4"; then
	echo -e "\e[31Files found for $round4. Skipping.\e[0m"
else
	echo -e "\e[31Files not found for $round4. Running RepeatMasker now.\e[0m"

	# Run RepeatMasker
	RepeatMasker -pa "$t" \
		-engine ncbi \
		-lib "$new_known_19snakes_library" \
		-a \
		-dir "$round4" \
		"$round3_genome" \
		-nolow 2>&1 | tee "Logs/$round4.log"

	# Rename some of the files to make it more clear what they are
	rename 's/Tetrapoda_mask.fasta/19snake_known_mask/g' "$round4"/*
	rename 's/.masked$/.fasta/g' "$round4"/*


	# Check if the files were created, exit otherwise
	if check_round2 "$round4"; then
		echo "Files found for $round4. Continuing to the next round."
	else
		echo "Error: Files for $round4 were not created."
		exit 1
	fi
fi


# Round 5: 05_19 Snake unknown masking
# Find the FASTA file from round 4
round4_genome=$(find "$round4" -type f -name "*.19snake_known_mask.fasta" | head -n 1)

# Set path for the 10th round of repclassifier
unknown_library="../1_Repclassifier/round-10_Self/round-10_Self.unknown"
# Define the custom known file for round 6 that Todd made containing 18 snake repeat elements
custom_18snake_unknown_repeats="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/known_repeat_elements/18Snakes.Unknown.clust.fasta"
# Set the new file path
new_unknown_19snakes_library="../1_Repclassifier/round-10_Self/19Snakes_with_$species_abbreviation.Unknown.clust.fasta"

# Concatenate the two files togehter
cat "$custom_18snake_unknown_repeats" "$unknown_library" > "$new_unknown_19snakes_library"

# run "19snake" unknown masking
echo -e "\e[31mRunning step 5: 18 snake (Schield 2022) + $species_abbreviation unknown repeats\e[0m"
if check_round "$round5"; then
	echo -e "\e[31Files found for $round5. Skipping.\e[0m"
else
	echo -e "\e[31Files not found for $round5. Running RepeatMasker now.\e[0m"

	# Run RepeatMasker
	RepeatMasker \
		-pa "$t" \
		-engine ncbi \
		-lib "$new_unknown_19snakes_library" \
		-a \
		-dir "$round5" \
		"$round4_genome" \
		-nolow 2>&1 | tee "Logs/$round5.log"
	
	# Rename some of the files to make it more clear what they are
	rename 's/19snake_known_mask.fasta/19snake_unknown_mask/g' "$round5"/*
	rename 's/.masked$/.fasta/g' "$round5"/*


	# Check if the files were created, exit otherwise
	if check_round2 "$round5"; then
		echo "Files found for $round5. Continuing to the next round."
	else
		echo "Error: Files for $round5 were not created."
		exit 1
	fi
fi


# Find the FASTA file from round 4
round5_genome=$(find "$round5" -type f -name "*.19snake_unknown_mask.fasta" | head -n 1)

if check_round "$round6"; then
	echo -e "\e[31Files found for $round6. Skipping.\e[0m"
else
	echo -e "\e[31Files not found for $round6. Concatenating outputs now.\e[0m"

	# Summarize/combine full output
	cp "$round5_genome" "$round6/$species_name.complex_hard-mask.masked.fasta"

	# Concatenate all of the .out files from each of the runs
	echo -e "\e[31mConcatenating .out outputs...\e[0m"
	cat <(cat "$round1/$species_name.simple_mask.out") \
		<(cat "$round2/$species_name.BovB_mask.out" | tail -n +4) \
		<(cat "$round3/$species_name.Tetrapoda_mask.out" | tail -n +4) \
		<(cat "$round4/$species_name.19snake_known_mask.out" | tail -n +4) \
		<(cat "$round5/$species_name.19snake_unknown_mask.out" | tail -n +4) \
		> "$round6/$species_name.Full_Mask.out"

	# Combine RepeatMasker tabular files for all repeats - .out
	echo -e "\e[31mConcatenating .cat.gz outputs...\e[0m"
	cat "$round1/$species_name.simple_mask.cat.gz" \
		"$round2/$species_name.BovB_mask.cat.gz" \
		"$round3/$species_name.Tetrapoda_mask.cat.gz" \
		"$round4/$species_name.19snake_known_mask.cat.gz" \
		"$round5/$species_name.19snake_unknown_mask.cat.gz" \
		> "$round6/$species_name.Full_Mask.cat.gz"

	# Combine RepeatMasker repeat alignments for all repeats - .align
	echo -e "\e[31mConcatenating .align outputs...\e[0m"
	cat "$round1/$species_name.simple_mask.align" \
		"$round2/$species_name.BovB_mask.align" \
		"$round3/$species_name.Tetrapoda_mask.align" \
		"$round4/$species_name.19snake_known_mask.align" \
		"$round5/$species_name.19snake_unknown_mask.align" \
		> "$round6/$species_name.Full_Mask.align"

	# Resummarize repeat compositions from combined analysis of all RepeatMasker rounds
	echo -e "\e[31mProcessing repeats...\e[0m"
	ProcessRepeats \
		-species tetrapoda \
		"$round6/$species_name.Full_Mask.cat.gz"
fi


# Create repeat landscape files
echo -e "\e[31mCreating repeat landscape files...\e[0m"
# Convert the reference genome to .2bit
if [ ! -f "$reference_genome_extra_dir/$species_name.2bit" ]; then
	echo -e "\e[31mConverting reference genome to .2bit format...\e[0m"
	faToTwoBit "$reference_genome" "$reference_genome_extra_dir/$species_name.2bit"
else
	echo -e "\e[31mReference genome already in .2bit format. Skipping.\e[0m"
fi

# Calculate divergence
if [ ! -f "$round6/$species_name.Full_Mask.landscape" ]; then
	echo -e "\e[31mCalculating divergence...\e[0m"
	calcDivergenceFromAlign.pl -s "$round6/$species_name.Full_Mask.landscape" "$round6/$species_name.Full_Mask.align"
else
	echo -e "\e[31mDivergence already calculated. Skipping.\e[0m"
fi

# Create repeat landscape
if [ ! -f "$round6/$species_name.Full_Mask.landscape.html" ]; then
	echo -e "\e[31mCreating repeat landscape...\e[0m"
	createRepeatLandscape.pl -div "$round6/$species_name.Full_Mask.landscape" \
		-twoBit "$reference_genome_extra_dir/$species_name.2bit" \
		> "$round6/$species_name.Full_Mask.landscape.html"
else
	echo -e "\e[31mRepeat landscape already created. Checking if it has any data...\e[0m"
	# Check if the landscape file has data in it or not
	if ! grep -q "<svg" "$round6/$species_name.Full_Mask.landscape.html" || ! grep -q "RepeatLandscape" "$round6/$species_name.Full_Mask.landscape.html"; then
		echo "Warning: The landscape HTML file was created but appears to be empty or missing expected content. Recreating..."

		# Recreate the landscape file
		createRepeatLandscape.pl -div "$round6/$species_name.Full_Mask.landscape" \
			-twoBit "$reference_genome_extra_dir/$species_name.2bit" \
			> "$round6/$species_name.Full_Mask.landscape.html"
	else
		echo "Repeat landscape HTML file created successfully with content."
	fi
fi

# Create GFFs
if [ ! -f "$round6/$species_name.Full_Mask.gff3" ]; then
	echo -e "\e[31mCreating GFF3...\e[0m"
	rmOutToGFF3.pl "$round6/$species_name.Full_Mask.out" > "$round6/$species_name.Full_Mask.gff3"
else
	echo -e "\e[31mGFF3 already created. Skipping.\e[0m"
fi

# Reformat GFF3
if [ ! -f "$round6/$species_name.Full_Mask.reformat.gff3" ]; then
	echo -e "\e[31mReformatting GFF3 (Daren Card method)...\e[0m"
	cat "$round6/$species_name.Full_Mask.gff3" \
	| perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
	> "$round6/$species_name.Full_Mask.reformat.gff3"
else
	echo -e "\e[31mGFF3 already reformatted. Skipping.\e[0m"
fi

# Beggin soft-masking the genome
echo -e "\e[31mSoft masking the reference genome...\e[0m"

# Find original Ns (come from scaffolding of the original pre-masked genome) and then find all Ns in the hard-masked genome
# Remove original Ns from the all N list
# Add to this the regions that were soft-masked in step 1
# Combine these sets of coordinates to get all masked regions and soft mask them

# TODO: Fix this to make it faster
# OPTIMIZE: This is slow as shit apparently
# seqkit takes extremely long and should be optimized

# Take the reference genome and find all the hard masked features where Ns are
seqkit locate \
	--bed \
	-rPp "N+" "$reference_genome" \
	> original_N_coords.bed
# Merge the book end/overlapping features into a bed file, which come from the temporary bedfile above
bedtools merge -i original_N_coords.bed > consolidated_original_N_coords.bed

# Take the hard mask genome and find all of the hard masked features where Ns are
seqkit locate \
	--bed -rPp "N+" "$round6/$species_name.complex_hard-mask.masked.fasta" \
	> masked_N_coords.bed
# Merge the overlapping features into a single bed file
bedtools merge -i masked_N_coords.bed > consolidated_masked_N_coords.bed

# Subtract the features in the original genome from those in the masked genome to get features found only in the masked genome
bedtools subtract -a consolidated_masked_N_coords.bed -b consolidated_original_N_coords.bed > consolidated_new_masked_N_coords.bed

# Take the soft masked genome from the first round and create a bed file of those features
seqkit locate \
	--bed -rPp "[acgtryswkmbdhvn]+" "$round1_genome" \
	> soft_masked_coords.bed
# Merge the soft maskeed overlapping features
bedtools merge -i soft_masked_coords.bed > consolidated_soft_masked_coords.bed

# Concatenate the masked files together
cat consolidated_new_masked_N_coords.bed consolidated_soft_masked_coords.bed | sort -k1,1 -k2,2n > all_masked_coords.bed

# Soft mask the reference genome from the coordinates here
bedtools maskfasta -soft -fi "$reference_genome" \
	-bed all_masked_coords.bed \
	-fo "$round6/$species_name.Full_Mask.soft.fasta"

# Remove the temporary bed files
rm *.bed

# Compress the outputs to save space
echo -e "\e[31mCompressing outputs...\e[0m"
gzip */*.out
gzip */*.fasta
gzip */*.align
