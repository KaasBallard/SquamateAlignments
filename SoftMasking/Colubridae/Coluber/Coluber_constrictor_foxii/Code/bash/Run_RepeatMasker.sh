#!/bin/bash
: <<'ScriptDescription'
Date: 2025/06/12
This script is designed to run RepeatMasker, and is based off of Sid's script and Daren's blog post.
It performs hard masking and then soft masking.
RepeatMasker can be found here:
https://www.repeatmasker.org/
https://github.com/Dfam-consortium/RepeatMasker

The aforementioned blog post can be found here:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
ScriptDescription

# ====================  SETUP ====================
# Set the start time of the script
start_time=$(date '+%Y-%m-%d %H:%M:%S')

# Set ntfy.sh topic for notifications
ntfy_topic="kaas-ballard-Robin-scripts-72724027625978"

# SET: Change this each time you run this script
# Set the family name for the species
family="Colubridae"

# SET: Change this each time you run this script
# Set the genus name for the species
genus="Coluber"

# SET: Change this each time you run this script
# Set the name for the species
species="Coluber_constrictor_foxii"

# SET: Change this when needed
# Species name abreviation
species_abbreviation="colConFoxi1"

# SET: Change this when needed
# Reference genome
reference_genome="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/Public_Genomes/Colubridae/Coluber/Coluber_constrictor_foxii/ncbi_dataset/data/GCA_038048745.1/GCA_038048745.1_rColCon1_p1.0_genomic.fna"

# NOTE: Change this when needed
# Set a directory for the genomic files create by this script
reference_genome_extra_dir="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/Public_Genomes/$family/$genus/$species/FromRepeatMaskerProcess"

# Make the above directory if it does not exist
[ ! -d "$reference_genome_extra_dir" ] && mkdir -p "$reference_genome_extra_dir"

# Genome basename
species_name=$(basename "$reference_genome" | sed -E 's/\.(fasta|fna|fa)$//') # This will allow for more filename extensions in the species name

# NOTE: Change this when needed
# Set the RepeatMasker directory
repeat_masker_dir="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/$family/$genus/$species/Results/2_RepeatMasker"

# Make the logs directory if it doesn't exist
[ ! -d "$repeat_masker_dir/Logs" ] && mkdir -p "$repeat_masker_dir/Logs"

# Set the script to create a log file for errors
log_file="$repeat_masker_dir/Logs/RepeatMasker_run_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$log_file") 2>&1

# NOTE: Change this when needed
# Add the FindRegionCoordinates program to the PATH
export PATH="$HOME/ExtraSSD2/Kaas/Projects/FindRegionCoordinates:$PATH"

# NOTE: Change this when needed
# Add the PERL5 path
export PERL5LIB=$HOME/mambaforge/envs/RepeatMaskAnnot/share/RepeatMasker:$PERL5LIB

# Change the directory to the RepeatMasker directory if it isn't there already
cd "$repeat_masker_dir" || exit 1

# Print the current working directory
echo "Current working directory: $(pwd)"

# Exit immediately if a command exits with a non-zero status
set -e

# ==================== FUNCTIONS ====================
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

# NOTE: May need to change this based on the environment set up of the system.
# Function to activate the correct environment once
activate_environment() {
    local env_name="$1"
    echo "Activating environment: $env_name"
    
    if [[ "$env_name" == "RepeatMaskAnnot" ]]; then
        source /home/administrator/mambaforge/bin/activate RepeatMaskAnnot
    elif [[ "$env_name" == "FindRegionCoordinates" ]]; then
        source "$HOME/ExtraSSD2/Kaas/Projects/FindRegionCoordinates/.venv/bin/activate"
    else
        echo "Unknown environment: $env_name"
        return 1
    fi
    
    # Verify environment activation
    if [[ "$env_name" == "RepeatMaskAnnot" && -z "$CONDA_DEFAULT_ENV" ]]; then
        echo "Failed to activate $CONDA_DEFAULT_ENV environment. Exiting."
        exit 1
    fi
    
    echo "Environment active: $CONDA_DEFAULT_ENV"
    return 0
}




# ==================== MAKE REPEAT MASKER ROUNDS ====================
# Activate the mamba environment
activate_environment "RepeatMaskAnnot"

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

# NOTE: Change this when needed
# Set the number of threads for RepeatMasker -pa to use
t=10

# ==================== REPEAT MASKER ROUND 1 ====================
# Round 1: SSR masking
# Run simple sequence repeat (SSR) masking
# Note that Daren likes to do this first as it sames time to mask all of the simple repeats, 
# and it keeps simple, complex, and interspersed rpeats (e.g. TEs) seperate
echo -e "\e[31mRunning step 1: SSR masking\e[0m"
if check_round "$round1"; then
	echo -e "\e[31Files found for $round1. Skipping.\e[0m"
else
	echo -e "\e[31Files not found for $round1. Running RepeatMasker now.\e[0m"

	# Send a notifcation that the first round is starting
	curl -d "🔔 Starting RepeatMasker round 1: SSR masking for $reference_genome at $(date). Check logs at $repeat_masker_dir/Logs/ for details." \
		ntfy.sh/"$ntfy_topic"

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
	rename 's/\.(fasta|fna)\./.simple_mask./g' "$round1"/*
	rename 's/.masked$/.fasta/g' "$round1"/*

	# Check if the files were created, exit otherwise
	if check_round "$round1"; then
		echo "Files found for $round1. Continuing to the next round."
	else
		echo "Error: Files for $round1 were not created."
		exit 1
	fi
fi


# ==================== REPEAT MASKER ROUND 2 ====================

# NOTE: Change this when needed
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
	if check_round "$round2"; then
		echo "Files found for $round2. Continuing to the next round."
	else
		echo "Error: Files for $round2 were not created."
		exit 1
	fi
fi



# ==================== REPEAT MASKER ROUND 3 ====================
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
	if check_round "$round3"; then
		echo "Files found for $round3. Continuing to the next round."
	else
		echo "Error: Files for $round3 were not created."
		exit 1
	fi
fi


# ==================== REPEAT MASKER ROUND 4 ====================
# Round 4: "19snake" known masking
# Find the FASTA from round 3
round3_genome=$(find "$round3" -type f -name "*.Tetrapoda_mask.fasta" | head -n 1)

# Set path for the 10th round of repclassifier
known_library="../1_Repclassifier/round-10_Self/round-10_Self.known"

# NOTE: Change these when needed
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
	if check_round "$round4"; then
		echo "Files found for $round4. Continuing to the next round."
	else
		echo "Error: Files for $round4 were not created."
		exit 1
	fi
fi


# ==================== REPEAT MASKER ROUND 5 ====================
# Round 5: 05_19 Snake unknown masking
# Find the FASTA file from round 4
round4_genome=$(find "$round4" -type f -name "*.19snake_known_mask.fasta" | head -n 1)

# Set path for the 10th round of repclassifier
unknown_library="../1_Repclassifier/round-10_Self/round-10_Self.unknown"
# Define the custom known file for round 6 that Todd made containing 18 snake repeat elements

# NOTE: Change these when needed
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
	if check_round "$round5"; then
		echo "Files found for $round5. Continuing to the next round."
	else
		echo "Error: Files for $round5 were not created."
		exit 1
	fi
fi


# ==================== REPEAT MASKER ROUND 6 ====================
# Find the FASTA file from round 5
round5_genome=$(find "$round5" -type f -name "*.19snake_unknown_mask.fasta" | head -n 1)

if check_round "$round6"; then
	echo -e "\e[31Files found for $round6. Skipping.\e[0m"
else
	# Send a notifcation that the final round is starting
	curl -d "🔔 Starting RepeatMasker round 6: Full Mask for $reference_genome at $(date). Check logs at $repeat_masker_dir/Logs/ for details." \
		ntfy.sh/"$ntfy_topic"

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


# ==================== POST-PROCESSING ====================
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

# ==================== SOFT MASKING ====================
# Beggin soft-masking the genome
echo -e "\e[31mSoft masking the reference genome...\e[0m"

# Activate the virtual environment
activate_environment "FindRegionCoordinates"

# Set path for the soft masked BED files
masked_coords_dir="Masked_Coordinates"

# Make a directory for the BED files
[ ! -d "$masked_coords_dir" ] && mkdir -p "$masked_coords_dir"

# Set the hard masked coordinates file name
original_hard_masked_coords_bed="$masked_coords_dir/$species_abbreviation.original_hard_masked_coordinates.bed"

# Find all of the hard masked features already in the genome from assembly
if [ -f "$original_hard_masked_coords_bed" ]; then
	echo -e "\e[31mHard masked coordinates file already exists. Skipping..\e[0m"
else
	echo -e "\e[31mHard masked coordinates file does not exist. Creating...\e[0m"
	# Take the reference genome and find all the hard masked features where Ns are
	FindRegionCoordinates.py \
		-f "$reference_genome" \
		-b "$original_hard_masked_coords_bed"
fi


# Set the hard masked coordinates file name for round 6
new_hard_masked_coords_bed="$masked_coords_dir/$species_abbreviation.new_hard_masked_coordinates.bed"

# Find all of the hard masked features found in the genome at round 6
if [ -f "$new_hard_masked_coords_bed" ]; then
	echo -e "\e[31mHard masked coordinates file already exists. Skipping..\e[0m"
else
	echo -e "\e[31mHard masked coordinates file does not exist. Creating...\e[0m"
	# Take the hard masked genome from round 6 and find all the hard masked features where Ns are
	FindRegionCoordinates.py \
		-f "$round6/$species_name.complex_hard-mask.masked.fasta" \
		-b "$new_hard_masked_coords_bed"
fi



# Activate the mamba environment
activate_environment "RepeatMaskAnnot"

# Set the hard masked coordinates file name
new_minus_original_coords_bed="$masked_coords_dir/$species_abbreviation.new_minus_original_coordinates.bed"

# Subtract the features in the original genome from those in the masked genome to get features found only in the masked genome
bedtools subtract -a "$new_hard_masked_coords_bed" -b "$original_hard_masked_coords_bed" > "$new_minus_original_coords_bed"
# bedtools subtract the hard masked N's plus scaffolding N's and subtract the scaffolding N's from the hard masked N's, leaving only the hard masked N's


# Activate the virtual environment
activate_environment "FindRegionCoordinates"

# Set the soft masked coordinates file name
soft_masked_coords_bed="$masked_coords_dir/$species_abbreviation.soft_masked_coordinates.bed"

# Find all of the soft masked features found in the genome at round 1
if [ -f "$soft_masked_coords_bed" ]; then
	echo -e "\e[31mSoft masked coordinates file already exists. Skipping..\e[0m"
else
	echo -e "\e[31mSoft masked coordinates file does not exist. Creating...\e[0m"
	# Take the soft masked genome from the first round and create a bed file of those features
	FindRegionCoordinates.py \
		-f "$round1_genome" \
		-r "[acgtryswkmbdhvn]+" \
		-b "$soft_masked_coords_bed"
fi

# Set the all masked coordinates file name
all_masked_coords_bed="$masked_coords_dir/$species_abbreviation.all_masked_coordinates.bed"

# Concatenate the masked files togther
cat "$new_minus_original_coords_bed" "$soft_masked_coords_bed" | sort -k1,1 -k2,2n > "$all_masked_coords_bed"

# Activate the mamba environment
activate_environment "RepeatMaskAnnot"

# Soft mask the reference genome from the coordinates here
bedtools maskfasta -soft -fi "$reference_genome" \
	-bed "$all_masked_coords_bed" \
	-fo "$round6/$species_name.Full_Mask.soft.fasta"

# Compress the outputs to save space
echo -e "\e[31mCompressing outputs...\e[0m"
gzip */*.out
gzip */*.fasta
gzip */*.align
gzip */*.gff3
gzip */*.tbl

# Set an end time for the script
end_time=$(date '+%Y-%m-%d %H:%M:%S')

# Calculate the duration of the script
runtime=$(date -u -d "$end_time" +"%s")-$(( $(date -u -d "$start_time" +"%s") ))

# Tell the user in the log file and term how long the script took to run
echo "Total runtime: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) seconds"

# Send a notification that the script has finished
# Check if the last gzip command was successful and if key output files exist
if [ $? -eq 0 ] && \
	[ -f "$round6/$species_name.Full_Mask.soft.fasta.gz" ] && \
	[ -f "$round6/$species_name.Full_Mask.out.gz" ] && \
	[ -f "$round6/$species_name.Full_Mask.cat.gz" ] && \
	[ -f "$round6/$species_name.Full_Mask.landscape.html" ] && \
	[ -f "$round6/$species_name.Full_Mask.reformat.gff3.gz" ] && \
	[ -f "$round6/$species_name.Full_Mask.gff3.gz" ] && \
	[ -f "$round6/$species_name.Full_Mask.align.gz" ] && \
	[ -f "$round6/$species_name.Full_Mask.tbl.gz" ]; then
	curl -d "✅ SUCCESS: RepeatMasker completed for $species_abbreviation at $(date). Total runtime was: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) seconds. Results in $repeat_masker_dir/$round6/" \
		ntfy.sh/"$ntfy_topic"
else
	curl -d "❌ FAILED: RepeatMasker failed or key output files are missing for $species_abbreviation at $(date). Check logs at $repeat_masker_dir/Logs/ for errors." \
		ntfy.sh/"$ntfy_topic"
fi