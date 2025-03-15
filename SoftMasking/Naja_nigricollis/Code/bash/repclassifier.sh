#!/bin/bash

: <<'ScriptDescription'
Date: 2025/03/11
This script is designed to repclassifier, which is a piece of softwart that uses Repbase database and RepeatMasker to try and identify unknown elements with sequence similarity to curated repeat elements in Repbase.
It has another mode that can uses a custom library, in the form of a fasta file. This script will use both modes.
The program can be found here:
https://github.com/darencard/GenomeAnnotation/blob/master/repclassifier

The blog that explains the program can be found here:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
ScriptDescription

# Set the number of threads for repclassifier
threads=40

# Activate the mamba environment before doing anything else
source /home/administrator/mambaforge/bin/activate RepeatMaskAnnot

# Set the RepeatModeler directory
repeat_modeler_dir="/home/administrator/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Naja_nigricollis/Results/0_RepeatModeller"

# Set the repclassifier directory for output
repclassifier_dir="/home/administrator/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Naja_nigricollis/Results/1_Repclassifier"

# Set the RepeatModeler families file
repeat_modeler_families=$(find "$repeat_modeler_dir" -name "*-families.fa" -print -quit)

# Ensure the file was found
if [[ -z "$repeat_modeler_families" ]]; then
    echo "Error: No '*-families.fa' file found in $repeat_modeler_dir"
    exit 1
fi

# Set the prefix for headers in the RepeatModeler families file
species_prefix="najNig1_"

# Set the name of the new file with prefixes added
repeat_modeler_families_prefix=$(basename "$repeat_modeler_families" .fa).prefix.fa

# Step #1: Prepare the RepeatModeler families file for repclassifier by altering the file slightly
cat "$repeat_modeler_families" | seqkit fx2tab | \
	awk -v prefix="$species_prefix" '{ if ($1 ~ /^>/) print prefix $0; else print $0 }' | \
	seqkit tab2fx > "$repeat_modeler_dir/$repeat_modeler_families_prefix"

# Set the file name for unknown elements
unknown_elements_file="$repeat_modeler_dir/$repeat_modeler_families_prefix.unknown"

# Set the file name for known elements
known_elements_file="$repeat_modeler_dir/$repeat_modeler_families_prefix.known"

# Step #2: Split elements into known and unknown classification files
# Find the known elements
cat "$repeat_modeler_families_prefix" | seqkit fx2tab | \
	grep -v "Unknown" | seqkit tab2fx > "$known_elements_file"
# Find the unkown elements
cat "$repeat_modeler_families_prefix" | seqkit fx2tab | \
	grep "Unknown" | seqkit tab2fx > "$unknown_elements_file"

# Make sure both files were created
# Ensure the known file was found
if [[ -z "$known_elements_file" ]]; then
    echo "Error: No '*-families.fa.known' file found in $repeat_modeler_dir"
    exit 1
fi
# Ensure the unknown file was found
if [[ -z "$unknown_elements_file" ]]; then
    echo "Error: No '*-families.fa.unknown' file found in $repeat_modeler_dir"
    exit 1
fi

# Step #3-12: Iterative repclassifier rounds
# Define the number of rounds to run the loop for
num_rounds=10
# Define the custom known file for round 6 that Todd made containing 18 snake repeat elements
custom_snake_repeats="/home/administrator/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/known_repeat_elements/18Snakes.Known.clust.fasta"

# Loop for all rounds (1-10)
for ((round=1; round <= num_rounds; round++)); do
	# If round 1, use the repclassifier with these settings
	if ((round==1)); then
		# Set the previous round directory
		prev_round_dir="$repeat_modeler_dir"
		# Set the current round output directory
		current_round_dir="$repclassifier_dir/round-01_RepbaseTetrapoda_Self"
		# Set the unknown elements file for this round
		unknown_file="$unknown_elements_file"
		# Set the known elements file for this round
		known_file="$known_elements_file"

		# Run repclassifier
		repclassifier \
			-t "$threads" \
			-d Tetrapoda \
			-u "$unknown_file" \
			-k "$known_file" \
			-a "$known_file" \
			-o "$current_round_dir"
	
	# If the round is between 2 and 5, use repclassifier with these settings
	elif ((2<=round<=5)); then
		# Special case for round 2, which uses output from round 1 with its special name
		if ((round==2)); then
			prev_round_dir="$repclassifier_dir/round-01_RepbaseTetrapoda_Self"
			unknown_file="$prev_round_dir/round-01_RepbaseTetrapoda_Self.unknown"
			known_file="$prev_round_dir/round-01_RepbaseTetrapoda_Self.known"
		else
			# For rounds 3-5, use the standard naming pattern
			prev_round_dir="$repclassifier_dir/round-$(printf "%02d" $((round-1)))_Self"
			unknown_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.unknown"
			known_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.known"
		fi

		# Set the current round output directory (same format for rounds 2-5)
		current_round_dir="$repclassifier_dir/round-$(printf "%02d" "$round")_Self"

		# Run repclassifier
		repclassifier \
			-t "$threads" \
			-u "$unknown_file" \
			-k "$known_file" \
			-a "$known_file" \
			-o "$current_round_dir"

	# If the round is 6, use repclassifier with the snake repeat elements file that Todd made long ago
	elif ((round == 6)); then
		# Set the previous round directory
		prev_round_dir="$repclassifier_dir/round-$(printf "%02d" $((round-1)))_Self"
		# Set the current round directory
		current_round_dir="$repclassifier_dir/round-06_18Snakes"
		# Set the unknown elements file
		unknown_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.unknown"
		# Set the file for the known elements
		known_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.known"

		# Run repclassifier
		repclassifier \
			-t "$threads" \
			-u "$unknown_file" \
			-k "$custom_snake_repeats" \
			-a "$known_file" \
			-o "$current_round_dir"

	# All other rounds
    else
		# Special case for round 7, as that round used the custom snake fasta library, and has a different name
		if ((round==7)); then
			# Set the previous round to round 6
			prev_round_dir="$repclassifier_dir/round-06_18Snakes"
			# Set the unknown file to the round 6 unknown elements
			unknown_file="$prev_round_dir/round-06_18Snakes.unknown"
			# Set the known file to the round 6 known elements
			known_file="$prev_round_dir/round-06_18Snakes.known"
		else
			# For rounds 8 and up
			prev_round_dir="$repclassifier_dir/round-$(printf "%02d" $((round-1)))_Self"
			unknown_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.unknown"
			known_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.known"
		fi

		# Set the current round output directory (same format for rounds 8+)
		current_round_dir="$repclassifier_dir/round-$(printf "%02d" "$round")_Self"

		# Run repclassifier
		repclassifier \
			-t "$threads" \
			-u "$unknown_file" \
			-k "$known_file" \
			-a "$known_file" \
			-o "$current_round_dir"
    fi

    echo "Completed round $round"
done

echo "All repclassifier rounds completed."