#!/bin/bash
: <<'ScriptDescription'
Date: 2025/06/12
This script is designed to repclassifier, which is a piece of softwart that uses Repbase database and RepeatMasker to try and identify unknown elements with sequence similarity to curated repeat elements in Repbase.
It has another mode that can uses a custom library, in the form of a fasta file. This script will use both modes.
The program can be found here:
https://github.com/darencard/GenomeAnnotation/blob/master/repclassifier

The blog that explains the program can be found here:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/
ScriptDescription

# Set the start time of the script
start_time=$(date '+%Y-%m-%d %H:%M:%S')

# Set ntfy.sh topic for notifications
ntfy_topic="kaas-ballard-Robin-scripts-72724027625978"

# NOTE: Change this each time you run the script
# Set the RepeatModeler directory
repeat_modeler_dir="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Elapidae/Hydrophis/Hydrophis_curtus/Results/0_RepeatModeler"

# NOTE: Change this each time you run the script
# Set the repclassifier directory for output
repclassifier_dir="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Elapidae/Hydrophis/Hydrophis_curtus/Results/1_Repclassifier"

# Create log directory under the output directory if it does not exist
[ ! -d "$repclassifier_dir/Logs" ] && mkdir -p "$repclassifier_dir/Logs"

# Set up logging
log_file="$repclassifier_dir/Logs/repclassifier_run_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$log_file") 2>&1

echo "Starting repclassifier script at $(date)"
echo "Log file created at: $log_file"

# NOTE: Change this when needed
# Set the number of threads for repclassifier
threads=40

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

# NOTE: Change this each time you run the script
# Define the path to repclassifier
repclassifier_cmd="$HOME/ExtraSSD2/Sid/Crotalus_atrox_genomics/9_NewRepeatAnnotationAndMasking/InstallationFiles/repclassifier"

# Set the RepeatModeler families file
repeat_modeler_families=$(find "$repeat_modeler_dir" -name "*-families.fa" -print -quit)
echo "RepeatModeler families file: $repeat_modeler_families"

# Ensure the file was found
if [[ -z "$repeat_modeler_families" ]]; then
    echo "Error: No '*-families.fa' file found in $repeat_modeler_dir"
    exit 1
fi

# NOTE: Change this each time you run the script
species_name="Hydrophis_curtus"

# NOTE: Change this each time you run the script
# Set the prefix for headers in the RepeatModeler families file
species_prefix="hyoCur1_"

# Send myself a notification that the script is starting
curl -d "🔔 Starting repclassifier for $species_name at $(date). Check logs at $repclassifier_dir/Logs/ for details." \
	ntfy.sh/"$ntfy_topic"

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
cat "$repeat_modeler_dir/$repeat_modeler_families_prefix" | seqkit fx2tab | \
	grep -v "Unknown" | seqkit tab2fx > "$known_elements_file"
# Find the unkown elements
cat "$repeat_modeler_dir/$repeat_modeler_families_prefix" | seqkit fx2tab | \
	grep "Unknown" | seqkit tab2fx > "$unknown_elements_file"

# Make sure both files were created
# Ensure the known file was found - CORRECTED ERROR MESSAGE
if [[ ! -f "$known_elements_file" ]]; then
	echo "Error: Known elements file not created: $known_elements_file"
	exit 1
fi
# Ensure the unknown file was found - CORRECTED ERROR MESSAGE
if [[ ! -f "$unknown_elements_file" ]]; then
	echo "Error: Unknown elements file not created: $unknown_elements_file"
	exit 1
fi

# Step #3-12: Iterative repclassifier rounds
# Define the number of rounds to run the loop for
num_rounds=10

# NOTE: Change this when needed
# Define the custom known file for round 6 that Todd made containing 18 snake repeat elements
custom_snake_repeats="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/known_repeat_elements/18Snakes.Known.clust.fasta"
echo "Custom snake repeats file: $custom_snake_repeats"

# Move to the repclassifier output directory
cd "$repclassifier_dir" || return

# Print the working directory to make sure I am in the the right one
echo "Current working directory: $(pwd)"

# Loop for all rounds (1-10)
for ((round=1; round<=num_rounds; round++)); do
	echo "-------------------------------------------------------"
	echo "Starting round $round at $(date)"

	# If round 1, use the repclassifier with these settings
	if ((round==1)); then
		# Set the previous round directory
		prev_round_dir="$repeat_modeler_dir"
		# Set the current round output directory
		current_round_dir="round-01_RepbaseTetrapoda_Self"
		# Set the unknown elements file for this round
		unknown_file="$unknown_elements_file"
		# Set the known elements file for this round
		known_file="$known_elements_file"

			echo "Round 1: Using Repbase Tetrapoda database"
			echo "Previous round directory: $prev_round_dir"
			echo "Current round directory: $current_round_dir"
			echo "Unknown file: $unknown_file"
			echo "Known file: $known_file"

		# Run repclassifier
		echo "Running repclassifier for round 1..."
		"$repclassifier_cmd" \
			-t "$threads" \
			-d Tetrapoda \
			-u "$unknown_file" \
			-k "$known_file" \
			-a "$known_file" \
			-o "$current_round_dir"
	
	# If the round is between 2 and 5, use repclassifier with these settings
	elif ((2<=round && round<=5)); then
		# Special case for round 2, which uses output from round 1 with its special name
		if ((round==2)); then
			prev_round_dir="round-01_RepbaseTetrapoda_Self"
			unknown_file="$prev_round_dir/round-01_RepbaseTetrapoda_Self.unknown"
			known_file="$prev_round_dir/round-01_RepbaseTetrapoda_Self.known"
				echo "Round 2: Using output from round 1 (special case)"
		else
			# For rounds 3-5, use the standard naming pattern
			prev_round_dir="round-$(printf "%02d" $((round-1)))_Self"
			unknown_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.unknown"
			known_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.known"
				echo "Round $round: Using standard naming pattern"
		fi

		# Set the current round output directory (same format for rounds 2-5)
		current_round_dir="round-$(printf "%02d" "$round")_Self"
			
			echo "Previous round directory: $prev_round_dir"
			echo "Current round directory: $current_round_dir"
			echo "Unknown file: $unknown_file"
			echo "Known file: $known_file"

		# Run repclassifier
		"$repclassifier_cmd" \
			-t "$threads" \
			-u "$unknown_file" \
			-k "$known_file" \
			-a "$known_file" \
			-o "$current_round_dir"

	# If the round is 6, use repclassifier with the snake repeat elements file that Todd made long ago
	elif ((round==6)); then
		# Set the previous round directory
		prev_round_dir="round-$(printf "%02d" $((round-1)))_Self"
		# Set the current round directory
		current_round_dir="round-06_18Snakes"
		# Set the unknown elements file
		unknown_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.unknown"
		# Set the file for the known elements
		known_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.known"

			echo "Round 6: Using custom snake repeats library"
			echo "Previous round directory: $prev_round_dir"
			echo "Current round directory: $current_round_dir"
			echo "Unknown file: $unknown_file"
			echo "Known file: $known_file"
			echo "Custom snake repeats: $custom_snake_repeats"

		# Run repclassifier
		echo "Running repclassifier for round 6..."
		"$repclassifier_cmd" \
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
			prev_round_dir="round-06_18Snakes"
			# Set the unknown file to the round 6 unknown elements
			unknown_file="$prev_round_dir/round-06_18Snakes.unknown"
			# Set the known file to the round 6 known elements
			known_file="$prev_round_dir/round-06_18Snakes.known"
				echo "Round 7: Using output from round 6 (special case)"
		else
			# For rounds 8 and up
			prev_round_dir="round-$(printf "%02d" $((round-1)))_Self"
			unknown_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.unknown"
			known_file="$prev_round_dir/round-$(printf "%02d" $((round-1)))_Self.known"
				echo "Round $round: Using standard naming pattern"
		fi

		# Set the current round output directory (same format for rounds 8+)
		current_round_dir="round-$(printf "%02d" "$round")_Self"

			echo "Previous round directory: $prev_round_dir"
			echo "Current round directory: $current_round_dir"
			echo "Unknown file: $unknown_file"
			echo "Known file: $known_file"

		# Run repclassifier
		echo "Running repclassifier for round $round..."
		"$repclassifier_cmd" \
			-t "$threads" \
			-u "$unknown_file" \
			-k "$known_file" \
			-a "$known_file" \
			-o "$current_round_dir"
	fi

	echo "Completed round $round at $(date)"

	# Check if the round completed successfully
	if [[ $? -eq 0 ]]; then
		echo "Round $round completed successfully"
	else
		echo "Error: Round $round failed with exit code $?"
		exit 1
	fi
done

echo "All repclassifier rounds completed successfully at $(date)"

# Set the end time of the script
end_time=$(date '+%Y-%m-%d %H:%M:%S')

# Calculate the total runtime
runtime=$(date -u -d "$end_time" +"%s")-$(( $(date -u -d "$start_time" +"%s") ))

# Tell the user in the log file and term how long the script took to run
echo "Total runtime: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) seconds"

# Send a notification that the script has finished successfully or unsuccessfully
# Check if the final known and unknown files exist
if [ -f "$repclassifier_dir/round-10_Self/round-10_Self.known" ] && [ -f "$repclassifier_dir/round-10_Self/round-10_Self.unknown" ]; then
	curl -d "✅ SUCCESS: repclassifier completed for $species_name at $(date). Total runtime was: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) seconds. Check logs at $repclassifier_dir/Logs/" \
		ntfy.sh/"$ntfy_topic"
# If the files do not exist, send a failure notification
else
	curl -d "❌ FAILED: repclassifier failed for $species_name at $(date). Check logs for errors." \
		ntfy.sh/"$ntfy_topic"
fi