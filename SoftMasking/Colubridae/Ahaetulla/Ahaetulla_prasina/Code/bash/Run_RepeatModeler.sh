#!/bin/bash

: <<'ScriptDescription'
Date: 2025/06/16

This script is designed to run RepeatModeler. I am following the method found in Daren Card's blog post:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

ScriptDescription

# Activate the mamba environment containing RepeatModeler
source "$HOME/mambaforge/bin/activate" RepeatMaskAnnot

# Set the start time of the script
start_time=$(date '+%Y-%m-%d %H:%M:%S')

# Set ntfy.sh topic for notifications
ntfy_topic="kaas-ballard-Robin-scripts-72724027625978"

# NOTE: Change this when needed
# Set the number of threads for RepeatModeler
threads=10

# SET: Change this each time you run this script
# Set the name for the RepeatModeler database
species_name="Ahaetulla_prasina"

# SET: Change this each time you run this script
# Set the family name for the species
family="Colubridae"

# SET: Change this each time you run this script
# Set the genus name for the species
genus="Ahaetulla"

# SET: Change this each time you run this script
# Reference genome
reference_genome="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/Public_Genomes/Colubridae/Ahaetulla/Ahaetulla_prasina/ncbi_dataset/data/GCF_028640845.1/GCF_028640845.1_ASM2864084v1_genomic.fna"

# NOTE: Change this when needed
# Set the output directory
output_directory="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/$family/$genus/$species_name/Results/0_RepeatModeler"

# Create log directory under the output directory if it does not exist
[ ! -d "$output_directory/Logs" ] && mkdir -p "$output_directory/Logs"

# Send myself a notification that the script is starting
curl -d "🔔 Starting RepeatModeler for $species_name at $(date). Reference genome: $reference_genome. Check logs at $output_directory/Logs/ for details." \
	ntfy.sh/"$ntfy_topic"

# Change the directory to the output directory so that the RM_ files are created there
cd "$output_directory" || { echo "Failed to change directory to $output_directory"; exit 1; }

# Step #1: Build a new RepeatModeler database for the reference genome
# Note: BuildDatabase can build directories that don't exist yet
if [ ! -f "$output_directory/$species_name.nsq" ]; then
	echo "Warning: Database at the specified path does not exist. Building a new database now."
	BuildDatabase -name "$output_directory/$species_name" \
		"$reference_genome" 2>&1 | tee "$output_directory/Logs/BuildDatabase.log"
else
	echo "Database already exists at the specified path. Skipping."
fi

# Step #2: Detect if the database was created and then run RepeatModeler
if [ ! -f "$output_directory/$species_name.nsq" ]; then
	echo "Error: Database was not created. Exiting."
	exit 1
else
	# Run RepeatModeler
	RepeatModeler \
		-threads "$threads" \
		-database "$output_directory/$species_name" 2>&1 | tee "$output_directory/Logs/RepeatModeler.log"
fi

# Set the end time of the script
end_time=$(date '+%Y-%m-%d %H:%M:%S')

# Calculate the duration of the script
runtime=$(date -u -d "$end_time" +"%s")-$(( $(date -u -d "$start_time" +"%s") ))

# Tell the user in the log file and term how long the script took to run
echo "Total runtime: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) seconds"

# Send a notification that the script has finished
if [ $? -eq 0 ] && \
	[ -f "$output_directory/$species_name.nsq" ] && \
	[ -f "$output_directory/$species_name-families.fa" ] && \
	[ -f "$output_directory/$species_name-families.stk" ] && \
	[ -f "$output_directory/$species_name.njs" ]; then # if [ $? -eq 0 ] checks if the last command was successful
	curl -d "✅ SUCCESS: RepeatModeler completed for $species_name at $(date). Total runtime: $((runtime / 3600)) hours, $(((runtime % 3600) / 60)) minutes, $((runtime % 60)) secondsCheck logs at $output_directory/Logs/" \
		ntfy.sh/"$ntfy_topic"
else
	curl -d "❌ FAILED: RepeatModeler failed for $species_name at $(date). Check logs for errors." \
		ntfy.sh/"$ntfy_topic"
fi