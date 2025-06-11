#!/bin/bash

: <<'ScriptDescription'
Date: 2025/06/11

This script is designed to run RepeatModeler. I am following the method found in Daren Card's blog post:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

ScriptDescription

# Activate the mamba environment containing RepeatModeler
source "$HOME/mambaforge/bin/activate" RepeatMaskAnnot

# NOTE: Change this when needed
# Set the number of threads for RepeatModeler
threads=44

# NOTE: Change this each time you run this script
# Set the output directory
output_directory="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Elapidae/Hydrophis/Hydrophis_curtus/Results/0_RepeatModeler"

# Create log directory under the output directory if it does not exist
[ ! -d "$output_directory/Logs" ] && mkdir -p "$output_directory/Logs"

# NOTE: Change this each time you run this script
# Set the name for the RepeatModeler database
database_name="Hydrophis_curtus"

# NOTE: Change this each time you run this script
# Reference genome
reference_genome="$HOME/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/Public_Genomes/Elapidae/Hydrophis/Hydrophis_curtus/ncbi_dataset/data/GCA_037043045.1/GCA_037043045.1_ASM3704304v1_genomic.fna"

# Send myself a notification that the script is starting
curl -d "🔔 Starting RepeatModeler for $database_name at $(date). Reference genome: $reference_genome. Check logs at $output_directory/Logs/ for details." \
	 ntfy.sh/kaas-ballard-Robin-scripts-72724027625978

# Change the directory to the output directory so that the RM_ files are created there
cd "$output_directory" || { echo "Failed to change directory to $output_directory"; exit 1; }

# Step #1: Build a new RepeatModeler database for the reference genome
# Note: BuildDatabase can build directories that don't exist yet
if [ ! -f "$output_directory/$database_name.nsq" ]; then
	echo "Warning: Database at the specified path does not exist. Building a new database now."
	BuildDatabase -name "$output_directory/$database_name" \
		"$reference_genome" 2>&1 | tee "$output_directory/Logs/BuildDatabase.log"
else
	echo "Database already exists at the specified path. Skipping."
fi

# Step #2: Detect if the database was created and then run RepeatModeler
if [ ! -f "$output_directory/$database_name.nsq" ]; then
	echo "Error: Database was not created. Exiting."
	exit 1
else
	# Run RepeatModeler
	RepeatModeler \
		-threads "$threads" \
		-database "$output_directory/$database_name" 2>&1 | tee "$output_directory/Logs/RepeatModeler.log"
fi

# Send a notification that the script has finished
if [ $? -eq 0 ]; then # if [ $? -eq 0 ] checks if the last command was successful
    curl -d "✅ SUCCESS: RepeatModeler completed for $database_name at $(date). Check logs at $output_directory/Logs/" \
         ntfy.sh/kaas-ballard-Robin-scripts-72724027625978
else
    curl -d "❌ FAILED: RepeatModeler failed for $database_name at $(date). Check logs for errors." \
         ntfy.sh/kaas-ballard-Robin-scripts-72724027625978
fi