#!/bin/bash

: <<'ScriptDescription'
Date: 2025/05/16

This script is designed to run RepeatModeler. I am following the method found in Daren Card's blog post:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

ScriptDescription

# Activate the mamba environment containing RepeatModeler
source "$HOME/miniforge3/bin/activate" RepeatMaskAnnot

# Set the number of threads for RepeatModeler
threads=30

# Set the output directory
output_directory="$HOME/Documents/Kaas/SquamateAlignments/SoftMasking/Viperidae/Vipera/Vipera_latastei/Results/0_RepeatModeler"

# Create log directory under the output directory if it does not exist
[ ! -d "$output_directory/Logs" ] && mkdir -p "$output_directory/Logs"

# Set the name for the RepeatModeler database
database_name="Vipera_latastei"

# Reference genome
reference_genome="$HOME/Documents/Kaas/SquamateAlignments/Reference_Genomes/Public_Genomes/Viperidae/Vipera/Vipera_latastei/ncbi_dataset/data/GCA_024294585.1/GCA_024294585.1_rVipLat1.pri_genomic.fna"

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