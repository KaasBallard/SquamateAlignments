#!/bin/bash

: <<'ScriptDescription'
Date: 2025/03/11

This script is designed to run RepeatModeler. I am following the method found in Daren Card's blog post:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

ScriptDescription

# Activate the mamba environment containing RepeatModeler
source /home/administrator/mambaforge/bin/activate RepeatMaskAnnot

# Set the number of threads for RepeatModeler
threads=40

# Set the output directory
output_directory="/home/administrator/ExtraSSD2/Kaas/Projects/SquamateAlignments/SoftMasking/Naja_nigricollis/Results/0_RepeatModeller"

# Create log directory under the output directory if it does not exist
[ ! -d "$output_directory/Logs" ] && mkdir -p "$output_directory/Logs"

# Set the name for the RepeatModeler database
database_name="Naja_nigricollis"

# Reference genome
reference_genome="/home/administrator/ExtraSSD2/Kaas/Projects/SquamateAlignments/Reference_Genomes/from_Sekar/_completeAssemblies/Naja_nigricollis_najNig1/Assembly/najNig2.ragtag.scaffold_naNa.REHEADER.MT.fasta"

# Step #1: Build a new RepeatModeler database for the reference genome
# Note: BuildDatabase can build directories that don't exist yet
if [ ! -f "$output_directory/$database_name.nsq" ]; then
	echo "Warning: Database at the specified path does not exist. Building a new database now."
	BuildDatabase -name "$output_directory/$database_name" \
		-engine rmblast \
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
	RepeatModeler -pa "$threads" \
		-engine rmblast \
		-database "$output_directory/$database_name" 2>&1 | tee "$output_directory/Logs/RepeatModeler.log"
fi