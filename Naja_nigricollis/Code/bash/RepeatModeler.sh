#!/bin/bash

: <<'ScriptDescription'
Date: 2025/03/11

This script is designed to run RepeatModeler. I am following the method found in Daren Card's blog post:
https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

ScriptDescription

# Activate the mamba environment containing RepeatModeler
source /home/administrator/mambaforge/bin/activate RepeatMaskAnnot

# Set the output directory
output_directory="/home/administrator/ExtraSSD2/Kaas/Projects/SquamateAlignments/Naja_nigricollis/Results/SoftMasking/0_RepeatModeller"

# Create log directory under the output directory if it does not exist
[ ! -d "$output_directory/Logs" ] && mkdir -p "$output_directory/Logs"

# Build a new RepeatModeler database for the reference genome
BuildDatabase -name \
	-engine 