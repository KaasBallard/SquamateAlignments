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
	seqkit tab2fx > "$repeat_modeler_families_prefix"

# Set the file name for unknown elements
unknown_elements_file="$repeat_modeler_families_prefix.unknown"

# Set the file name for known elements
known_elements_file="$repeat_modeler_families_prefix.known"

# Step #2: Split elements into known and unknown classification files
# Find the known elements
cat "$repeat_modeler_families_prefix" | seqkit fx2tab | \
	grep -v "Unknown" | seqkit tab2fx > "$known_elements_file"
# Find the unkown elements
cat "$repeat_modeler_families_prefix" | seqkit fx2tab | \
	grep "Unknown" | seqkit tab2fx > "$unknown_elements_file"

# Step #1: Run repclassifier with the Repbase database
repclassifier -t "$threads" \
	-d Tetrapoda \
	-u ../00_InitialRepeatModellerRun/croAtr2-families.prefix.fa.unknown \
	-k ../00_InitialRepeatModellerRun/croAtr2-families.prefix.fa.known \
	-a ../00_InitialRepeatModellerRun/croAtr2-families.prefix.fa.known \
	-o round-01_RepbaseTetrapoda-Self

../../../InstallationFiles/repclassifier -t 3 -u round-01_RepbaseTetrapoda-Self/round-01_RepbaseTetrapoda-Self.unknown -k round-01_RepbaseTetrapoda-Self/round-01_RepbaseTetrapoda-Self.known -a round-01_RepbaseTetrapoda-Self/round-01_RepbaseTetrapoda-Self.known -o round-02_Self

../../../InstallationFiles/repclassifier -t 3 -u round-02_Self/round-02_Self.unknown -k round-02_Self/round-02_Self.known -a round-02_Self/round-02_Self.known -o round-03_Self

../../../InstallationFiles/repclassifier -t 3 -u round-03_Self/round-03_Self.unknown -k round-03_Self/round-03_Self.known -a round-03_Self/round-03_Self.known -o round-04_Self

../../../InstallationFiles/repclassifier -t 3 -u round-04_Self/round-04_Self.unknown -k round-04_Self/round-04_Self.known -a round-04_Self/round-04_Self.known -o round-05_Self

../../../InstallationFiles/repclassifier -t 3 -u round-05_Self/round-05_Self.unknown -k 18Snakes.Known.clust.fasta -a round-05_Self/round-05_Self.known -o round-06_18Snakes

../../../InstallationFiles/repclassifier -t 3 -u round-06_18Snakes/round-06_18Snakes.unknown -k round-06_18Snakes/round-06_18Snakes.known -a round-06_18Snakes/round-06_18Snakes.known -o round-07_Self

../../../InstallationFiles/repclassifier -t 3 -u round-07_Self/round-07_Self.unknown -k round-07_Self/round-07_Self.known -a round-07_Self/round-07_Self.known -o round-08_Self

../../../InstallationFiles/repclassifier -t 3 -u round-08_Self/round-08_Self.unknown -k round-08_Self/round-08_Self.known -a round-08_Self/round-08_Self.known -o round-09_Self

../../../InstallationFiles/repclassifier -t 3 -u round-09_Self/round-09_Self.unknown -k round-09_Self/round-09_Self.known -a round-09_Self/round-09_Self.known -o round-10_Self
