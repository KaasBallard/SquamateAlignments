#!/bin/bash

: <<'ScriptDescription'
Date: 2025/03/11
This script is designed to repclassifier, which is a piece of softwart that uses Repbase database and RepeatMasker to try and identify unknown elements with sequence similarity to curated repeat elements in Repbase.
It has another mode that can uses a custom library, in the form of a fasta file. This script will use both modes.
The program can be found here:
https://github.com/darencard/GenomeAnnotation/blob/master/repclassifier

The blog that explains the program can be found here:

ScriptDescription

../../../InstallationFiles/repclassifier -t 3 -d Tetrapoda -u ../00_InitialRepeatModellerRun/croAtr2-families.prefix.fa.unknown -k ../00_InitialRepeatModellerRun/croAtr2-families.prefix.fa.known -a ../00_InitialRepeatModellerRun/croAtr2-families.prefix.fa.known -o round-01_RepbaseTetrapoda-Self

../../../InstallationFiles/repclassifier -t 3 -u round-01_RepbaseTetrapoda-Self/round-01_RepbaseTetrapoda-Self.unknown -k round-01_RepbaseTetrapoda-Self/round-01_RepbaseTetrapoda-Self.known -a round-01_RepbaseTetrapoda-Self/round-01_RepbaseTetrapoda-Self.known -o round-02_Self

../../../InstallationFiles/repclassifier -t 3 -u round-02_Self/round-02_Self.unknown -k round-02_Self/round-02_Self.known -a round-02_Self/round-02_Self.known -o round-03_Self

../../../InstallationFiles/repclassifier -t 3 -u round-03_Self/round-03_Self.unknown -k round-03_Self/round-03_Self.known -a round-03_Self/round-03_Self.known -o round-04_Self

../../../InstallationFiles/repclassifier -t 3 -u round-04_Self/round-04_Self.unknown -k round-04_Self/round-04_Self.known -a round-04_Self/round-04_Self.known -o round-05_Self

../../../InstallationFiles/repclassifier -t 3 -u round-05_Self/round-05_Self.unknown -k 18Snakes.Known.clust.fasta -a round-05_Self/round-05_Self.known -o round-06_18Snakes

../../../InstallationFiles/repclassifier -t 3 -u round-06_18Snakes/round-06_18Snakes.unknown -k round-06_18Snakes/round-06_18Snakes.known -a round-06_18Snakes/round-06_18Snakes.known -o round-07_Self

../../../InstallationFiles/repclassifier -t 3 -u round-07_Self/round-07_Self.unknown -k round-07_Self/round-07_Self.known -a round-07_Self/round-07_Self.known -o round-08_Self

../../../InstallationFiles/repclassifier -t 3 -u round-08_Self/round-08_Self.unknown -k round-08_Self/round-08_Self.known -a round-08_Self/round-08_Self.known -o round-09_Self

../../../InstallationFiles/repclassifier -t 3 -u round-09_Self/round-09_Self.unknown -k round-09_Self/round-09_Self.known -a round-09_Self/round-09_Self.known -o round-10_Self
