#!/home/administrator/mambaforge/envs/biopython/bin/python

# region Import Modules
# Import required modules
import os
import sys
import re
import argparse
from Bio import SeqIO
import numpy as np
import pandas as pd
import polars as pl
import biopython as bp
import pyfasta as pf
# endregion


# region Function Definitions

# Function #1: fasta_to_dataframe
def fasta_to_dataframe(fasta):
	"""
	Converts a FASTA file to a polars DataFrame object. 

	Args:
		fasta (str): Path to the FASTA file.

	Returns:
		pl.DataFrame: A DataFrame with 'header' and 'sequence' columns.
	"""

	# Parse the FASTA and extract headers and sequences
	records = SeqIO.parse(fasta, 'fasta')
	# Initialize empty lists to store headers and sequences
	headers = []
	sequences = []

	# Iterate over the records and append the headers and sequnences to the respective lists
	for record in records:
		headers.append(str(record.id))
		sequences.append(str(record.seq))
	
	# Create a polars data frame from the headers and sequences
	df = pl.DataFrame({'headers': headers, 'sequences': sequences})

	return df