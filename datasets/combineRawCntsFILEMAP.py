#!/usr/bin/python

import sys
import glob
import pandas as pd

from collections import defaultdict

# Funstions to run the script==================================
def get_patient_id():

    patientID = defaultdict(str)

    #""" Read the FILE_SAMPLE_MAP.txt FILE to store filename and id value"""
    with open('FILE_SAMPLE_MAP.txt') as f:

        lines = f.readlines()

        for line in lines[1:]:
#            print(len(line))
            file_id = line.split("\t")[0]
            sample_id = line.split("\t")[1]

            #print(file_id, 'and', sample_id.rstrip())
            patientID[file_id] = sample_id.rstrip()
        # To view the first line
#    print(patientID.items()[0])
    return patientID

def get_gene_cnts_dict(file_name):
    # Read the file and store the contents in a list
    with open(file_name) as file_object:

        # Create a geneCounts dictionary to store gene:raw_counts
        geneCounts = defaultdict(int)

        # Read and store each line in a list
        lines = file_object.readlines()

        # for each line store the gene and raw_counts in geneCounts dictionary
        # skip the first line i.e. header
        for line in lines[1:]:
            gene = line.split("\t")[0]
            raw_counts = line.split("\t")[1]

            geneCounts[gene] = raw_counts

    # To view the first line
#   print(geneCounts.items()[1])
    return geneCounts

#================================================================

# Variables
patientCounts = {}
pat_ID = {}
pat_ID = get_patient_id()

for filename in glob.glob('*.rsem.genes.results'):
    print 'Current file being read:', str(filename)

    # Get the patient id from the filename
    id = pat_ID.get(filename, "not found")
    print(id)

    # Get the gene counts dictionary for the patient
    geneCountsid = get_gene_cnts_dict(filename)

    # Print the first line of the geneCounts dict returned
    print(geneCountsid.items()[1])

    # Store the each geneCounts dict per patient
    patientCounts[id] = geneCountsid
    #print(patientCounts.items()[0])

# Create a table of the results with the missing values replaced with 0
df = pd.DataFrame(patientCounts)
df.fillna(0, inplace=True)
print(df)


# Write the data frame to a file
outfile = sys.argv[1]

df.to_csv(outfile)

