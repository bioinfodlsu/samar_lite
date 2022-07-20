#%%
import itertools
import csv
import argparse
import sys
import ast
from Bio import SeqIO

class Tab_alns:
    def __init__(self,tab_line):
        self.ref = tab_line[1]
        self.ref_aln_width = int(tab_line[3])

    def __str__(self):
        return 'ref:{self.ref} ref_aln_width:{self.ref_aln_width}'.format(self=self)

class Counts:
    def  __init__(self,ref_length): 
        self.length = ref_length
        self.unique_count = 0.0
        self.unique_count_norm = 0.0
        self.final_count = 0.0
        self.final_count_norm = 0.0
        self.tpm = 0.0
    
    def __str__(self):
        return 'length:{self.length}    unique_count:{self.unique_count}    unique_count_norm:{self.unique_count_norm}   final_count:{self.final_count}    final_count_norm:{self.final_count_norm}    tpm:{self.tpm}'.format(self=self)

def update_unique_single_end(aln,counts_dict):
    counts_dict[aln].unique_count += 1 #update count

def unique_pass_SE(input_alns, counts_dict):

    with open(input_alns) as infile:
        for line in infile:
            
            alns_str = line.split("\t")[1]
            alns = ast.literal_eval(alns_str)  # alignments of read 1

            alns = list(alns)
            

            # if it's a unique mapping
            if len(alns) == 1 :
                update_unique_single_end(alns[0],counts_dict)


    # done with first pass. 1. set final_count to be unique_count and 2.normalize unique counts by length (of non-zero counts)
    for k, v in counts_dict.items():

        v.final_count = v.unique_count

        if v.unique_count != 0:
            v.unique_count_norm = v.unique_count / v.length
        else:
            v.unique_count_norm = 0

def update_rescue_single_end(alns,counts_dict):
    unique_norm_counts = [counts_dict[ref_id].unique_count_norm for ref_id in alns]
    denom = sum(unique_norm_counts)
    if denom != 0 :
        props = [c/denom for c in unique_norm_counts]
        for i,ref_id in enumerate(alns):
            counts_dict[ref_id].final_count += props[i]
    else:
        denom = len(alns)
        for ref_id in alns:
            counts_dict[ref_id].final_count += 1/denom

def rescue_pass_SE(input_alns,counts_dict):

    with open(input_alns) as infile:
        for line in infile:
            alns_str = line.split("\t")[1]

            alns = ast.literal_eval(alns_str)
            alns = list(alns)

            if len(alns) > 1 :
                update_rescue_single_end(alns,counts_dict)
            else:
                continue

def main(reference,input_alns, out_counts, single_end, frag_len_mean=0, frag_len_std=0):
    '''
    Does 2 passes over the alignments.
    In the first pass, we only consider reads with unique alignments. Counts are recorded in the dict unique.
    In the second pass, we consider the remaining reads. Counts are distributed based on the proportion of uniquely mapped reads normalized by length.
    '''

    # initiate counts_dict
    counts_dict = {}  # key = peptide ID, value = object of class Counts
    for seq_record in SeqIO.parse(reference, "fasta"):
        counts_dict[seq_record.id] = Counts(len(seq_record))

    if single_end == "True" : #Looking at you, argparse.
        unique_pass_SE(input_alns,counts_dict)
        rescue_pass_SE(input_alns,counts_dict)

    else:
        sys.exit("Unexpected string for endendness")

    #compute TPM
    #read/length
    scaling_factor = 0.0
    for counts in counts_dict.values():
        counts.final_count_norm = counts.final_count/counts.length
        scaling_factor += counts.final_count_norm

    for counts in counts_dict.values():
        if scaling_factor != 0:
            counts.tpm = counts.final_count_norm * 1000000/scaling_factor


    with open(out_counts,"w") as tab_file:
        writer = csv.writer(tab_file, delimiter='\t')

        #adding Header
        writer.writerow(["Name","Length","EffectiveLength","TPM", "NumReads"])

        for key, counts in sorted(counts_dict.items()):
            writer.writerow([key,counts.length,counts.length,counts.tpm,counts.final_count])

    return counts_dict

#%%
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('reference',help="fasta file containing the reference protein sequences")
    parser.add_argument('input_file', help="alignments in tab format")
    parser.add_argument('output_file',help="output file containing counts")
    parser.add_argument('single_end',help="true if reads are single-end and not paired-end")
    parser.add_argument('--frag_len_mean',type=float,help="mean of fragment length, when translated. Use the same value as last-pair-probs")
    parser.add_argument('--frag_len_std',type=float,help="standard deviation of fragment length, when translated. Use the same value as last-pair-probs") 

    args = parser.parse_args()
    main(args.reference,args.input_file,args.output_file, args.single_end, args.frag_len_mean, args.frag_len_std )