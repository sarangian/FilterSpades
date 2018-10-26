#!/usr/bin/env python

import re
import argparse, sys
import pandas as pd

def estimated_genome_range(log_file="spades.log"):

    log = open(log_file, "r")
    lower_bound = 0
    upper_bound = 0
    genome_length = 0
    kmer = 0
    genome_sizes = []
    for line in log:
        if re.search("Estimated genome size", line):
            words = line.rstrip().split()
            last_word = words[-1]
            size = int(last_word)
            genome_sizes.append(size)

        if re.search("k:", line):
            kmers = [x.rstrip() for x in line.split(':')]
            kmer=kmers[-1]
            print ("List of K-Mers Used:   ", kmer)

    sorted_g_sizes = sorted(genome_sizes)
    mindex = (len(sorted_g_sizes) / 2)
    #lastindex = len(sorted_g_sizes) - 2
    lower_bound_index=mindex -3
    upper_bound_index=mindex +3
    genome_length=sorted_g_sizes[mindex]
    lower_bound=sorted_g_sizes[lower_bound_index]
    upper_bound=sorted_g_sizes[upper_bound_index]

    print ("Estimated Genome Sizes: ", genome_sizes)

    print ("   Sorted Genome Sizes: ", sorted_g_sizes)
    print ("\n")
    print ("Median genome size index: ", mindex)
    print ("\n")
    print ("             lower bound: ", sorted_g_sizes[lower_bound_index])
    print ("    Median Genome Length: ", genome_length)
    print ("             Upper Bound: ", sorted_g_sizes[upper_bound_index])

    print ("\n")
    return lower_bound, upper_bound, genome_length



def process(args):

    keep=True
    l1=0
    lengths=[]
    if args.input == sys.stdin:
        print("Reading from stdin...")
        #print("Reading from stdin...", file = sys.stderr)
    for l in args.input:
        if l[0] == '>':
            parts = l.split('_')
            if parts[2]!='length' or parts[4]!='cov':
                raise RuntimeError("Invald syntax %s"%parts)
            keep = int(parts[3])>=args.length and float(parts[5])>=args.cov
            if keep:
                lengths.append(int(parts[3]))
        if keep and args.output:
            args.output.write(l)

    #Calculate N50
    tot = sum(lengths)
    s=0
    n50 =0
    for l in sorted(lengths):
        s+=l
        if s>tot/2:
            n50 = l
            break

    info_out = sys.stderr if args.output else sys.stdout
    data_string = "K-Mer Coverage: {}\t Min Contig Length: {}K\t Number of Contigs: {} \t Genome Size: {:.2f} mb\t N50: {}".format(args.cov, args.length/1000.0, len(lengths), tot/1000000.0, n50)

    contig_list=sorted(lengths)
    
    print(data_string)
    return tot, len(lengths)


def arguments():
    parser = argparse.ArgumentParser(description='Filter spades output file (scaffolds.fasta)')
    parser.add_argument('--spadeslog', default="spades.log",
                        help="spades.log file to fetch estimated genome size")
    parser.add_argument('--sample_out', default="sample_",
                        help="Output file name of filtered scaffolds")

    parser.add_argument('input', type=argparse.FileType('r'),
                        nargs='?', default='-',
                        help="Input fasta file (default stdin)")
    parser.add_argument('--cov', default=0, type=float,
                        help="Filter contigs with coverage less than this")
    parser.add_argument('--length', default=0, type=int,
                        help="Filter contigs with length less than this")
    parser.add_argument('--output', default=None, type=argparse.FileType('w'),
                        help="Output file.  Default is to just print summary")



    return parser

if __name__ == '__main__':
    parser = arguments()
    args = parser.parse_args()

    lower_bound, upper_bound, genome_length = estimated_genome_range(args.spadeslog)
    #print(lower_bound, genome_length, upper_bound)
    print ("Filtrartion of spades scaffold started based on combination of k-mer coverage starting from 1 to 20 (with step size 0.1)\nand minimum contig length 1000 to 5000 (with step size 500) \n") 

    contig_df = pd.DataFrame(columns=['cov', 'length', 'contig_count'])
    for c in range(10, 201, 1):
        cov = c / 10.0
        for length in range(1000, 5001, 500):
            args.cov = cov
            args.length = length
            genome_length,contig_count = process(args)
            args.input.seek(0)
            if lower_bound <= genome_length <= upper_bound:
                contig_df = contig_df.append({'cov': cov, 'length': length, 'contig_count': contig_count}, ignore_index=True)

    median_contig_ln = (contig_df['contig_count'].median())
    flt_med_df = contig_df[contig_df['contig_count'] <= median_contig_ln]

    temp_df_1 = flt_med_df[flt_med_df['cov'].isin(flt_med_df.groupby('contig_count').max()['cov'].values)]
    temp_df_2 = temp_df_1[temp_df_1['length'].isin(temp_df_1.groupby('contig_count').max()['length'].values)]
    final_df  = temp_df_2[temp_df_2['contig_count'].isin(temp_df_1.groupby('length').min()['contig_count'].values)]
    print ("\n")
    print ("Median Contig Length:", median_contig_ln)
    #print (final_df)
    print ("\n")
    print ("Filtered Contigs")
    filtered_sacffold_count=1
    for index, row in final_df.iterrows():
        cov = row['cov']
        length = row['length']
        contig_count = row['contig_count']
        filename = "{}_{}.fna".format(args.sample_out, filtered_sacffold_count)
        # filename = "COV_{}_MLN_{:.0f}_NOC_{:.0f}.fna".format(cov, length, contig_count)
        out_file = open(filename, "w")
        args.length = length
        args.cov = cov
        args.output = out_file
        process(args)
        args.output = None
        args.input.seek(0)
        filtered_sacffold_count=filtered_sacffold_count+1





