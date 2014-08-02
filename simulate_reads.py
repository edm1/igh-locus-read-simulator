#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
#

import argparse
from src.bio_file_parsers import fasta_parser
from src.bio_file_parsers import write_fasta
from src.bio_file_parsers import reverse_complement
import re
import sys
import numpy
import subprocess

def main():

    # Parse args
    global args
    args = parse_arguments()

    # Check proportions sum to less than 100
    if sum(args.Proportions) >= 100:
        sys.exit('Proportions must sum to <100.')

    # Load germline sequences
    germ_seqs = {}
    germ_seqs['V'] = load_germline_seqs(args.Vgenes)
    germ_seqs['D'] = load_germline_seqs(args.Dgenes)
    germ_seqs['J'] = load_germline_seqs(args.Jgenes, is_j=True)

    # Create the fasta file of simulated reads
    fasta_file = args.OutPrefix + '.fasta'
    with open(fasta_file, 'w') as out_h:
        cluster_num = 1
        remaining_reads = args.NumReads

        # Do clonal reads
        for prop in args.Proportions:
            # Calc number of reads to make
            clus_num_reads = int((float(prop)/100) * args.NumReads)
            remaining_reads = remaining_reads - clus_num_reads
            # Make the clone
            clone = Clone(germ_seqs, cluster_num, rev_comp=True)
            # Write it to the fasta file multiple times
            for i in range(clus_num_reads):
                title = 'Cluster_{0}_{1}'.format(cluster_num, i)
                write_fasta(out_h, title, clone.seq)
            cluster_num += 1

        # Fill remaining with unique reads
        for i in range(remaining_reads):
            clone = Clone(germ_seqs, cluster_num, rev_comp=True)
            title = 'Cluster_{0}_{1}'.format(cluster_num, 0)
            write_fasta(out_h, title, clone.seq)
            cluster_num += 1

    # Run ART to simulate MiSeq reads
    fastq_file = args.OutPrefix + '.fastq'
    cmd = [args.ARTbin,
           '-amp',
           '-na',
        #    '-sam',
           '-f 1',
           '-l {0}'.format(args.TotalReadLen),
           '-i {0}'.format(fasta_file),
           '-o {0}'.format(args.OutPrefix)
           ]
    cmd = ' '.join(cmd)
    subprocess.call(cmd, shell=True)

    print cmd


# A clone is a somatic recombination of the VDJ sequences
class Clone:
    def __init__(self, germ_seqs, title, rev_comp=False):
        """ Generates random recombination and sequence
        """
        self.title = title
        self.rev_comp = rev_comp
        # Select random V, D and J
        self.J = numpy.random.choice(germ_seqs['J'].keys())
        self.D = numpy.random.choice(germ_seqs['D'].keys())
        self.V = numpy.random.choice(germ_seqs['V'].keys())
        # Generate del sizes from poisson distribution
        self.J5del = numpy.random.poisson(args.J5del)
        self.D3del = numpy.random.poisson(args.D3del)
        self.D5del = numpy.random.poisson(args.D5del)
        self.V3del = numpy.random.poisson(args.V3del)
        # Generate insertion sizes from poisson distribution
        self.VDins = numpy.random.poisson(args.VDins)
        self.DJins = numpy.random.poisson(args.DJins)

        # Create the sequence
        self.make_seq(germ_seqs, rev_comp)
        self.write_log()

    def make_seq(self, germ_seqs, rev_comp):

        # Start with J sequence with 5' deleted
        seq = germ_seqs['J'][self.J][self.J5del:]
        # Add DJ random insert to 5'
        ins = gen_random_insertion(self.DJins)
        seq =  ins + seq
        # Add D seq with deletions at both ends
        d_seq = germ_seqs['D'][self.D]
        seq = d_seq[self.D5del:len(d_seq)-self.D3del] + seq
        # Add VD random insert to 5'
        ins = gen_random_insertion(self.VDins)
        seq = ins + seq
        # Add V with 3' del
        v_seq = germ_seqs['V'][self.V]
        v_seq = v_seq[:len(v_seq)-self.V3del]
        seq = v_seq + seq
        # ART requires the seq to be at least the max length
        if len(seq) < args.TotalReadLen:
            ins = gen_random_insertion(args.TotalReadLen - len(seq) + 1)
            seq = ins + seq
        # Clip the sequence to total read length
        if not len(seq) == args.TotalReadLen:
            seq = seq[len(seq)-args.TotalReadLen-1:]
        # Save to self
        if rev_comp == True:
            self.seq = reverse_complement(seq)
        else:
            self.seq = seq

    def write_log(self):
        log_file = args.OutPrefix + '.log'
        with open(log_file, 'a') as out_h:
            out = ['> Cluster {0}'.format(self.title),
                   self.seq,
                   'V gene: {0}'.format(self.V),
                   'D gene: {0}'.format(self.D),
                   'J gene: {0}'.format(self.J),
                   "V3' del: {0}".format(self.V3del),
                   "VD ins: {0}".format(self.VDins),
                   "D5' del: {0}".format(self.D5del),
                   "D3' del: {0}".format(self.D3del),
                   "DJ ins: {0}".format(self.DJins),
                   "J5' del: {0}".format(self.J5del),
                   "Reverse complement: {0}".format(self.rev_comp)
                  ]
            out_h.write('\n'.join(out) + '\n\n')

def gen_random_insertion(length):
    bases = ['A', 'T', 'G', 'C']
    return ''.join(numpy.random.choice(bases, length))


def load_germline_seqs(filen, is_j=False):
    """ Loads the fasta file containing the VDJ germline sequences. If it is J,
        the 3' after CTGGGG is removed.
    """
    seqs = {}

    # Parse fasta
    with open(filen, 'r') as in_h:
        for title, seq in fasta_parser(in_h):

            title = title.split(' ')[0]

            # Clip off 3' if J seq
            if is_j == True:
                m = re.search(r'([A|T|G|C]+CTGGGG[A|T|G|C]{5}).*', seq)
                if m:
                    seq = m.group(1)
                else:
                    sys.exit('CTGGGG could not be found in {0}'.format(title))

            seqs[title] = seq

    return seqs

def parse_arguments():
    """ Load command line args.
    """

    parser = argparse.ArgumentParser(description='Simulate IGH MiSeq reads.')

    # Required args
    parser.add_argument('Proportions',
                        metavar='N',
                        type=int,
                        nargs='+',
                        help=('Percentage of reads that will make up each' +
                              'clonal group. Must sum to <100. e.g. 20 10 5'))
    parser.add_argument('--OutPrefix',
                        metavar='<str>',
                        type=str,
                        required=True,
                        help='Output prefix for simulated fasta, fastq and log')

    # Optional arguments
    parser.add_argument('--NumReads',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=100000,
                        help="Number of reads to simulate. (100,000)")
    parser.add_argument('--Vgenes',
                        metavar='<fasta>',
                        type=str,
                        required=False,
                        default='VDJ_germline_sequences/V_genes.fasta',
                        help='Location of fasta containing V germline sequences.')
    parser.add_argument('--Dgenes',
                        metavar='<fasta>',
                        type=str,
                        required=False,
                        default='VDJ_germline_sequences/D_genes.fasta',
                        help='Location of fasta containing D germline sequences.')
    parser.add_argument('--Jgenes',
                        metavar='<fasta>',
                        type=str,
                        required=False,
                        default='VDJ_germline_sequences/J_genes.fasta',
                        help='Location of fasta containing J germline sequences.')
    parser.add_argument('--J5del',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=2,
                        help="Mean J 5' deletion size. (2)")
    parser.add_argument('--D3del',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=2,
                        help="Mean D 3' deletion size. (2)")
    parser.add_argument('--D5del',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=2,
                        help="Mean D 5' deletion size. (2)")
    parser.add_argument('--V3del',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=2,
                        help="Mean V 3' deletion size. (2)")
    parser.add_argument('--VDins',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=4,
                        help="Mean V -> D insertion size. (4)")
    parser.add_argument('--DJins',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=4,
                        help="Mean D -> J insertion size. (4)")
    parser.add_argument('--TotalReadLen',
                        metavar='<int>',
                        type=int,
                        required=False,
                        default=250,
                        help="Total read length, surplus bases will be clipped from V 5'. Max 250. (250)")
    parser.add_argument('--ARTbin',
                        metavar='<str>',
                        type=str,
                        required=False,
                        default='/home/ed/Programs/art_bin_VanillaIceCream/art_illumina',
                        help="Location of ART Illumina binary file")
    # parser.add_argument('--numCPU',
    #                     metavar='<int>',
    #                     type=int,
    #                     required=False,
    #                     default=1,
    #                     help='Number of CPUs to use. (1)')


    return parser.parse_args()

if __name__ == '__main__':
    main()
