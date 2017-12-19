#! /usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function
from itertools import cycle
import sys
from anarci import anarci
import multiprocessing
import scipy


gencode = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
           'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
           'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
           'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',

           'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
           'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
           'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
           'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',

           'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
           'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
           'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
           'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',

           'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
           'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
           'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
           'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

stop_codons = {'TAA': 1, 'TGA': 1, 'TAG': 1}

ACGT = ['A', 'C', 'G', 'T']
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
AA_INDEX = {aa: i for i, aa in enumerate(AA_LIST)}

AHO_L = 149


def run_anarci(sequences):
    '''Run ANARCI annotation to get Aho numbering.'''
    allowed_species = 'human'
    allow = 'H'
    ncpu = 1
    scheme = 'aho'
    numbered, alignment_details, hit_tables = anarci(sequences, scheme=scheme, output=False, allow=allow, ncpu=ncpu, allowed_species=allowed_species)
    return numbered


def convert2profile(numbers):
    '''Convert the ANARCI output to a profile.'''
    profile = [[0]*21 for fi in range(AHO_L)]
    for ai in numbers[0][0]:  # The last 0 is skipping the start, end from "validate_numbering"
        i = int(ai[0][0]) - 1
        assert(int(ai[0][0]) <= AHO_L)  # max AHO_L in aho numbering
        try:
            aa_idx = AA_INDEX[ai[1]]
        except:
            print(numbers)
            print(ai)
            print(info_i_j[ii])
            sys.exit()
        profile[i][aa_idx] += 1
    return profile


class MutationModel():
    '''
    A class for a mutation model, and functions to mutate sequences.
    Based on GCtree code, modified to run 20x faster and tested
    to reproduce results from the original GCtree code.
    '''
    def __init__(self, sequence, mutability_file=None, substitution_file=None, mutation_order=False, with_replacement=False):
        self.sequence_length = len(sequence)
        self.mutability_p_master = scipy.zeros((self.sequence_length))
        self.substitution_p_master = scipy.zeros((self.sequence_length, 4))

        self.mutation_order = mutation_order
        self.with_replacement = with_replacement
        if mutability_file is not None and substitution_file is not None:
            self.context_model = {}
            with open(mutability_file, 'r') as f:
                f.readline()  # Eat header
                for line in f:
                    motif, score = line.replace('"', '').split()[:2]
                    self.context_model[motif] = float(score)

            # kmer k
            self.k = None
            with open(substitution_file, 'r') as f:
                # eat header
                f.readline()
                for line in f:
                    fields = line.replace('"', '').split()
                    motif = fields[0]
                    if self.k is None:
                        self.k = len(motif)
                        assert self.k % 2 == 1
                    else:
                        assert len(motif) == self.k
                    self.context_model[motif] = (self.context_model[motif], [float(x) for x in fields[1:5]])
        else:
            self.context_model = None


    def mutabilities(self, sequence, update_pos=None, mutability_p=None, substitution_p=None):
        '''returns the mutability of a sequence at each site, along with nucleotide biases'''
        if update_pos is None:
            mutability_p = self.mutability_p_master.copy()
            substitution_p = self.substitution_p_master.copy()
            update_pos = {i:1 for i in range(self.sequence_length)}

        # ambiguous left end motifs
        for i in range(self.k//2):
            if i not in update_pos:
                continue
            kmer_suffix = sequence[:(i+self.k//2+1)]
            matches = [value for key, value in self.context_model.iteritems() if key.endswith(kmer_suffix)]
            len_matches = len(matches)
            # use mean over matches
            mutability = sum(match[0] for match in matches)/len_matches
            substitution = [sum(d[1][n] for d in matches)/len_matches for n in range(4)]
            mutability_p[i] = mutability
            substitution_p[i] = substitution[:]
        # unambiguous internal kmers
        for i in range(self.k//2, self.sequence_length - self.k//2):
            if i not in update_pos:
                continue
            mutability_p[i] = self.context_model[sequence[(i-self.k//2):(i+self.k//2+1)]][0]
            substitution_p[i] = self.context_model[sequence[(i-self.k//2):(i+self.k//2+1)]][1]
        # ambiguous right end motifs
        for i in range(self.sequence_length - self.k//2, self.sequence_length):
            if i not in update_pos:
                continue
            kmer_prefix = sequence[(i-self.k//2):]
            matches = [value for key, value in self.context_model.iteritems() if key.startswith(kmer_prefix)]
            len_matches = len(matches)
            # use mean over matches
            mutability = sum(match[0] for match in matches)/len_matches
            substitution = [sum(d[1][n] for d in matches)/len_matches for n in range(4)]
            mutability_p[i] = mutability
            substitution_p[i] = substitution[:]

        return mutability_p, substitution_p


    def mutate(self, sequence, mutability_p, substitution_p, m=1):
        """
        Mutate a sequence, with lambda0 the baseline mutability
        Cannot mutate the same position multiple times
        @param sequence: the original sequence to mutate
        @param m: number of mutations to perform
        @param frame: the reading frame index
        """
        trials = 50
        if '*' in [gencode[sequence[si:(si+3)]] for si in range(0, len(sequence), 3)]:
            raise RuntimeError('sequence contains stop codon!')

        for i in range(m):
            sequence_list = list(sequence)
            for trial in range(1, trials+1):
                # Determine the position to mutate from the mutability matrix
                mut_pos = scipy.random.multinomial(1, mutability_p/mutability_p.sum()).argmax()
                # Now draw the target nucleotide using the substitution matrix
                chosen_target = scipy.random.multinomial(1, substitution_p[mut_pos]).argmax()
                original_base = sequence_list[mut_pos]
                sequence_list[mut_pos] = ACGT[chosen_target]
                mut_pos_frame = mut_pos % 3
                if ''.join(sequence_list[(mut_pos-mut_pos_frame):(mut_pos-mut_pos_frame+3)]) not in stop_codons:
                    sequence = ''.join(sequence_list) # reconstruct our sequence
                    update_pos = {up: 1 for up in range((mut_pos-self.k//2), (mut_pos+self.k//2+1))}
                    mutability_p, substitution_p = self.mutabilities(sequence, update_pos=update_pos, mutability_p=mutability_p, substitution_p=substitution_p)
                    break
                if trial == trials:
                    raise RuntimeError('stop codon in simulated sequence on '+str(trials)+' consecutive attempts')
                sequence_list[mut_pos] = original_base # <-- we only get here if we are retrying
            sequence = ''.join(sequence_list) # reconstruct our sequence
        return sequence


    def simulate_AAprofile(self, sequence, numb_profile, muts_iter, N=None, S=None, verbose=False):
        '''
        Simulate neutral amino acid substitution profile under a k-mer motif based mutation model.
        '''
        mutability_p, substitution_p = self.mutabilities(sequence)
        aa = ''.join([gencode[sequence[si:(si+3)]] for si in range(0, len(sequence), 3)])
        profile = [[0]*21 for fi in range(len(numb_profile))]
        S_cont = True
        i = 0
        while S_cont and (N is None or N > i):
            m = next(muts_iter)
            mut_seq = self.mutate(sequence, mutability_p.copy(), substitution_p.copy(), m=m)
            aa_mut = [gencode[mut_seq[si:(si+3)]] for si in range(0, len(mut_seq), 3)]
            for j, obs in enumerate(numb_profile):
                if obs == 0:
                    continue
                aa = aa_mut.pop(0)
                aa_idx = AA_INDEX[aa]
                profile[j][aa_idx] += 1
            i += 1
            if i % 100000 == 0:
                flat_profile = [pj for j, obs in enumerate(numb_profile) if obs != 0 for pj in profile[j]]
                missing_obs = len(flat_profile) - sum([(p>0)*1 for p in flat_profile])
                print('Still missing', missing_obs)
                if missing_obs < 50:
                    S_cont = True

        return profile


def sim_profile(row):
    cols = row.split(',')
    Nmuts_idx = args.profile_header.index('Nmuts')
    muts_iter = [int(m) for m in cols[Nmuts_idx].split(':')]
    muts_iter = cycle(muts_iter)  # <-- circular iterator
    seq_idx = args.profile_header.index('naive')
    sequence = cols[seq_idx]
    assert('*' not in [gencode[sequence[si:(si+3)]] for si in range(0, len(sequence), 3)])
    aa = ''.join([gencode[sequence[si:(si+3)]] for si in range(0, len(sequence), 3)])
    sequence_i = [['naive_seq', aa]]
    numbers_i = run_anarci(sequence_i)
    assert(numbers_i is not None)
    profile = convert2profile(numbers_i[0])
    numb_profile = [sum(pp) for pp in profile]
    clID_idx = args.profile_header.index('clusterID')
    if len(aa) == sum(numb_profile):
        mutation_model = MutationModel(sequence, args.mutability, args.substitution)
        profile = [cols[clID_idx], mutation_model.simulate_AAprofile(sequence, numb_profile, muts_iter, N=args.N, S=args.S, verbose=args.verbose)]
    else:
        profile = [cols[clID_idx], False]
    return profile


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Simulate an amino acid substitution profile given a starting sequence,'
                                                 'a substituion model and the number of mutations.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mutability', type=str, help='Path to mutability model file in S5F like format.')
    parser.add_argument('substitution', type=str, help='Path to substitution model file, in S5F like format.')
    parser.add_argument('--sequence', type=str, required=False, default=None, metavar="[ATGC]", help='Seed naive nucleotide sequence for the clonal family.')
    parser.add_argument('--N', type=int, required=False, default=None, help='Simulation size.')
    parser.add_argument('--nproc', type=int, required=False, default=1, help='Number of processes to start.')
    parser.add_argument('--S', type=int, required=False, default=None, help='Stopping criterium.')  ####### To be implemented e.g. 1) after all substitution have been observed, 2) after some convergence criterium
    parser.add_argument('--m', type=int, required=False, default=None, help='Number of substitutions to use in the simulation.')
    parser.add_argument('--m_file', type=str, required=False, help='Path to file with comma separated substituions for each clone from the clonal family.')
    parser.add_argument('--profile_file', type=str, required=False, help='Path to file with observed count profiles to calculate expected profiles for.')  #### To be implemented by looping through the list of profiles
    parser.add_argument('--verbose', type=bool, default=False, help='Print progress during simulation.')
    parser.add_argument('--outfile', type=str, required=False, help='Path to output file for simulated count profiles.')
    global args
    args = parser.parse_args()

    if args.profile_file:
        with open(args.profile_file) as fh:
            header = fh.readline()
            setattr(args, 'profile_header', header.strip().split(','))
            rows = fh.readlines()
        # Paralellise the process:
        pool = multiprocessing.Pool(args.nproc)
        profiles = pool.map(sim_profile, rows)
        # profiles = map(sim_profile, rows)  # No multiprocessing
    elif args.sequence is not None:
        assert('*' not in [gencode[args.sequence[si:(si+3)]] for si in range(0, len(args.sequence), 3)])

        muts_iter = list()
        if args.m_file is not None:
            with open(args.m_file) as fh:
                muts_iter = [int(m) for m in fh.readline().split(',')]
        elif args.m:
            try:
                assert(int(args.m) > 0)
                muts_iter.append(int(args.m))
            except:
                raise RuntimeError('Is the m parameter a positive interger')
        else:
            raise RuntimeError('Either "--m" or "--m_file" argument needs to be provided.')
        muts_iter = cycle(muts_iter)  # <-- circular iterator

        aa = ''.join([gencode[args.sequence[si:(si+3)]] for si in range(0, len(args.sequence), 3)])
        sequence_i = [['naive_seq', aa]]
        numbers_i = run_anarci(sequence_i)
        assert(numbers_i is not None)
        profile = convert2profile(numbers_i[0])
        numb_profile = [sum(pp) for pp in profile]
        assert(len(aa) == sum(numb_profile))
        mutation_model = MutationModel(args.sequence, args.mutability, args.substitution)
        profile = mutation_model.simulate_AAprofile(args.sequence, numb_profile, muts_iter, N=args.N, S=args.S, verbose=args.verbose)
    else:
        raise RuntimeError('Either --sequence or --profile_file argument must be provided.')

    if args.outfile:
        with open(args.outfile, 'w') as fhout:
            for p in profiles:
                if p[1] is False:
                    p[1] = [[0]*21 for fi in range(AHO_L)]
                pf = [p[0]] + [j for i in p[1] for j in i]
                ps = ','.join(map(str, pf))
                print(ps, file=fhout)


if __name__ == '__main__':
    main()

