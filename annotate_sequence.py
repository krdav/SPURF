#!/usr/bin/env python

import sys, os, glob, csv, random, copy, time, shutil, pickle
import gc
from itertools import cycle
csv.field_size_limit(sys.maxsize)
from anarci import anarci
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
global PATH2FILE
PATH2FILE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1, PATH2FILE)
from neutral_profile import MutationModel
partis_path = PATH2FILE + '/partis'
sys.path.insert(1, partis_path + '/python')
import utils
import glutils

# Notice a gap is added as the 21th amino acid:
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
AA_INDEX = {aa:i for i, aa in enumerate(AA_LIST)}
AHO_L = 149



def hamming_dist(seq1, seq2):
    '''Hamming distance between two sequences of equal length'''
    return sum(x != y for x, y in zip(seq1, seq2))


def repair_seq(seq, naiveDNA):
    # Convert to mutable:
    naiveDNA = list(naiveDNA)
    trim_seq = list(seq)
    # Repair all Ns from N padding or ambiguous bases:
    for i in range(len(naiveDNA)):
        if trim_seq[i] == 'N':
            trim_seq[i] = naiveDNA[i]
    assert('N' not in trim_seq)
    # Join to string and return:
    return ''.join(trim_seq)


def run_partis(seq):
    '''
    Infer VDJ genes and the naive sequence using partis.
    '''
    # Specify filenames:
    pretty_random_fnam = str(random.randint(1, 10**100))
    inpf = pretty_random_fnam + '_input'
    outf = pretty_random_fnam + '_output'
    # Write input fasta file for partis:
    with open(TMPDIR+'/'+inpf+'.fa', 'w') as fho:
        fho.write('>{}\n{}\n'.format('input_sequence', seq))
    # Run partis:
    cmd = '{}/bin/partis annotate --locus {} --species {} --infname {}/{}.fa --outfname {}/{}.csv'.format(partis_path, args.LOCUS, args.SPECIES, TMPDIR, inpf, TMPDIR, outf)
    os.system('{} > {}/{}.log'.format(cmd, TMPDIR, pretty_random_fnam))

    try:
        # Read the partis output file and extract the naive sequence:
        with open(TMPDIR+'/'+outf+'.csv') as fh:
            reader = csv.DictReader(fh)
            data = list(reader)
        ann = data[0]
        # Extract germline bounds info and trim the naive DNA sequence:
        try:
            utils.process_input_line(ann)       # Process dataframe row
            utils.add_implicit_info(glfo, ann)  # Adding germline infor
        except Exception as e:
            print e
            raise e

        if ann['stops'] is True:
            raise Exception('Input sequence contain stop codon. This is no valid.')
        elif ann['v_5p_del'] > 30 or ann['j_3p_del'] > 12:
            raise Exception('Incomplete input sequence error. 5-prime end missing {} nt and 3-prime missing {} nt. Max allowed is 30 and 12, respectively.'.format(ann['v_5p_del'], ann['j_3p_del']))
        elif ann['indelfos'][0]['indels']:
            raise Exception('Input sequence contains indels, this is currently not supported.')

        # Extract full size VDJ sequence for both the inferred naive and the input:
        full_gl_v = glfo['seqs']['v'][ann['v_gene']]  # Germline V
        full_gl_j = glfo['seqs']['j'][ann['j_gene']]  # Germline J

        gl_v_5p_del = full_gl_v[:ann['v_5p_del']]                      # 5-prime not included in input
        gl_j_3p_del = full_gl_j[(len(full_gl_j) - ann['j_3p_del']):]   # 3-prime not included in input
        assert full_gl_v[ann['v_5p_del']:] == ann['v_gl_seq']
        naiveDNA = gl_v_5p_del + ann['naive_seq'] + gl_j_3p_del  # Add the missing positions
        full_input_seq = 'N' * ann['v_5p_del'] + ann['input_seqs'][0] + 'N' * ann['j_3p_del']  # N pad the input sequence
        assert(len(naiveDNA) == len(full_input_seq))

        # Remove the untranslated end:
        if len(naiveDNA)%3 != 0:
            naiveDNA = naiveDNA[0:-(len(naiveDNA)%3)]
        if len(full_input_seq)%3 != 0:
             full_input_seq = full_input_seq[0:-(len(seq)%3)]
        if len(naiveDNA) != len(full_input_seq):
            raise Exception('Sequences not equally long after trimming.\nInput: {}\nNaive: {}\n.'.format(full_input_seq, naiveDNA))

        # Replace Ns in input sequence with naive DNA bases:
        full_input_seq = repair_seq(full_input_seq, naiveDNA[:])

        # If the inferred naive sequence contains a stop codon replace it by the input sequence codon:
        if '*' in str(Seq(naiveDNA, generic_dna).translate()):
            print 'Found stop codon in inferred naive sequnce, will replace with input sequence codon.'
            print 'Before replacement:', naiveDNA
            naiveDNA_l = list(naiveDNA[:])
            for codon in range(0, len(naiveDNA), 3):
                if '*' == str(Seq(naiveDNA[codon:codon+3], generic_dna).translate()):
                    naiveDNA_l[codon:codon+3] = full_input_seq[codon:codon+3]
            naiveDNA = ''.join(naiveDNA_l)
            print 'After replacement:', naiveDNA
        if '*' in str(Seq(naiveDNA, generic_dna).translate()):
            raise Exception('Naive sequence could not be repaired.')
        if naiveDNA == full_input_seq:
            print 'Warning: input sequence is identical to the inferred naive sequence.'
    finally:
        # Clean up:
        os.system('rm -r {}/{}* _output/*{}*'.format(TMPDIR, pretty_random_fnam, pretty_random_fnam))
    return(naiveDNA, full_input_seq, (ann['v_gene'], ann['d_gene'], ann['j_gene']))


def run_anarci(sequences):
    '''Run ANARCI annotation to get AHo numbering.'''
    allowed_species = 'human'
    allow = 'H'
    ncpu = 1
    scheme = 'aho'
    numbered, alignment_details, hit_tables = anarci(sequences, scheme=scheme, output=False, allow=allow, ncpu=ncpu, allowed_species=allowed_species)
    return numbered


def simulate_profile(muts, naiveDNA, numb_profile, mutability, substitution):
    '''
    Make a simulation under a neutral S5F motif model
    to make the expected neutral substitution profile
    for a given naive sequence and mutation burden.
    '''
    mutation_model = MutationModel(naiveDNA, mutability, substitution)
    muts_iter = cycle(muts)  # Cycle through the list of mutations in the input
    profile = mutation_model.simulate_AAprofile(naiveDNA, numb_profile, muts_iter, N=args.SIM_SIZE, S=None, verbose=True)
    gap = [0]*20 + [args.SIM_SIZE]
    profile = [pos if sum(pos) > 0 else gap for pos in profile]
    return profile


def make_dataframe(input_p, naive_p, neut_p, VDJ):
    '''Make an easy to print dataframe.'''
    header = ['Nseqs', 'v_gene', 'd_gene', 'j_gene']
    # Flatten AHo numbers:
    profile_header = ['p_{}_a_{}'.format(i, j) for i in range(1, 150) for j in range(1, 22)]
    header.extend(profile_header)
    df = [header]

    # Input sequence:
    flat_profile_input = [ai for si in input_p for ai in si]
    cols = [1]
    cols.extend(VDJ)
    cols.extend(flat_profile_input)
    df.append(cols)

    # Naive sequence:
    flat_profile_naive = [ai for si in naive_p for ai in si]
    cols = [1]
    cols.extend(VDJ)
    cols.extend(flat_profile_naive)
    df.append(cols)

    # Neutral profile:
    flat_profile_neut = [ai for si in neut_p for ai in si]
    cols = [args.SIM_SIZE]
    cols.extend(VDJ)
    cols.extend(flat_profile_neut)
    df.append(cols)
    return(df)


def AHo_annotate_naive(naiveAA):
    naive_list = [('naiveAA', naiveAA)]
    AHo_out = run_anarci(naive_list)

    if None in AHo_out or len(AHo_out) != 1:
        raise Exception('AHo numbering failed. Here is the output from ANARCI:', AHo_out)

    # Initialize the empty profile:
    profile = [[0]*21 for fi in range(AHO_L)]
    # Count each position in all input sequences:
    numbers = AHo_out[0][0]
    npop = numbers[0][:] + [[[-1]]]       # Add dummy to the end
    AAseq_pop = list(naiveAA[:]) + ['END']  # Add dummy to the end
    if len(AAseq_pop) != len(npop):
        Exception('AHo numbering failed. Here is the output from ANARCI:', AHo_out)
    ai = npop.pop(0)
    AA = AAseq_pop.pop(0)
    previous_p = -1
    AHo_seq = ''
    numb_profile = list()
    # Counts for each AHo position:
    for p in range(AHO_L):
        i = int(ai[0][0]) - 1
        # Weirdness in the junction region can make positions
        # assigned multiple times with insertion letters:
        if i == previous_p:
            Exception('AHo numbering failed. Looks like insertion numbers which are not supported. Here is the output from ANARCI:', AHo_out)
        else:
            previous_p = p
        if i == p:  # If the sequence has this AHo position defined assign the amino acid
            aa_idx = AA_INDEX[AA]
            AHo_seq += AA
            # Move to the next position:
            ai = npop.pop(0)
            AA = AAseq_pop.pop(0)
            numb_profile.append(1)
        else:  # If the sequence doesn't have this AHo position defined assign a gap character
            aa_idx = AA_INDEX['-']
            AHo_seq += '-'
            numb_profile.append(0)
        profile[p][aa_idx] += 1  # Add the observation count
    # Break out of the loop if the numbering failed:
    if len(AHo_seq) != AHO_L:
        raise Exception('len(AHo_seq) != AHO_L\nAHo_seq:\n{}\nExiting...'.format(AHo_seq))
    assert(len(npop) == 0)
    assert(len(AAseq_pop) == 0)
    assert(len(numb_profile) == AHO_L and sum(numb_profile) == len(naiveAA))
    return(profile, numb_profile)


def AHo_annotate_input(inputAA, numb_profile):
    AHo_input = [[0]*21 for fi in range(len(numb_profile))]
    aa_list = list(inputAA)
    for j, obs in enumerate(numb_profile):
        if obs == 0:
            aa = '-'
        else:
            aa = aa_list.pop(0)
        aa_idx = AA_INDEX[aa]
        AHo_input[j][aa_idx] = 1
    return(AHo_input)


def write_dataframe(df, outfile):
    fh_out = open(outfile, 'w')
    for row in df:
        fh_out.write(','.join(list(map(str, row))))
        fh_out.write('\n')
    fh_out.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Annotate BCR sequence for AMP.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--sequence', type=str, required=True, help='Sequence for annotation.')
    parser.add_argument('--outfile', type=str, default='out.csv', help='Output csv filename.')
    parser.add_argument('--SIM_SIZE', type=int, required=False, default=10000, help='Number of random draws to simulate the neutral profile.')
    parser.add_argument('--LOCUS', type=str, required=False, default='igh', help='Locus, either igh, igk or igl.')
    parser.add_argument('--SPECIES', type=str, required=False, default='human', help='Species, either human.')
    global args
    args = parser.parse_args()
    mutability = PATH2FILE + '/S5F/Mutability.csv'
    substitution = PATH2FILE + '/S5F/Substitution.csv'

    # Read default germline info:
    global glfo
    glfo = glutils.read_glfo(partis_path + '/data/germlines/human', locus=args.LOCUS)

    naive, fixed_input_seq, VDJ = run_partis(args.sequence)
    naiveAA = str(Seq(naive, generic_dna).translate())
    fixed_input_seqAA = str(Seq(fixed_input_seq, generic_dna).translate())

    # AHo annotate on the naive amino acid sequence:
    AHo_naive, numb_profile = AHo_annotate_naive(naiveAA)
    print(numb_profile)

    # Use the AHo annotation to make a profile over the input sequence:
    AHo_input = AHo_annotate_input(fixed_input_seqAA, numb_profile)

    # Simulate a profile under a neutral substitution process:
    Nmuts = hamming_dist(naive, fixed_input_seq)
    sim_profile = simulate_profile([Nmuts], naive, numb_profile, mutability, substitution)
    print(sim_profile)

    df = make_dataframe(AHo_input, AHo_naive, sim_profile, VDJ)
    write_dataframe(df, args.outfile)    


if __name__ == '__main__':
    # Make a tmp dir to dump crap:
    pretty_random_fnam = str(random.randint(1, 10**100))
    global TMPDIR
    TMPDIR = '/tmp/kd_tmp_' + pretty_random_fnam
    os.mkdir(TMPDIR)
    try:
        main()
    finally:
        shutil.rmtree(TMPDIR)  # rm -rf tmp dir

