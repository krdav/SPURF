#!/usr/bin/env python

import sys, os, glob, csv, random, copy, time, shutil, pickle
import gc
# random.seed(666)
from itertools import cycle
csv.field_size_limit(sys.maxsize)
from anarci import anarci
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
partis_path = '/fh/fast/matsen_e/kdavidse/partis'
sys.path.insert(1, partis_path + '/python')
import utils
import glutils
sys.path.insert(1, '/fh/fast/matsen_e/kdavidse/aammp')
from simulate_AAprofile_opti import MutationModel

# Notice a gap is added as the 21th amino acid:
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
AA_INDEX = {aa:i for i, aa in enumerate(AA_LIST)}
AHO_L = 149


class BadNaive(Exception):
    '''When the naive sequence is bad.'''


class CustomCry(Exception):
    '''Print exception.'''


def Most_Common(lst):
    '''Return the most used list element.'''
    data = Counter(lst)
    return data.most_common(1)[0][0]


def repair_seq(seq, naiveDNA, vj_bounds, keep_check=False):
    '''
    Repair unobserved sequence with the naive sequence, codon by codon,
    then strip the Ns if the sequence is still N padded and return the in-frame sequence.
    This function also checks the validity of a sequence if the keep_check option is set true.
    '''
    naiveDNA = naiveDNA[vj_bounds[0]:vj_bounds[1]]  # Slice from the beginning of V to the end of J
    # Remove the untranslated end:
    if len(naiveDNA)%3 != 0:
        naiveDNA = naiveDNA[0:-(len(naiveDNA)%3)]
    # Same the the input sequence:
    seq = seq[vj_bounds[0]:vj_bounds[1]]
    if len(seq)%3 != 0:
        seq = seq[0:-(len(seq)%3)]
    assert(len(seq) == len(naiveDNA))

    # Convert to mutable:
    naiveDNA = list(naiveDNA)
    trim_seq = list(seq)
    # Repair codon by codon, starting by the 3-prime end.
    # Reverse the codons to the the 3-prime end:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    # Repair until no more consecutive Ns from N padding:
    for i in range(3, args.MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
        else:
            break
    # Back reverse again, then trim the 5-prime end the same way:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    for i in range(3, args.MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
    # Join to string and check for inconsistencies:
    trim_seq = ''.join(trim_seq)
    # If in "check for keeping" mode, return either true or false:
    if keep_check is True:
        if 'N' in trim_seq or '*' in str(Seq(trim_seq, generic_dna).translate()) or 3*args.MIN_LEN > len(trim_seq):
            return False
        else:
            return True
    # Else expect no inconsistencies and return trimmed sequence if none:
    if 'N' in trim_seq or '*' in str(Seq(trim_seq, generic_dna).translate()) or 3*args.MIN_LEN > len(trim_seq):
        raise BadNaive('Problem with the naive sequence. Either contains N or stop codon, or is too short: {}'.format(trim_seq))
    else:
        return trim_seq


def repair_new_naive(seq, naiveDNA, vj_bounds):
    '''
    Check for consistencies in the re-inferred partis naive sequences.
    Most functions similar to the "repair_seq" function.
    '''
    naiveDNA = naiveDNA[vj_bounds[0]:vj_bounds[1]]
    if len(naiveDNA)%3 != 0:
        naiveDNA = naiveDNA[0:-(len(naiveDNA)%3)]

    seq = seq[vj_bounds[0]:vj_bounds[1]]
    if len(seq)%3 != 0:
        seq = seq[0:-(len(seq)%3)]  # Remove the untranslated end
    if len(seq) != len(naiveDNA) and len(seq) < vj_bounds[1]:
        raise BadNaive('Sequences not equally long when trimming.\nInput: {}\nNaive: {}\nBounds {}:{}\nIt appears that the input sequences is shorter than the bounds.'.format(seq, naiveDNA, str(vj_bounds[0]), str(vj_bounds[1])))
    elif len(seq) != len(naiveDNA):
        raise BadNaive('Sequences not equally long when trimming.\nInput: {}\nNaive: {}\nBounds {}:{}'.format(seq, naiveDNA, str(vj_bounds[0]), str(vj_bounds[1])))

    # Convert to mutable:
    naiveDNA = list(naiveDNA)
    trim_seq = list(seq)
    # Repair codon by codon, starting by the 3-prime end.
    # Reverse the codons to the the 3-prime end:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    # Repair until no more consecutive Ns from N padding:
    for i in range(3, args.MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
    # Back reverse again, then trim the 5-prime end the same way:
    naiveDNA = naiveDNA[::-1]
    trim_seq = trim_seq[::-1]
    for i in range(3, args.MAX_REP, 3):
        if 'N' in trim_seq[i-3:i]:
            trim_seq[i-3:i] = naiveDNA[i-3:i]
        else:
            break
    # Join to string and return:
    return ''.join(trim_seq)


def partis_naive_seq(lseq, fnam):
    '''
    Given a number of sequences infer the naive sequence using partis.
    '''
    # Specify filenames:
    pretty_random_fnam = str(random.randint(1, 10**100))
    inpf = pretty_random_fnam + '_input'
    outf = pretty_random_fnam + '_output'
    # Write input fasta file for partis:
    with open(TMPDIR+'/'+inpf+'.fa', 'w') as fho:
        for i, s in enumerate(lseq):
            fho.write('>{}\n{}\n'.format(str(i), s))
    # Run partis:
    cmd = '{}/bin/partis partition --locus {} --species {} --infname {}/{}.fa --outfname {}/{}.csv'.format(partis_path, args.LOCUS, args.SPECIES, TMPDIR, inpf, TMPDIR, outf)
    # os.system(cmd)  # Print partis STDOUT to screen
    os.system('{} > {}/{}.log'.format(cmd, TMPDIR, pretty_random_fnam))

    try:
        # Read the partis output file and extract the naive sequence:
        with open(TMPDIR+'/'+outf+'-cluster-annotations.csv') as fh:
            reader = csv.DictReader(fh)
            data = list(reader)
        # assert(len(data) == 1)  # There should really only be one clonal family, but there often are, so just take the first (largest)
        # Extract germline bounds info and trim the naive DNA sequence:
        try:
            utils.process_input_line(data[0])       # Process dataframe row
            fnam_base = fnam.split('_partitions')[0].split('/')
            #glfo = glutils.read_glfo('{}/_output/{}/hmm/germline-sets'.format(fnam_base[0], fnam_base[-1]), locus=args.LOCUS)
            glfo = glutils.read_glfo(partis_path + '/data/germlines/human', locus=args.LOCUS)
            utils.add_implicit_info(glfo, data[0])  # Adding germline infor
        except Exception as e:
            print e
            raise e

        naiveDNA = data[0]['naive_seq'][:]
        first_lseq = data[0]['input_seqs'][:][0]
        vj_bounds = (data[0]['regional_bounds']['v'][0], data[0]['regional_bounds']['j'][1])
        naiveDNA = repair_new_naive(naiveDNA[:], naiveDNA[:], vj_bounds)
        first_lseq = repair_new_naive(first_lseq, naiveDNA[:], vj_bounds)
        try:
            assert(len(first_lseq) == len(naiveDNA))
        except:
            print 'len(first_lseq) != len(data[0]["naive_seq"])'
            print len(first_lseq)
            print first_lseq
            print len(naiveDNA)
            print naiveDNA
        # If the inferred naive sequence contains a stop codon replace it by the input sequence codon:
        if '*' in str(Seq(naiveDNA, generic_dna).translate()):
            print 'Found stop codon in inferred naive sequnce, will replace with input sequence codon.'
            print 'Before replacement:', naiveDNA
            naiveDNA_l = list(naiveDNA[:])
            for codon in range(vj_bounds[0], vj_bounds[1], 3):
                if '*' == str(Seq(naiveDNA[codon:codon+3], generic_dna).translate()):
                    naiveDNA_l[codon:codon+3] = first_lseq[codon:codon+3]
            naiveDNA = ''.join(naiveDNA_l)
            print 'After replacement:', naiveDNA
        if naiveDNA == first_lseq:
            print 'Complaining to say naiveDNA == first_lseq (nothing bad just to be sure the repair is not just replacing the naive sequence with the input entirely)'

        return(naiveDNA)
    finally:
        # Clean up:
        os.system('rm -r {}/{}* _output/*{}*'.format(TMPDIR, pretty_random_fnam, pretty_random_fnam))


def extract_seqs(fnam):
    '''
    Reads a partis cluster-annotations file and extracts relevant information and sequences.
    '''
    # Read cluster annotations into a data list of dictionaries:
    with open(fnam) as fh:
        reader = csv.DictReader(fh)
        data = list(reader)

    sequences_i = list()
    info_i = list()

    if args.allele_finding:
        fnam_base = fnam.split('_partitions')[0].split('/')
        glfo = glutils.read_glfo('{}/_output/{}/hmm/germline-sets'.format(fnam_base[0], fnam_base[-1]), locus=args.LOCUS)
    else:
        glfo = glutils.read_glfo(partis_path + '/data/germlines/human', locus=args.LOCUS)
    for row in data:
        # Process the partis data row and add germline information:
        try:
            utils.process_input_line(row)
            # Read default germline info
            utils.add_implicit_info(glfo, row)
        except Exception as e:  # Skip rows that cannot be processed
            if 'failed annotation' not in e:
                pass
                # print('First skip')
                # print(e)
            else:
                print 'Reading from'
                print '{}/_output/{}/hmm/germline-sets'.format(fnam_base[0], fnam_base[-1])
                print e
            continue

#        # Process the partis data row and add germline information:
#        try:
#            utils.process_input_line(row)
#            utils.add_implicit_info(glfo, row)
#        except:  # Skip rows that cannot be processed
#            continue

        # Extract the full N padded naive sequence,
        # and find the v -and j gene bound on this naive sequence:
        cdr3_bounds = (row['codon_positions']['v'], row['codon_positions']['j'] + 3)
        vj_bounds = (row['regional_bounds']['v'][0], row['regional_bounds']['j'][1])
        naiveDNA = row['naive_seq']
        # Skip naive sequences too short or with stop codons:
        if repair_seq(naiveDNA, naiveDNA, vj_bounds, keep_check=True) is False:
            continue
        trimmed_naiveDNA = repair_seq(naiveDNA[:], naiveDNA[:], vj_bounds)
        naiveAA = str(Seq(trimmed_naiveDNA, generic_dna).translate())

        # There has been a name change and this try/except is meant to provide backwards compatability:
        try:
            lseq = row['input_seqs'][:]
        except:
            lseq = row['seqs'][:]
        ir_lseq = row['indel_reversed_seqs']
        stop_seq = row['stops']
        assert(len(lseq) == len(ir_lseq))
        assert(len(lseq) == len(stop_seq))
        # Only keep sequences without indels and stop codons and minimum length amino acid length (QC):
        ### ir_lseq[i] == '' or lseq[i] == ir_lseq[i]  <-- No indels
        ### stop_seq[i]  <-- No partis annotated stops (there seems still to be stops after these are removed though)
        ### repair_seq(lseq[i], naiveDNA, vj_bounds, keep_check=True)  <-- Checks whether the sequence is long enougth or have stop codons
        keep_idx = [1 if ((ir_lseq[i] == '' or lseq[i] == ir_lseq[i]) and stop_seq[i] is False and repair_seq(lseq[i], naiveDNA, vj_bounds, keep_check=True)) else 0 for i in range(len(lseq))]

        # Now only keep those sequences that passed QC:
        lseq = [s for s, keep in zip(lseq, keep_idx) if keep == 1]
        # Get amino acid sequences:
        lAAseq = [str(Seq(repair_seq(s[:], naiveDNA[:], vj_bounds), generic_dna).translate()) for s in lseq]
        # And mutation frequencies:
        mut_freqs = [s for s, keep in zip(row['mut_freqs'], keep_idx) if keep == 1]
        assert(len(mut_freqs) == len(lseq))
        # Convert frequency to counts:
        Nmuts = [int(round(float(t[0])*len(t[1].strip('N')))) for i, t in enumerate(zip(mut_freqs, lseq))]

        # Deduplicate AAseqs and lseq according to the duplications on amino acid level:
        lAAseq_dict = dict()
        lseq_unique = list()
        for i, aa in enumerate(lAAseq):
            if aa in lAAseq_dict:
                lAAseq_dict[aa].append(i)
            else:
                lAAseq_dict[aa] = [i]
                lseq_unique.append(repair_seq(lseq[i][:], naiveDNA[:], vj_bounds))
        assert(len(lAAseq_dict) == len(lseq_unique))
        # Make the deduplicated sequence list and the mutation rates:
        lAAseq_dedup = list()
        Nmuts_dedup = list()
        for aa, idxs in lAAseq_dict.items():
            lAAseq_dedup.append(aa)
            Nmut_list = [float(Nmuts[i]) for i in idxs]
            Nmuts_dedup.append(int(round(sum(Nmut_list)/len(Nmut_list))))
        assert(len(lAAseq_dedup) == len(Nmuts_dedup))
        assert(len(lAAseq_dedup) == len(lseq_unique))

        # Exclude small clonal families after all the QC and deduplication:
        if len(lAAseq_dedup) < args.MIN_OBS:
            continue

        # Store the results in a list:
        sequences_i.append(['naive_seq', naiveAA])  # This format is for ANARCI numbering
        info_i.append({'fnam': fnam, 'v_gene': row['v_gene'], 'd_gene': row['d_gene'], 'j_gene': row['j_gene'],
                       'naive_seq': naiveAA, 'naive_seq_DNA': trimmed_naiveDNA, 'Nmuts': Nmuts_dedup[:],
                       'AAseqs': lAAseq_dedup[:], 'DNAseqs': lseq_unique[:]})
    return(sequences_i, info_i)


def run_anarci(sequences):
    '''Run ANARCI annotation to get AHo numbering.'''
    allowed_species = 'human'
    allow = 'H'
    ncpu = 1
    scheme = 'aho'
    numbered, alignment_details, hit_tables = anarci(sequences, scheme=scheme, output=False, allow=allow, ncpu=ncpu, allowed_species=allowed_species)
    return numbered


def flat_counts(all_numbered, info_i_j):
    '''Create a flat AHo profile of counts.'''
    assert(len(all_numbered) == len(info_i_j))
    profiles = list()
    fail_count = 0
    ii = -1  # Counter for look-up in info_i_j
    for numbers_i in all_numbered[:]:  # All clonal families
        ii += 1
        # Check that the numbering went successfully:
        if numbers_i is None:
            fail_count += 1
            del info_i_j[ii]
            del all_numbered[ii]
            ii -= 1
            continue
        numbers = numbers_i[0][:]
        if numbers is None or len(numbers[0]) != len(info_i_j[ii]['naive_seq']):
            fail_count += 1
            del info_i_j[ii]
            del all_numbered[ii]
            ii -= 1
            continue

        # Initialize the empty profile:
        profile = [[0]*21 for fi in range(AHO_L)]
        AHo_align = list()
        # Count each position in all input sequences:
        for AAseq in info_i_j[ii]['AAseqs'][:]:
            assert('*' not in AAseq)
            npop = numbers[0][:] + [[[-1]]]       # Add dummy to the end
            AAseq_pop = list(AAseq[:]) + ['END']  # Add dummy to the end
            try:
                assert(len(AAseq_pop) == len(npop))
            except:
                print 'assert(len(AAseq_pop) == len(npop))'
                print len(info_i_j[ii]['AAseqs'])
                print len(AAseq_pop)
                print len(npop)
                # print info_i_j[ii]['AAseqs']
                sys.exit()
            ai = npop.pop(0)
            AA = AAseq_pop.pop(0)
            previous_p = -1
            failed_numbering = False
            AHo_seq = ''
            # Counts for each AHo position:
            for p in range(AHO_L):
                i = int(ai[0][0]) - 1
                # Weirdness in the junction region can make positions
                # assigned multiple times with insertion letters:
                if i == previous_p:
                    failed_numbering = True
                    break
                else:
                    previous_p = p
                if i == p:  # If the sequence has this AHo position defined assign the amino acid
                    aa_idx = AA_INDEX[AA]
                    AHo_seq += AA
                    # Move to the next position:
                    ai = npop.pop(0)
                    AA = AAseq_pop.pop(0)
                else:  # If the sequence doesn't have this AHo position defined assign a gap character
                    aa_idx = AA_INDEX['-']
                    AHo_seq += '-'
                profile[p][aa_idx] += 1  # Add the observation count
            # Break out of the loop if the numbering failed:
            if failed_numbering is True:
                break
            if len(AHo_seq) != AHO_L:
                raise CustomCry('len(AHo_seq) != AHO_L\nAHo_seq:\n{}\nExiting...'.format(AHo_seq))
            # Check that both sequence and numbering is empty:
            assert(len(npop) == 0)
            assert(len(AAseq_pop) == 0)
            AHo_align.append(AHo_seq)

        # Delete entry if wrongly numbered:
        if failed_numbering is True:
            fail_count += 1
            del info_i_j[ii]
            del all_numbered[ii]
            ii -= 1
            continue

        nseq = len(info_i_j[ii]['AAseqs'])
        # Now validate the profile:
        assert(sum([sum(p) for p in profile]) == AHO_L*nseq)
        # Append results:
        profiles.append((profile[:], nseq))
        info_i_j[ii]['AHo_align'] = AHo_align[:]

    if fail_count:
        print 'There were clonal families that failed validation and were discared:', fail_count
    return profiles, info_i_j, all_numbered


def flat_counts_subsample_one(package):
    '''
    Create a flat AHo profile of counts for the subsamples.
    The numbered sequences have already gone through different checks in the flat_counts function,
    and therefore this is more loose in its assertions.
    '''
    ii, numbers_i, info, args = copy.deepcopy(package)
    print 'Making subsample', ii
    profiles = list()
    new_info = list()
    numbers = numbers_i[0][:]
    ### First make the full profiles
    # Initialize the empty profile:
    profile = [[0]*21 for fi in range(AHO_L)]
    for AAseq in info['AAseqs'][:]:
        numb_profile = list()
        npop = numbers[0][:] + [[[-1]]]       # Add dummy to the end
        AAseq_pop = list(AAseq[:]) + ['END']  # Add dummy to the end
        assert(len(AAseq_pop) == len(npop))
        ai = npop.pop(0)
        AA = AAseq_pop.pop(0)
        for p in range(AHO_L):
            i = int(ai[0][0]) - 1
            if i == p:  # If the sequence has this AHo position defined assign the amino acid
                aa_idx = AA_INDEX[AA]
                # Move to the next position:
                ai = npop.pop(0)
                AA = AAseq_pop.pop(0)
                numb_profile.append(1)
            else:  # If the sequence doesn't have this AHo position defined assign a gap character
                aa_idx = AA_INDEX['-']
                numb_profile.append(0)
            profile[p][aa_idx] += 1  # Add the observation count
        # Check that both sequence and numbering is empty:
        assert(len(npop) == 0)
        assert(len(AAseq_pop) == 0)
        assert(len(numb_profile) == AHO_L)

    nseq = len(info['AAseqs'])
    # Now validate the profile:
    assert(sum([sum(p) for p in profile]) == AHO_L*nseq)
    profiles.append((profile[:], nseq))

    # Add new info about the full sample:
    new_info.append(copy.deepcopy(info))
    new_info[-1]['idx'] = ii
    new_info[-1]['SubSize'] = nseq

    ### Then add the subsampled profiles:
    if args.totN_sub is None:
        nsub_seqs = int(args.fsub * nseq)
    else:
        nsub_seqs = args.totN_sub
    # Make a random subsample, nboots times:
    for nb in range(args.nboots):
        tries = 100
        for t in range(tries):
            # Make a copy of the full info:
            new_info.append(copy.deepcopy(info))
            # Initialize the profile:
            profile = [[0]*21 for fi in range(AHO_L)]
            # Make a random sample:
            rand_idx = random.sample(list(range(nseq)), nsub_seqs)
            try:
                # Subsample the full profile:
                new_info[-1]['AAseqs'] = [copy.deepcopy(info['AAseqs'][ri]) for ri in rand_idx]
                new_info[-1]['DNAseqs'] = [copy.deepcopy(info['DNAseqs'][ri]) for ri in rand_idx]
            except:
                print len(info['AAseqs'])
                print len(info['DNAseqs'])
                raise
            # Re-infer the naive sequence:
            try:
                naiveDNA = partis_naive_seq(new_info[-1]['DNAseqs'], new_info[-1]['fnam'])
                naiveAA = str(Seq(naiveDNA, generic_dna).translate())
                if len(naiveAA) != len(new_info[-1]['naive_seq']):
                    raise BadNaive('Problem with the naive sequence. '
                                   'The subsampled sequence has an inferred naive sequence of different length '
                                   'than its parent clonal family: \n{}\n'.format(naiveAA, new_info[-1]['naive_seq']))
                elif 'N' in naiveDNA:
                    raise BadNaive('Problem with the naive sequence. It contains N.')
            except BadNaive as e:
                print 'Got message:', e
                print 'Used this input sequence:', new_info[-1]['DNAseqs']
                print 'Bad naive sequence. Trying again.'
                new_info.pop()  # Remove the subsample of the bad naive sequence
                continue  # Try again if bad naive
            new_info[-1]['naive_seq'] = naiveAA
            new_info[-1]['naive_seq_DNA'] = naiveDNA
            assert(sum(numb_profile) == len(new_info[-1]['naive_seq']))

            # Make the profile from the random sample:
            for AAseq in new_info[-1]['AAseqs'][:]:
                npop = numbers[0][:] + [[[-1]]]       # Add dummy to the end
                AAseq_pop = list(AAseq[:]) + ['END']  # Add dummy to the end
                assert(len(AAseq_pop) == len(npop))
                ai = npop.pop(0)
                AA = AAseq_pop.pop(0)
                for p in range(AHO_L):
                    i = int(ai[0][0]) - 1
                    if i == p:  # If the sequence has this AHo position defined assign the amino acid
                        aa_idx = AA_INDEX[AA]
                        # Move to the next position:
                        ai = npop.pop(0)
                        AA = AAseq_pop.pop(0)
                    else:  # If the sequence doesn't have this AHo position defined assign a gap character
                        aa_idx = AA_INDEX['-']
                    profile[p][aa_idx] += 1  # Add the observation count
                # Check that both sequence and numbering is empty:
                assert(len(npop) == 0)
                assert(len(AAseq_pop) == 0)

            # The sequence specific information needs to be subsampled as well:
            new_info[-1]['idx'] = ii
            new_info[-1]['SubSize'] = nsub_seqs
            new_info[-1]['Nmuts'] = [copy.deepcopy(info['Nmuts'][ri]) for ri in rand_idx]
            # Add the subsampled profile:
            assert(sum([sum(p) for p in profile]) == AHO_L*nsub_seqs)
            profiles.append((profile[:], nseq))
            # Simulate the expected neutral evolution profile:
            sim_profile = simulate_profile(new_info[-1]['Nmuts'], new_info[-1]['naive_seq_DNA'], numb_profile)
            new_info.append(copy.deepcopy(new_info[-1]))
            new_info[-1]['SubSize'] = args.SIM_SIZE
            profiles.append((sim_profile[:], nseq))
            break  # If getting here break out of the loop
        if t == (tries - 1):
            raise CustomCry('Exceded maximum number of tries, quitting.')

    return profiles, new_info


def simulate_profile(muts, naiveDNA, numb_profile):
    '''
    Make a simulation under a neutral S5F motif model
    to make the expected neutral substitution profile
    for a given naive sequence and mutation burden.
    '''
    mutability = '/fh/fast/matsen_e/kdavidse/gctree/S5F/Mutability.csv'
    substitution = '/fh/fast/matsen_e/kdavidse/gctree/S5F/Substitution.csv'
    mutation_model = MutationModel(naiveDNA, mutability, substitution)
    muts_iter = cycle(muts)  # Cycle through the list of mutations in the input
    profile = mutation_model.simulate_AAprofile(naiveDNA, numb_profile, muts_iter, N=args.SIM_SIZE, S=None, verbose=True)
    return profile


def sub_par(all_numbered_trim_sub, info_i_j_sub, args):
    '''Submit the jobs to a job pool.'''
    assert(len(all_numbered_trim_sub) == len(info_i_j_sub))
    # Prep. data packages for pool:
    packages = [(i, all_numbered_trim_sub[i], info_i_j_sub[i], args) for i in range(len(all_numbered_trim_sub))]
    if args.nproc == 1:  # Run without subprocess
        results = map(flat_counts_subsample_one, packages)
    else:  # Paralellize the processes by forking
        import multiprocessing
        pool = multiprocessing.Pool(processes=args.nproc)  # Start the pool
        results = pool.map_async(flat_counts_subsample_one, packages, chunksize=1)  # Run subprocesses
        pool.close()
        pool.join()
        results = results.get()
    # Oneliner to expand and flatten the two lists to be returned:
    # results = [[[1...nboot], 1...nboot]... nprofiles]  <-- This is how results look
    return [[se for e in li for se in e] for li in zip(*results)]


def make_dataframe(profiles, info_i_j):
    '''Make an easy to print dataframe.'''
    assert(len(profiles) == len(info_i_j))
    header = ['clusterID', 'naiveAA', 'naive', 'Nseqs', 'v_gene', 'd_gene', 'j_gene', 'filename', 'Nmuts']
    # Flatten AHo numbers:
    profile_header = ['p_{}_a_{}'.format(i, j) for i in range(1, 150) for j in range(1, 22)]
    header.extend(profile_header)
    df = [header]
    for i, info in enumerate(info_i_j):
        cols = [i, info['naive_seq'], info['naive_seq_DNA'], profiles[i][1], info['v_gene'], info['d_gene'], info['j_gene'], info['fnam'], ':'.join(map(str, info['Nmuts']))]
        flat_profile = [ai for si in profiles[i][0] for ai in si]
        cols.extend(flat_profile)
        try:
            assert(len(cols) == len(header))
        except:
            raise CustomCry('More data columns than header elements.')
        df.append(cols)
    return df


def make_dataframe_sub(profiles, info_i_j):
    '''Make an easy to print dataframe.'''
    assert(len(profiles) == len(info_i_j))
    header = ['clusterID', 'naiveAA', 'naive', 'SubSize', 'Nseqs', 'v_gene', 'd_gene', 'j_gene', 'filename', 'Nmuts']
    # Flatten AHo numbers:
    profile_header = ['p_{}_a_{}'.format(i, j) for i in range(1, 150) for j in range(1, 22)]
    header.extend(profile_header)
    df = [header]
    for i, info in enumerate(info_i_j):
        cols = [info['idx'], info['naive_seq'], info['naive_seq_DNA'], info['SubSize'], profiles[i][1], info['v_gene'], info['d_gene'], info['j_gene'], info['fnam'], ':'.join(map(str, info['Nmuts']))]
        flat_profile = [ai for si in profiles[i][0] for ai in si]
        cols.extend(flat_profile)
        try:
            assert(len(cols) == len(header))
        except:
            raise CustomCry('More data columns than header elements.')
        df.append(cols)
    return df


def write_dataframe(df, outfile):
    fh_out = open(outfile, 'w')
    for row in df:
        fh_out.write(','.join(list(map(str, row))))
        fh_out.write('\n')
    fh_out.close()


def write_seqs(info_i_j, outfile):
    header = ['clusterID', 'AAseqs']
    fh_out = open(outfile, 'w')
    fh_out.write(','.join(list(map(str, header))))
    fh_out.write('\n')
    for i, info in enumerate(info_i_j):
        if 'idx' in info:
            cols = [info['idx'], ':'.join(info['AHo_align'])]
        else:
            cols = [i, ':'.join(info['AHo_align'])]
        fh_out.write(','.join(list(map(str, cols))))
        fh_out.write('\n')
    fh_out.close()


def write_seqs_sub(info_i_j, outfile):
    header = ['clusterID', 'AAseqs']
    fh_out = open(outfile, 'w')
    fh_out.write(','.join(list(map(str, header))))
    fh_out.write('\n')
    for i, info in enumerate(info_i_j):
        if 'idx' in info:
            cols = [info['idx'], ':'.join(info['AAseqs'])]
        else:
            cols = [i, ':'.join(info['AAseqs'])]
        fh_out.write(','.join(list(map(str, cols))))
        fh_out.write('\n')
    fh_out.close()


def run_file(f):
    '''
    Wrapping the sequence extraction etc. for each cluster file.
    This function can then be called outside by a pool of subprocesses.
    '''
    sequences_i, info_i = extract_seqs(f)
    assert(len(sequences_i) == len(info_i))
    print 'Reading {} containing {} GCs.'.format(f, len(sequences_i))
    if len(sequences_i) == 0:
        return (False, False, False)
    all_numbered = run_anarci(sequences_i)  # Run ANARCI on each naiveAA sequence
    profiles, info_i, all_numbered_trim = flat_counts(all_numbered, info_i)
    assert(len(all_numbered_trim) == len(info_i))
    assert(len(profiles) == len(info_i))
    if len(profiles) == 0:
        return (False, False, False)
    else:
        return (info_i, all_numbered_trim, profiles)


def make_profile(glob_res, saved_variables):
    if args.nproc == 1:
        results = map(run_file, glob_res)  # Run without subprocess
    else:
        import multiprocessing
        # Paralellize the process:
        pool = multiprocessing.Pool(args.nproc)  # Start the pool
        results = pool.map_async(run_file, glob_res, chunksize=1)  # Run subprocesses
        pool.close()
        pool.join()
        results = results.get()

    # Unpack and merge the results:
    info_i_j = list()           # GC data for j in file i extra information beloning to the sequences
    profiles = list()           # Count profiles
    all_numbered_trim = list()  # AHo numbered naive sequences
    for t in results:
        if False not in t:
            info_i, all_numbered_trim_i, profiles_i = t
            info_i_j.extend(info_i)
            all_numbered_trim.extend(all_numbered_trim_i)
            profiles.extend(profiles_i)

    # Collapse all clonal families with same naive amino acid sequnce:
    collapse_count = 0
    if args.collapse_clusters != None:
        if args.collapse_clusters == 'DNA':
            DNAorAA = 'naive_seq_DNA'
        elif args.collapse_clusters == 'AA':
            DNAorAA = 'naive_seq'
        else:
            raise CustomCry('Did not understand collapse_clusters input:', args.collapse_clusters)
        naiveAA_set = dict()
        cf = 0
        for cf_info in info_i_j[:]:
            if cf_info[DNAorAA] in naiveAA_set:  # Duplicated naiveAA sequence
                first_obs = naiveAA_set[cf_info[DNAorAA]]
                info_dup = info_i_j.pop(cf)  # Extract info
                profiles.pop(cf)             # Discard profile. This needs recalculation
                all_numbered_trim.pop(cf)    # Discard numbering since this is the same

                ### Update info for info_i_j[first_obs]
                # First extract all the info:
                lAAseq = info_i_j[first_obs]['AAseqs'][:] + info_dup['AAseqs'][:]
                lDNAseq = info_i_j[first_obs]['DNAseqs'][:] + info_dup['DNAseqs'][:]
                Nmuts = info_i_j[first_obs]['Nmuts'][:] + info_dup['Nmuts'][:]
                # Deduplicate AAseqs:
                lAAseq_dict = dict()
                lDNAseq_dedup = list()
                for i, aa in enumerate(lAAseq):
                    if aa in lAAseq_dict:
                        lAAseq_dict[aa].append(i)
                    else:
                        lAAseq_dict[aa] = [i]
                        lDNAseq_dedup.append(lDNAseq[i])
                # Make the deduplicated list and take the mutation rates,
                #  as the mutation rate for the deduplicated sequence:
                lAAseq_dedup = list()
                Nmuts_dedup = list()
                for aa, idxs in lAAseq_dict.items():
                    lAAseq_dedup.append(aa)
                    Nmut_list = [float(Nmuts[i]) for i in idxs]
                    Nmuts_dedup.append(int(round(sum(Nmut_list)/len(Nmut_list))))
                assert(len(lAAseq_dedup) == len(lDNAseq_dedup))
                assert(len(lAAseq_dedup) == len(Nmuts_dedup))
                assert(len(info_dup['AAseqs']) <= len(lAAseq_dedup))

                # Update the lists of info elements:
                info_i_j[first_obs]['AAseqs'] = lAAseq_dedup[:]
                info_i_j[first_obs]['DNAseqs'] = lDNAseq_dedup[:]
                info_i_j[first_obs]['Nmuts'] = Nmuts_dedup[:]

                # Rerun the profile calculation:
                one_profile, _, _ = flat_counts(all_numbered_trim[first_obs:(1+first_obs)], info_i_j[first_obs:(1+first_obs)])
                assert(len(one_profile) == 1)
                # Update the profile:
                profiles[first_obs] = one_profile.pop()
                collapse_count += 1
            else:
                naiveAA_set[cf_info[DNAorAA]] = cf
                cf += 1

    if collapse_count > 0:
        print 'Clonal families were collapsed due to identical {}: {}'.format(DNAorAA, collapse_count)
    print 'Total GCs in files:', len(profiles)
    print 'Total sequences in all GCs:', sum([len(ll['AAseqs']) for ll in info_i_j])

    # Subset the most frequent V gene:
    if args.highest_Vgene:
        # Find the most frequent V gene:
        vgenes = [info['v_gene'].split('*')[0] for info in info_i_j]
        MC_vgene = Most_Common(vgenes)
        # Slice the subset out:
        keep_idx = [idx for idx, info in enumerate(info_i_j) if info['v_gene'].split('*')[0] == MC_vgene]
        profiles = [profiles[i] for i in keep_idx]
        info_i_j = [info_i_j[i] for i in keep_idx]
        all_numbered_trim = [all_numbered_trim[i] for i in keep_idx]
        # Add extention to output names:
        outfile_p += '_' + MC_vgene
        outfile_s += '_' + MC_vgene
        outfile_sub_p += '_' + MC_vgene
        outfile_sub_s += '_' + MC_vgene

        print 'After subsetting to V gene', MC_vgene
        print 'Total GCs in files:', len(profiles)
        print 'Total sequences in all GCs:', sum([len(ll['AAseqs']) for ll in info_i_j])

    if len(profiles) < args.nsubs:
        raise CustomCry('Requested more subsamples than profiles. Lower the number of subsamples or add more clonal families to the data.')

    # Sort lists according to number of sequence in each clonal family:
    sort_idx = [t[0] for t in sorted(enumerate(info_i_j), key=lambda x: len(x[1]['AAseqs']), reverse=True)]
    profiles = [profiles[i] for i in sort_idx]
    info_i_j = [info_i_j[i] for i in sort_idx]
    all_numbered_trim = [all_numbered_trim[i] for i in sort_idx]

    min_nseq = len(info_i_j[args.nsubs-1]['AAseqs'])
    if args.totN_sub is None and (int(args.fsub * min_nseq) == 0 or args.minN_mixing > min_nseq):
        raise CustomCry('The {}th largest clonal family of sequences only has {} sequences and therefore this cannot be subsampled by a fraction of {}.'.format(args.nsubs, min_nseq, args.fsub))
    elif args.totN_sub is not None and (args.totN_sub > min_nseq or args.minN_mixing > min_nseq):
        raise CustomCry('The {}th largest clonal family of sequences only has {} sequences and therefore this cannot be subsampled by {} sequences.'.format(args.nsubs, min_nseq, args.totN_sub))

    # Enforce equal sampling from the different dataset files,
    # requiring a minimum of minN observations:
    if args.sub_mixing:
        # Make a list of list filenames and indices for clonal families with minN observations:
        allowed_subs = [(t[1]['fnam'], t[0]) for t in enumerate(info_i_j) if len(t[1]['AAseqs']) >= args.minN_mixing]
        # Convert to dict with list of indices:
        fnam_alleowed_subs = dict()
        for t in allowed_subs:
            if t[0] not in fnam_alleowed_subs:
                fnam_alleowed_subs[t[0]] = [t[1]]
            else:
                fnam_alleowed_subs[t[0]].append(t[1])
        subs_idxs = list()  # The list indices to use for subsampling
        fnams = cycle(list(fnam_alleowed_subs.keys()))  # <-- circular iterator of filenames
        j = 0
        while j < args.nsubs:   # Loop until all args.nsubs indices for the subsampled profiles are found
            fnam = next(fnams)  # Cycle through the filenames
            if fnam in fnam_alleowed_subs:
                sub_idx = fnam_alleowed_subs[fnam].pop(0)  # Pop from the top to get highest frequency first
                subs_idxs.append(sub_idx)
                # When no more allowed subsamples for a given file, delete the entry
                if len(fnam_alleowed_subs[fnam]) == 0:
                    del fnam_alleowed_subs[fnam]
                j += 1
        # Convert the list to a set for fast lookup:
        subs_idxs = set(subs_idxs)
        # Make the subsets:
        all_numbered_trim_sub = [el for i, el in enumerate(all_numbered_trim[:]) if i in subs_idxs]
        info_i_j_sub = [el for i, el in enumerate(info_i_j[:]) if i in subs_idxs]
        all_numbered_trim = [el for i, el in enumerate(all_numbered_trim[:]) if i not in subs_idxs]
        profiles = [el for i, el in enumerate(profiles[:]) if i not in subs_idxs]
        info_i_j = [el for i, el in enumerate(info_i_j[:]) if i not in subs_idxs]
    else:  # If no mixing just take the highest frequency profiles
        all_numbered_trim_sub = all_numbered_trim[0:args.nsubs]
        info_i_j_sub = info_i_j[0:args.nsubs]
        all_numbered_trim = all_numbered_trim[args.nsubs:]
        profiles = profiles[args.nsubs:]
        info_i_j = info_i_j[args.nsubs:]
    assert((len(info_i_j) + len(info_i_j_sub)) == len(sort_idx))
    assert(len(info_i_j_sub) == args.nsubs)

    # Dump important variables:
    with open(saved_variables, 'wb') as f:
        pickle.dump((profiles, all_numbered_trim, info_i_j, all_numbered_trim_sub, info_i_j_sub), f)
        print 'Dumped variables to pickle file:', saved_variables

    return profiles, all_numbered_trim, info_i_j, all_numbered_trim_sub, info_i_j_sub


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Create a BCR sequences profile from partis clonal family partitions.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--glob_pat', type=str, required=False, default='*cluster-annotations.csv', help='Glob pattern to find the partis "cluster annotation" files.')
    parser.add_argument('--nsubs', type=int, required=False, default=200, help='How many subsampled profiles should be taken out?')
    parser.add_argument('--fsub', type=float, required=False, default=None, help='Fraction of sequences in the full profile to include in the subsample.')
    parser.add_argument('--totN_sub', type=int, required=False, default=10, help='Total number of sequences in the full profile to include in the subsample (argmuent overlaps with --args.fsub).')
    parser.add_argument('--nboots', type=int, required=False, default=10, help='"Bootstrap" replicate subsamples of the same full profile.')
    parser.add_argument('--nproc', type=int, required=False, default=1, help='Number of processes to start.')
    parser.add_argument('--dataset_name', type=str, required=False, default='not_specified', help='Dataset name.')
    parser.add_argument('--outfile_p', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles.txt', help='Output name for profile.')
    parser.add_argument('--outfile_s', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles_seqs.txt', help='Output name for seqeunces belonging to each profile.')
    parser.add_argument('--outfile_sub_p', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles_subYYY.txt', help='Output name for subsampled profile.')
    parser.add_argument('--outfile_sub_s', type=str, required=False, default='DATASET_NAMEtopXXX_aammp_profiles_seqs_subYYY.txt', help='Output name for sequences belonging to each subsampled profile.')
    parser.add_argument('--highest_Vgene', type=bool, required=False, default=False, help='Only report the highest frequency V gene.')
    parser.add_argument('--allele_finding', type=bool, required=False, default=False, help='If partis clustering was run using allele finding mode use this gene set instead.')
    parser.add_argument('--sub_mixing', type=bool, required=False, default=True, help='Prefer taking the subsamples from different filenames.')
    parser.add_argument('--collapse_clusters', type=str, required=False, default=None, help='Collapse clusters with the same naive sequence, either on amino acid level (AA) or DNA level (DNA) or None.')
    parser.add_argument('--minN_mixing', type=int, required=False, default=100, help='Minimum cluster size to be included in the filename mixing of subsamples.')
    parser.add_argument('--MIN_OBS', type=int, required=False, default=5, help='Minimum number of sequences in each clonal family.')
    parser.add_argument('--MIN_LEN', type=int, required=False, default=50, help='Minimum length of any given amino acid sequence in the dataset.')
    parser.add_argument('--MAX_REP', type=int, required=False, default=30, help='Max number of nt. repaired at each end.')
    parser.add_argument('--SIM_SIZE', type=int, required=False, default=10000, help='Number of random draws to simulate the neutral profile.')
    parser.add_argument('--LOCUS', type=str, required=False, default='igh', help='Locus, either igh, igk or igl.')
    parser.add_argument('--SPECIES', type=str, required=False, default='human', help='Species, either human.')

    global args
    args = parser.parse_args()
    assert((args.fsub, args.totN_sub).count(None) == 1)  # One parameter must be specified for the subsample

    # Read default germline info:
#    global glfo
#    glfo = glutils.read_glfo(partis_path + '/data/germlines/human', locus=args.LOCUS)

    # Insert dataset name in filename:
    outfile_p = args.outfile_p.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    outfile_s = args.outfile_s.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    outfile_sub_p = args.outfile_sub_p.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    outfile_sub_s = args.outfile_sub_s.replace('DATASET_NAME', '{}'.format(args.dataset_name))
    # Insert testset size in filename:
    outfile_p = outfile_p.replace('topXXX', '{}'.format(str(args.nsubs)))
    outfile_s = outfile_s.replace('topXXX', '{}'.format(str(args.nsubs)))
    outfile_sub_p = outfile_sub_p.replace('topXXX', '{}'.format(str(args.nsubs)))
    outfile_sub_s = outfile_sub_s.replace('topXXX', '{}'.format(str(args.nsubs)))
    # Insert subsample size or fraction in filename:
    if args.totN_sub is None:
        outfile_p = outfile_p.replace('YYY', '{}'.format(str(args.fsub)))
        outfile_s = outfile_s.replace('YYY', '{}'.format(str(args.fsub)))
        outfile_sub_p = outfile_sub_p.replace('YYY', '{}'.format(str(args.fsub)))
        outfile_sub_s = outfile_sub_s.replace('YYY', '{}'.format(str(args.fsub)))
    else:
        outfile_p = outfile_p.replace('YYY', 'N{}'.format(str(args.totN_sub)))
        outfile_s = outfile_s.replace('YYY', 'N{}'.format(str(args.totN_sub)))
        outfile_sub_p = outfile_sub_p.replace('YYY', 'N{}'.format(str(args.totN_sub)))
        outfile_sub_s = outfile_sub_s.replace('YYY', 'N{}'.format(str(args.totN_sub)))

    # Find all the files to process:
    glob_res = glob.glob(args.glob_pat)
    print 'These are the files that are going to be parsed:\n{}'.format('\n'.join(glob_res))
    # Check if there are saved results from previous run:
    saved_variables = '{}_saved_variables.pickle'.format(outfile_sub_p)
    if os.path.isfile(saved_variables):
        with open(saved_variables, 'rb') as f:
            profiles, all_numbered_trim, info_i_j, all_numbered_trim_sub, info_i_j_sub = pickle.load(f)
        print 'Loaded saved variables from pickled file:', saved_variables
    else:  # No pre-generated results so generate them now
        profiles, all_numbered_trim, info_i_j, all_numbered_trim_sub, info_i_j_sub = make_profile(glob_res, saved_variables)

    # Make subsampled profile:
    profiles_sub, info_i_j_sub = sub_par(all_numbered_trim_sub, info_i_j_sub, args)

    # Generate dataframes with the profiles and extra info (notice the slice from args.nsubs excluded):
    df = make_dataframe(profiles, info_i_j)
    df_sub = make_dataframe_sub(profiles_sub, info_i_j_sub)

    # Write the profiles:
    write_dataframe(df, outfile_p)
    write_dataframe(df_sub, outfile_sub_p)

    # Write the sequences that are the basis of the profiles:
    write_seqs(info_i_j, outfile_s)
    write_seqs_sub(info_i_j_sub, outfile_sub_s)


if __name__ == '__main__':
    # Make a tmp dir to dump crap:
    pretty_random_fnam = str(random.randint(1, 10**100))
    global TMPDIR
    TMPDIR = '/tmp/kd_tmp_' + pretty_random_fnam
    os.mkdir(TMPDIR)
    try:
        main()
        print 'Done'
    finally:
        shutil.rmtree(TMPDIR)  # rm -rf tmp dir


# ANARCI --sequence some_fasta.fa --outfile some_fasta.anno --scheme aho --restrict H --ncpu 1 --use_species human
# unique_ids,v_gene,d_gene,j_gene,cdr3_length,mut_freqs,input_seqs,indel_reversed_seqs,naive_seq,indelfos,duplicates,v_per_gene_support,d_per_gene_support,j_per_gene_support,v_3p_del,d_5p_del,d_3p_del,j_5p_del,v_5p_del,j_3p_del,vd_insertion,dj_insertion,fv_insertion,jf_insertion,mutated_invariants,in_frames,stops
# ((117, ' '), 'V')
# unique_ids,v_gene,d_gene,j_gene,cdr3_length,mut_freqs,input_seqs,indel_reversed_seqs,naive_seq,indelfos,duplicates,v_per_gene_support,d_per_gene_support,j_per_gene_support,v_3p_del,d_5p_del,d_3p_del,j_5p_del,v_5p_del,j_3p_del,vd_insertion,dj_insertion,fv_insertion,jf_insertion,mutated_invariants,in_frames,stops

# partis
# /fh/fast/matsen_e/kdavidse/partis/bin/partis run-viterbi --chain h --species human --infname inp.fasta --outfname out.csv
# cmd = '{}/bin/partis partition --chain {} --species {} --infname tmp/{}.fa --outfname tmp/{}.csv'.format(partis_path, chain, species, inpf, outf)
# out-cluster-annotations.csv
# out.csv
