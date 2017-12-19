#!/usr/bin/env python

import sys, os, glob, csv, random, copy, time, shutil, pickle
import gc
from itertools import cycle
csv.field_size_limit(sys.maxsize)
from anarci import anarci
from collections import Counter
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
this_file = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(1, this_file)
from neutral_profile import MutationModel
partis_path = this_file + '/partis/python'
sys.path.insert(1, partis_path)
import utils
import glutils

# Notice a gap is added as the 21th amino acid:
AA_LIST = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
AA_INDEX = {aa:i for i, aa in enumerate(AA_LIST)}
AHO_L = 149


