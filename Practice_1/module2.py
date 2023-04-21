import matplotlib.pyplot as plt
import numpy as np
from Levenshtein import distance
from Bio import Align, pairwise2, SeqIO
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from Bio.SeqUtils import seq3, gc_fraction, MeltingTemp as melting_temp, GC123, GC_skew, xGC_skew, nt_search
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, MMCIFParser
import py3Dmol
from collections import Counter
import os

here = os.path.dirname(os.path.abspath(__file__))

dna_seq = Seq('ATGATCTCGTAA')

record = SeqRecord(dna_seq, id="test", annotations={"molecule_type": "DNA"})

# print(Counter(seq))

dna_frequency = Counter(dna_seq)

# plt.bar(dna_frequency.keys(), dna_frequency.values())


# Transcribe dna to mRna and then to an Amino Acid
dna_seq.complement()

m_rna = dna_seq.transcribe()

amino_acid = m_rna.translate()

# Direct Translation:
aa = dna_seq.translate()

# transcribe back to DNA:
m_rna.back_transcribe()

# Get full Amino Acid name
seq3(aa)

# Dna Codontable
CodonTable.unambiguous_dna_by_name['Standard']

# Rna CodonTable
CodonTable.unambiguous_rna_by_name['Standard']


# GC Contents in DNA
# higher GC content indicates a relatively higher melting temperature
gc_fraction(dna_seq)


melting_temp.Tm_Wallace(dna_seq)
melting_temp.Tm_GC(dna_seq)

def highest_gc(seq1, seq2):
    gc1 = gc_fraction(seq1) * 100
    gc2 = gc_fraction(seq2) * 100
    at1 = 100 - gc1
    at2 = 100 - gc2
    melting1 = melting_temp.Tm_GC(seq1)
    melting2 = melting_temp.Tm_GC(seq2)
    result1 = "Sequence: {} , GC: {} , AT: {} , Temp: {}".format(seq1, gc1, at1, melting1)
    result2 = "Sequence: {} , GC: {} , AT: {} , Temp: {}".format(seq2, gc2, at2, melting2)

    if gc1 > gc2:
        return result1
    else:
        return result2
    
# GC Skew
# print(GC123(dna_seq))

# print(GC_skew(dna_seq, 10))

# Sequence Alighnment
seq1 = Seq('ACTCGT')
seq2 = Seq('ATTCG')

# Global Alighnments
alignments = pairwise2.align.globalxx(seq1,seq2)
# for a in alignments:
#     print(format_alignment(*a))

# local Alighnments
loc_alignments = pairwise2.align.localxx(seq1,seq2)
# for a in loc_alignments:
#     print(format_alignment(*a))


test_alignment = pairwise2.align.globalms(seq1,seq2, 2, -1, -0.5, -0.1)
# for a in test_alignment:
    #  print(format_alignment(*a))

seqA = Seq('AAGGCTT')
seqB = Seq('AAGGC')
seqC = Seq('AAGGCAT')

alignerLocal = Align.PairwiseAligner()
alignerLocal.mode = "local"

AvB = alignerLocal.align(seqA, seqB)
BvC = alignerLocal.align(seqB, seqC)
AvC = alignerLocal.align(seqA, seqC)
# print('AvB: ' ,AvB.score/len(seqB)*100)
# print('BvC: ' ,BvC.score/len(seqB)*100)
# print('AvC: ' ,AvC.score/len(seqC)*100)

# Hamming Distance fxn

seq1 = Seq('ACTAT')
seq2 = Seq('ACTTA')
seq3 = Seq('ACTT')
def hamming_distance(lhs,rhs):
    return( len([(x,y) for x,y in zip(lhs,rhs) if x!=y]))
hamming_distance(seq1, seq2)

# 0 if match
hamming_distance(seq1, seq1)
hamming_distance(seq1, seq3)

# Levenshtein Distance
distance(str(seq1), str(seq2))

# print("Hamming Distance", hamming_distance(seq1, seq3))
# print("Levenshtein Distance", distance(str(seq1), str(seq3)))


# stack overflow code here https://stackoverflow.com/questions/40822400/how-to-create-a-dotplot-of-two-dna-sequence-in-python
def delta(x,y):
    return 0 if x == y else 1

def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))

def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M,t, seq1, seq2, nonblank = chr(10000), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)

def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice

# Sequence Alignment with Dot Plot
seq1 = Seq('ACTTAG')
seq2 = Seq('AC')


def dotplotx(seq1,seq2,k = 1,t = 1):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    plt.xticks(np.arange(len(list(seq1))),list(seq1))
    plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()

dna1 = Seq('ATGATCTCGTAA')
dna2 = Seq('ATTATGTCGTAA')

# dotplotx(dna1, dna2)

# Query Data with file formats 
# for record in SeqIO.parse("sequence.fasta", "fasta"):
#     print(record)

# Reading FASTA
filename = os.path.join(here, 'sequence.fasta')
dna_record = SeqIO.read(filename, "fasta")
dna_seq = dna_record.seq

# Reading Genbank
filename = os.path.join(here, 'sequence.gb')
gb_dna_record = SeqIO.parse(filename, "genbank")
    

# Working with 3D structures
parser = PDBParser()
filename = os.path.join(here, '6lu7.pdb')
structure = parser.get_structure('6LU7', filename)
model = structure[0]

# structure => model => chain => residue => atom
# check for chains
# for chain in model:
#     for atom in chain:
#         print(atom)

# Using py3DMol
view1 = py3Dmol.view(query='pdb:6LU7')
view1.setStyle({'cartoon': {'color':'spectrum'}})
view1.addModel()
view1.show()

# Create a Py3Dmol view object
view = py3Dmol.view()

# Define a small molecule in SMILES format
smiles = 'C[C@@H](O)[C@H](N)C(=O)O'

# Add the molecule to the viewer
view.addModel(smiles, 'smiles')

# Set the style of the molecule
view.setStyle({'stick':{}})

# Display the viewer
view.show()
# not displaying

