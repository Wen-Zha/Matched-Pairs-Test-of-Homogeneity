import numpy as np
import itertools as ite
from scipy.stats import chi2
from Bio.Nexus import Nexus
from Bio import AlignIO

#import scipy.stats as sp

def matrix(s1,s2):
    m= np.zeros((4,4))
    for i in range(len(s1)):
        #cases of A->A
        if s1[i]==s2[i]:
            if s1[i]=="A":
                m[0,0]=m[0,0]+1
            elif s1[i]=="C":
                m[1,1]=m[1,1]+1
            elif s1[i]=="G":
                m[2,2]=m[2,2]+1
            elif s1[i]=="T":
                m[3,3]=m[3,3]+1
        #other cases
        elif s1[i]=="A":
            if s2[i]=="C":
                m[1,0]=m[1,0]+1
            elif s2[i]=="G":
                m[2,0]=m[2,0]+1
            elif s2[i]=="T":
                m[3,0]=m[3,0]+1
        elif s1[i]=="C":
            if s2[i]=="A":
                m[0,1]=m[0,1]+1
            elif s2[i]=="G":
                m[2,1]=m[2,1]+1
            elif s2[i]=="T":
                m[3,1]=m[3,1]+1
        elif s1[i]=="G":
            if s2[i]=="A":
                m[0,2]=m[0,2]+1
            elif s2[i]=="C":
                m[1,2]=m[1,2]+1
            elif s2[i]=="T":
                m[3,2]=m[3,2]+1
        elif s1[i]=="T":
            if s2[i]=="C":
                m[1,3]=m[1,3]+1
            elif s2[i]=="G":
                m[2,3]=m[2,3]+1
            elif s2[i]=="A":
                m[0,3]=m[0,3]+1           
        else:
            print("error")
    return m
def MPTS(m):
    s = 0.0
    for (i,j) in ite.product(range(0,4),range(0,4)):
        if i<j:
            n = (m[i,j]-m[j,i])**2
            d = m[i,j]+m[j,i]
            s = s+(float(n)/float(d))
    return s

def Alignment(aln_path):

    dat = Nexus.Nexus()

    dat.read(aln_path)

    aln = AlignIO.read(open(aln_path), "nexus")

    aln_array = np.array([list(rec) for rec in aln], np.character)

    dat.charsets.keys() #these are the names to the CHARSETS in the .nex file, which you can iterate over in a for loop

    results = {} # it would be better to use a Pandas data frame
    for name in dat.charsets:
        sites = dat.charsets[name]
        # slice the alignment to get the sub-alignment for the CHARSET
 
        charset_aln = aln_array[:, sites]

        results[name] = MPTS_aln(charset_aln)
    
    return results[name]
    
def MPTS_aln(charset_aln):
    # make iterator for all pairs of rows in the numpy array called charset_aln

    # get list of p values where each entry is a p-value for a single pair of rows

    # return that list
    
    pvalues={}
    pvalues.append
    
    return pvalues

if __name__ == '__main__': 
    print("hello")
    s1 = "AAACTGTTCTACAGGACTGATG"
    s2 = "CATCTGTAATAGGCCACTGATA"
    
    seq1= "TTTCTGTAGACTACAGCCGAACTGATACAATACAAGCACAAACAATTCACCGCGTCGCGCACAGTCGTCAAAGCGGCATTCCATAAAAGTTCATCCATACCCCGAGGTAACCTCACGTCGTCACGGGCTGACGTAATCACGAAAGCACCGCCCGACCGGTCAAGCCTCAGAAGGGTCGAACACGGACTCAGTCTCAAGTGCTCCTCCACAAACGTCATACTTAGTTCACCATCCCCGAGCCTATTTCCCTTAAAATGCGGTAACCCGGCCAGGGAGGAGAGAAAGAGTGG"
    seq2= "CGTCTGGGATCTTTTGCCGGGCTGGGTCGCTACACGAACGCAGAGTTCTACTCCGGTCGCACTTGCGGATGAGTTGGTTACGGAGAGTGCGGGTCTTTTCCCAAAGTTCATTTCCCGTCGTTTCGGCCTGTTGTAATCATGTGTGCTCCGCCCCATCGGTGAAGCCCCGCTAGCGTATTACTCGGAATGTGTATCTAGTGCCAATTCATATACGTATAAGTTAGTTTATAATCTCCTCGCCTATTTCCTTAGAAATAGTATTATCGATCTTTGACGGAGTGAACTATGGG"
    seq3= "CTACAGTTAAGTTCTGCAGAGCTGCTTGACTATACGATCAACGAATACAAGACGGGGCGCACAGGCATATAAGTGGGATTCCGTAAGATCATGTCTCTACCCAAAGGGTACATGTTGTCTTCACGGCCAAACCTAATCACGGGTCACCTGCCCAACAGTTGAAGGCGCGCCAGGCCGGCCCACGCATACAGACTCCAGAGCAACTCCATCAACGTCAAACTTACCTCAAAGTCTCCGCGGCTAGTCCCATTGAAATACGATTATCTCACCTTGCAAGAGTGAAAAAATCG"
    seq4= "CTTCGGTATAGTTCTGCCGAGCTGGTTCGCTACATGATCAATGATTACGACCCTGGGCCCTCTGGCGTATGAGTGGGATGGTGTCAAATTTCTTCTTGACCGGCAGGTCACCTCTTGTCCTGAGGGCCGGGCGGCAGCAGGTCTGTGCTGTTTTGCCTGTAATGCCTCGTCAGGCGGGAGCACGGTTTTAGTATCCTTGCCTACTCTATTATTCTGAAATTTAGCTGATAATCTCTTCAGCTAATTCTTTAGAAATAGGCTTATCGTCCCGGGTTGGTGCGAAACATCCG"
    
    m1 = matrix(seq2,seq3)
    m = np.array([[191, 71, 68, 57],[14, 142, 22, 33],[16, 12, 144, 29],[26, 19, 18, 138]]) 
    MA = np.array([[248, 1, 1, 1],[1, 251, 1, 1],[1, 1, 252, 1],[1, 1, 1, 249]])
    #m=np.array([[248,0,0,0],[0,251,0,0],[0,0,252,0],[0,0,0,249]])
    s = MPTS(m1) 
    print("HELLO")
    
    print(1 - chi2.cdf(s,6.0))
    
    aln = Alignment('/Users/user/Documents/! ! 2017 ANU/Semester 1/SCNC2103/data reader/alignment.nex')
    #use itertools to loop over all sequences (itertools all pairs)