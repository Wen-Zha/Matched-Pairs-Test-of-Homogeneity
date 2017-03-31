import numpy as np
import itertools as ite
from scipy.stats import chi2

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
            print "error"
    return m
def MPTS(m):
    s = 0.0
    for (i,j) in ite.product(range(0,4),range(0,4)):
        if i<j:
            n = (m[i,j]-m[j,i])**2
            d = m[i,j]+m[j,i]
            s = s+(float(n)/float(d))
    return s

if __name__ == '__main__': 
    #some model data
    s1 = "AAACTGTTCTACAGGACTGATG"
    s2 = "CATCTGTAATAGGCCACTGATA"
    
    seq1= "TTTCTGTAGACTACAGCCGAACTGATACAATACAAGCACAAACAATTCACCGCGTCGCGCACAGTCGTCAAAGCGGCATTCCATAAAAGTTCATCCATACCCCGAGGTAACCTCACGTCGTCACGGGCTGACGTAATCACGAAAGCACCGCCCGACCGGTCAAGCCTCAGAAGGGTCGAACACGGACTCAGTCTCAAGTGCTCCTCCACAAACGTCATACTTAGTTCACCATCCCCGAGCCTATTTCCCTTAAAATGCGGTAACCCGGCCAGGGAGGAGAGAAAGAGTGG"
    seq2= "CGTCTGGGATCTTTTGCCGGGCTGGGTCGCTACACGAACGCAGAGTTCTACTCCGGTCGCACTTGCGGATGAGTTGGTTACGGAGAGTGCGGGTCTTTTCCCAAAGTTCATTTCCCGTCGTTTCGGCCTGTTGTAATCATGTGTGCTCCGCCCCATCGGTGAAGCCCCGCTAGCGTATTACTCGGAATGTGTATCTAGTGCCAATTCATATACGTATAAGTTAGTTTATAATCTCCTCGCCTATTTCCTTAGAAATAGTATTATCGATCTTTGACGGAGTGAACTATGGG"
    seq3= "CTACAGTTAAGTTCTGCAGAGCTGCTTGACTATACGATCAACGAATACAAGACGGGGCGCACAGGCATATAAGTGGGATTCCGTAAGATCATGTCTCTACCCAAAGGGTACATGTTGTCTTCACGGCCAAACCTAATCACGGGTCACCTGCCCAACAGTTGAAGGCGCGCCAGGCCGGCCCACGCATACAGACTCCAGAGCAACTCCATCAACGTCAAACTTACCTCAAAGTCTCCGCGGCTAGTCCCATTGAAATACGATTATCTCACCTTGCAAGAGTGAAAAAATCG"
    seq4= "CTTCGGTATAGTTCTGCCGAGCTGGTTCGCTACATGATCAATGATTACGACCCTGGGCCCTCTGGCGTATGAGTGGGATGGTGTCAAATTTCTTCTTGACCGGCAGGTCACCTCTTGTCCTGAGGGCCGGGCGGCAGCAGGTCTGTGCTGTTTTGCCTGTAATGCCTCGTCAGGCGGGAGCACGGTTTTAGTATCCTTGCCTACTCTATTATTCTGAAATTTAGCTGATAATCTCTTCAGCTAATTCTTTAGAAATAGGCTTATCGTCCCGGGTTGGTGCGAAACATCCG"
    #make the matrix of proportions
    m = matrix(seq2,seq3)
    #generate a symmetry score
    s = MPTS(m) 
    
    print(1 - chi2.cdf(s,6.0))