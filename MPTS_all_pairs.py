import numpy as np
import itertools as ite
from scipy.stats import chi2
from Bio.Nexus import Nexus
from Bio import AlignIO

#import scipy.stats as sp

def matrix(s1,s2):
    m= np.zeros((4,4))
    elog = [] #error log
    for i in range(len(s1)):
        #empty cases
        if s1[i] == "?":
            m = m
        elif s2[i] == "?":
            m = m
        elif s1[i] == "-":
            m = m
        elif s2[i] == "-":
            m = m
        #cases of A->A
        elif s1[i]==s2[i]:
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
            elog.append([s1[i],s2[i],"error in matrix"])
    return m, elog
def MPTS(m):
    """ inputs
            m: a 4x4 matrix of proportions
        outputs
            p: is a p-value for the matched pairs test of symmetry
    """
    s = 0.0
    for (i,j) in ite.product(range(0,4),range(0,4)):
        if i<j:
            n = (m[i,j]-m[j,i])**2
            d = m[i,j]+m[j,i]
            if float(d) != 0.:
                s = s+(float(n)/float(d)) 
                p = 1 - chi2.cdf(s,6.0)
            else:
                p = 'NA' #note p gets overridden later so better to find a fix here.

    return p
    
def MPTS_aln(aln):
    """inputs 
            charset_aln = alignment array of sites
        output
            p = array of pvalues with corresponding charset pair, (0,1),(0,2),...,(0,43),(1,2)...
    
    """
    #use itertools to loop over all sequences (itertools all pairs)
    # make iterator for all pairs of rows in the numpy array called charset_aln

    # get list of p values where each entry is a p-value for a single pair of rows

    # return that list
    aln_array = np.array([list(rec) for rec in aln], np.character)

    dat.charsets.keys() #these are the names to the CHARSETS in the .nex file, which you can iterate over in a for loop
    for name in dat.charsets:
        sites = dat.charsets[name]
        # slice the alignment to get the sub-alignment for the CHARSET
        charset_aln = aln_array[:, sites]
    p = [] #empty list of p values
    for n in dat.charsets.keys():
        for q in ite.combinations(list(range(len(aln))),2): #iterating over all taxa for sites
            m, elog = matrix(aln_array[:,dat.charsets[n]][q[0]].tostring().upper().decode(),aln_array[:,dat.charsets[n]][q[1]].tostring().upper().decode())
        
            p.append([n,q,MPTS(m)]) 
    
    return p, elog  

if __name__ == '__main__': 
    aln_path = '/Users/user/Documents/! ! 2017 ANU/Semester 1/SCNC2103/data reader/alignment.nex'
    
    dat = Nexus.Nexus()

    dat.read(aln_path)

    aln = AlignIO.read(open(aln_path), "nexus")
    p, elog = MPTS_aln(aln)
    print(p)
    print(len(p))