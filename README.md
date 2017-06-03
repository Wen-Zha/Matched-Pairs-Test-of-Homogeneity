# Matched-Pairs-Test-of-Homogeneity
Does 3.3.2 of Chptr 15 Lars S (IOME)

This python code runs the matched pairs test of symmetry (MPTS), matched pairs test of marginal symmetry (MPTMS) and matched pairst test of internal symmetry (MPTIS).

To run the software, feed a nexus alignment when prompted. It should return a csv file containing the data for each pairwise test (example file Cognato_2001data.csv), and a summary table (example file Cognato_2001tablebinom.csv). These files will be found in the root directory.

Example: to be added later

For sufficiently small datasets, plot(df,dset) will return a plot of the results for each gene (example Cognato_2001chart.png), and save this plot as an image to the root directory.

# Quick Start Guide for Jupyter
Change working directory to file path of Symmetry_Tests.py

    e.g. cd 'file path of Symmetry_Tests.py'
  
run program

    %run 'Symmetry_Tests.py'
  
when prompted input path of alignment

    input nex file here: file path of alignment/alignment.nex (no '' needed)
  
Your results (summary table and all results) will be outputed as a csv file in the working direcotry
