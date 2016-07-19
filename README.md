# Readme file for cloning primer designer script

This program designs your cloning primers for traditional cloning with the desired meltig temperature and without unwanted secondary structure.

* First you provide the gene name and sufficient length of homology regions of the beginning and end of your gene of interest. Recommended length is 30-40 bp. Both sequences are provided in forward direction.
* Then the homology regions get selected by finding the shortest primers that have a Tm of 60C or more using NEB's Q5 polymerase Tm calculator. This will be the homology Tm.
* Next the restriction endonuclease sites and any other additional sequences (provided in forward orientation) are added onto the primers.
* Then a 6 bp ATATTA sequence is added on the 5' end of the primers, and its GC content is optimized such that the two primer Tms match up.
* Finally, the 6 bp flanking regions are permuted and the ones with the least secondary structure -smallest absolute deltaG values of homo- and heterodimers using IDT's oligocalc algorithm- are selected.
* The output is a text file containing the sequence of your primers, the homology region Tms and the common Tm to be used in PCR, deltaG values and the full primer Tms and the common Tm to be used in PCR.

It is recommended that you use these primers with a two-step PCR where the first 5 cycles have an annealing temperature at the common homology Tm, and the next 30 cycles' annealing temperature is the common Tm for the full length primers.

Requirements:

Chrome with chromedriver

The program opens up two Chrome windows to perform the necessary calculations on NEB's and IDT's websites.

