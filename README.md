# DiagnoSeq
Detect diagnostic alleles between hybrid species and its parental species with a DNA-based data.

The script "daignoSeq.py" is given paternal and maternal fasta files (in this order) including all data of all individuals, and it takes the fast files of the hybrid individuals with ".fa" extension from the current working directory one-by-one as inputs, and it gives two outputs for each hybrid individual with names "<NameOfTheHybridIndv>_P" for diagnostic alleles from paternal side, and "<NameOfTheHybridIndv>_M" for the diagnostic alleles from the maternal side. These files are made up of two columns, first column includes number of locus, second column includes the sequence of the diagnostic allele.

