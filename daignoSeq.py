#DIAGNOSTIC ALLELE DETECTION SCRIPT COMPREHENSION

# Get the libraries
import pandas as pd
import sys


# Take the parental objects (FASTA files)
with open(sys.argv[1], "r") as p1:
    with open(sys.argv[2], "r") as m1:

        # Load the data
        # Create Pandas object from input fasta files
        p0 = pd.read_csv(p1, sep="\t", header=None)
        m0 = pd.read_csv(m1, sep="\t", header=None)


        # Convert the single-columns fasta format to two-column data frame
        rowIndexp0 = p0.index[:]
        idxp0 = rowIndexp0.values.tolist()

        even_idxp0 = [x for x in idxp0 if x%2 == 0]
        odd_idxp0 = [x for x in idxp0 if x%2 != 0]

        p0DF = pd.DataFrame(columns=['SeqNames', 'Seq'])

        p0DF['SeqNames'] = list(p0.loc[even_idxp0, 0])
        p0DF['Seq'] = list(p0.loc[odd_idxp0, 0])


        rowIndexm0 = m0.index[:]
        idxm0 = rowIndexm0.values.tolist()

        even_idxm0 = [x for x in idxm0 if x%2 == 0]
        odd_idxm0 = [x for x in idxm0 if x%2 != 0]

        m0DF = pd.DataFrame(columns=['SeqNames', 'Seq'])

        m0DF['SeqNames'] = list(m0.loc[even_idxm0, 0])
        m0DF['Seq'] = list(m0.loc[odd_idxm0, 0])


        # Change the seq names as Locus Numbers in 1st column of
        # the data frame of parentals
        for seqName in range(0, len(m0DF)):
            m0DF.at[seqName, 'SeqNames'] = m0DF.iloc[seqName, 0].split('_')[5]

        for seqName in range(0, len(p0DF)):
            p0DF.at[seqName, 'SeqNames'] = p0DF.iloc[seqName, 0].split('_')[5]


        # Export the maternal file
        m0DF.to_csv('m0DF_2col', sep = '\t', index=False, header=None)

        # Export the paternal file
        p0DF.to_csv('p0DF_2col', sep = '\t', index=False, header=None)


## Hybrid code part
# Get hybrid files from the current working directory as "*.fa" extension
import glob
h0_list = glob.glob("*.fa")
for h_file in h0_list:
    h0 = pd.read_csv(h_file, sep="\t", header=None)

    # Convert the single-columns fasta format to two-column data frame
    rowIndexh0 = h0.index[:]
    idxh0 = rowIndexh0.values.tolist()

    even_idxh0 = [x for x in idxh0 if x%2 == 0]
    odd_idxh0 = [x for x in idxh0 if x%2 != 0]

    h0DF = pd.DataFrame(columns=['SeqNames', 'Seq'])

    h0DF['SeqNames'] = list(h0.loc[even_idxh0, 0])
    h0DF['Seq'] = list(h0.loc[odd_idxh0, 0])

    # Change the SeqNames in the 1st column as a locus number in the hybrid data
    for seqName in range(0, len(h0DF)):
        h0DF.at[seqName, 'SeqNames'] = h0DF.iloc[seqName, 0].split('_')[5]

    # Export the hybrid indv. data frame
    h0DF.to_csv('h0DF_2col', sep = '\t', index=False, header=None)

    
    # Extract number of locus of the hybrid sample
    locusNumUniqH = h0DF['SeqNames'].unique().tolist()

    # Create paternal and maternal output files for diagnostic alleles
    PatDF = pd.DataFrame(columns=['SeqNames', 'Seq'])
    MatDF = pd.DataFrame(columns=['SeqNames', 'Seq'])


    for line in locusNumUniqH:
        # Get index numbers of locus "line" in hybrid, pataernal and maternal data
        hidx = h0DF[h0DF.SeqNames == line].SeqNames.index.tolist()
        pidx = p0DF[p0DF.SeqNames == line].SeqNames.index.tolist()
        midx = m0DF[m0DF.SeqNames == line].SeqNames.index.tolist()

        if len(pidx) > 0 and len(midx) > 0 and len(hidx) > 0:
            # Check uniqueness of parental seqs
            if ~p0DF.loc[pidx, 'Seq'].isin(m0DF.loc[midx, 'Seq']).any().any() and ~m0DF.loc[midx, 'Seq'].isin(p0DF.loc[pidx, 'Seq']).any().any():
                # Inspect if hybrid 1st allle is in paternal side
                inspectorPh0 = pd.DataFrame([h0DF.loc[hidx[0],'Seq']]).isin(list(p0DF.loc[pidx, 'Seq'])).any().any()
                # Inspect if hybrid 1st allle is in maternal side
                inspectorMh0 = pd.DataFrame([h0DF.loc[hidx[0],'Seq']]).isin(list(m0DF.loc[midx, 'Seq'])).any().any()
                # Inspect if hybrid 2nd allle is in paternal side
                inspectorPh1 = pd.DataFrame([h0DF.loc[hidx[1],'Seq']]).isin(list(p0DF.loc[pidx, 'Seq'])).any().any()
                # Inspect if hybrid 2nd allle is in maternal side
                inspectorMh1 = pd.DataFrame([h0DF.loc[hidx[1],'Seq']]).isin(list(m0DF.loc[midx, 'Seq'])).any().any()
                # Separate the hybrid alleles as diagnostic paternal and maternal
                # and save them in the output files
                if inspectorPh0 and ~inspectorMh0 and ~inspectorPh1 and inspectorMh1:
                    PatDF = PatDF.append(pd.DataFrame(h0DF.loc[hidx[0],:]).T, ignore_index=True)
                    MatDF = MatDF.append(pd.DataFrame(h0DF.loc[hidx[1],:]).T, ignore_index=True)
                elif ~inspectorPh0 and inspectorMh0 and inspectorPh1 and ~inspectorMh1:
                    PatDF = PatDF.append(pd.DataFrame(h0DF.loc[hidx[1],:]).T, ignore_index=True)
                    MatDF = MatDF.append(pd.DataFrame(h0DF.loc[hidx[0],:]).T, ignore_index=True)

    # Export the paternal and maternal diagnostic allele files of the hybrid indv.
    PatDF.to_csv("./{0}_P".format(h_file), sep = '\t', index=False, header=None)
    MatDF.to_csv("./{0}_M".format(h_file), sep = '\t', index=False, header=None)
