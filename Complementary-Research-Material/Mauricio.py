#!/usr/bin/env python
# coding: utf-8

# In[24]:


#Adapter from here: http://www.cbs.dtu.dk/courses/27624/Exercise_1.php
from Bio import SeqIO

#Reading a genbank file
#record = SeqIO.read("C:/Users/Mauricio/Desktop/7BUF_A.fasta" , "fasta" )
#print (record.id)
#print (record.description)
#print (record.name)
#print (record.seq)


# In[27]:


#reading from FASTA file
record = SeqIO.read ("C:/Users/Mauricio/Desktop/files/7BUF_A.fasta" , "fasta"  )
print (record.description)
print (record.seq)


# In[30]:


#Multiple alignment sequence
from Bio import SeqIO
import os
#Converting multiple SwissProt files to a single Fasta file
records = []
for filename in os.listdir(r"C:\Users\Mauricio\Desktop\files"):
    handle = open(r"C:/Users/Mauricio/Desktop/files" + "/" + filename)
    record = SeqIO.read( handle, "fasta" )
    records.append ( record )

print(len(records))


# In[31]:


SeqIO.write(records, r"C:/Users/Mauricio/Desktop/files_variants.fasta", "fasta" )


# In[33]:


#Performing a Multiple sequence alignment using ClustalOmega
from Bio.Align.Applications import ClustalOmegaCommandline
in_file = "C:/Users/Mauricio/Desktop/files_variants.fasta"
out_file= "LHBs_variants.aln"
newtree_file="LHBs_variants.dnd"
clustalo_exe = r"C:\Users\Mauricio\Desktop\clustal-omega-1.2.2-win64\clustalo.exe"
clustalo_cline = ClustalOmegaCommandline(clustalo_exe, infile=in_file,outfile= out_file, verbose= True, auto=True, outfmt="clustal",guidetree_out=newtree_file,force=True)
assert os.path.isfile(clustalo_exe), "Clustal O executable missing"
stdout, stderr = clustalo_cline()
print(clustalo_cline)


# In[34]:


#import the AlignIO module from BioPython
from Bio import AlignIO
alignment = AlignIO.read( open(out_file) , "clustal" )
print("Alignment length %i" % alignment.get_alignment_length())
print(alignment[0])


# In[35]:


#Calculating summary information
from Bio.Align import AlignInfo
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()
print(consensus)


# In[36]:


#Position Specific Score Matrices
my_pssm = summary_align.pos_specific_score_matrix( consensus , chars_to_ignore=['N'])
print(my_pssm[0]["M"])
print(my_pssm)


# In[14]:


#Since SubsMat might be removed in the future. Here is another way. Not explained in video.
from Bio.Align import substitution_matrices
names = substitution_matrices.load()

observed_frequencies = alignment.substitutions
print(observed_frequencies)
import numpy as np
observed_frequencies /= np.sum(observed_frequencies)
residue_frequencies = np.sum(observed_frequencies, 0)
print(format(residue_frequencies, "%.4f"))

expected_frequencies = np.dot(residue_frequencies[:, None], residue_frequencies[None, :])
print(format(expected_frequencies, "%.4f"))

m = np.log2(observed_frequencies/expected_frequencies)
print(m)
#a warning will show since log2 operations of 0 are not possible, this will be placed as inf or ND.


# In[39]:


#Construct a Phylogenetic Tree using Bio.Phylo
from Bio import Phylo
tree = Phylo.read("LHBs_variants.dnd", "newick")
Phylo.draw_ascii(tree)


# In[ ]:




