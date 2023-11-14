from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import pandas as pd
# Read SNP data from CSV
snp_data = pd.read_csv('snp_data.csv', index_col=0)  # Assume the first column contains individual names

# Convert SNP data to a format suitable for phylogenetic analysis
sequences = []
for index, row in snp_data.iterrows():
    sequence = ''.join(row)
    record = SeqRecord(Seq(sequence), id=index)
    sequences.append(record)

# Create a MultipleSeqAlignment object
alignment = MultipleSeqAlignment(sequences)

# Calculate distance matrix using 'blosum62' model
calculator = DistanceCalculator('blosum62')
dm = calculator.get_distance(alignment)

# Construct the phylogenetic tree
constructor = DistanceTreeConstructor(calculator)
tree = constructor.upgma(dm)

# Save the tree to a Newick file
Phylo.write(tree, 'phylogenetic_tree.nwk', 'newick')

#Visualise .nwk using https://itol.embl.de/
