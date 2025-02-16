from Bio import Phylo
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load your tree from a Newick file (adjust the filename accordingly)
tree = Phylo.read("your_tree_SARS.newick", "newick")

# Group tips by variant.
# This assumes tip names are in the format: "Variant_identifier" (e.g. "Alpha_001", "Beta_002", etc.)
variants = {}
for tip in tree.get_terminals():
    variant = tip.name.split("_")[0]
    variants.setdefault(variant, []).append(tip.name)

# Compute the MRCA for each variant.
variant_mrca = {}
for variant, tip_names in variants.items():
    variant_mrca[variant] = tree.common_ancestor(tip_names)
    print(f"MRCA for {variant}: {variant_mrca[variant].name if variant_mrca[variant].name else 'Unnamed Node'}")

# Create a matrix (DataFrame) for distances between variant MRCAs.
variants_list = list(variant_mrca.keys())
distance_matrix = pd.DataFrame(index=variants_list, columns=variants_list, dtype=float)

for var1 in variants_list:
    for var2 in variants_list:
        # Compute distance between the MRCA of var1 and the MRCA of var2.
        distance_matrix.loc[var1, var2] = tree.distance(variant_mrca[var1], variant_mrca[var2])

# Save the matrix as a text file (tab-separated).
distance_matrix.to_csv("variant_mrca_distance_matrix.txt", sep="\t")
print("Distance matrix saved to variant_mrca_distance_matrix.txt")

# Generate a heatmap of the distance matrix.
plt.figure(figsize=(10, 8))
sns.heatmap(distance_matrix.astype(float), annot=True, fmt=".2f", cmap="viridis")
plt.title("Distance Between MRCA of Each Variant")
plt.xlabel("Variant")
plt.ylabel("Variant")
plt.tight_layout()
plt.savefig("variant_mrca_distance_heatmap.png")
plt.show()