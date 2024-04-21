import argparse
from Bio.PDB import PDBParser, Superimposer, PDBIO
from pathlib import Path
from fastcluster import linkage
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform


def assign_chain_identifiers(structures):
   # Assign unique chain identifiers (A, B, C, ...) to each protein
   chain_id = 'A'
   for structure in structures:
       for chain in structure.get_chains():
           chain.id = chain_id
           chain_id = chr(ord(chain_id) + 1)


def count_common_contacts(structure1, structure2, cutoff=4.0):
   common_contacts = 0
   atoms1 = list(structure1.get_atoms())
   atoms2 = list(structure2.get_atoms())
   for atom1 in atoms1:
       for atom2 in atoms2:
           if atom1 - atom2 < cutoff:
               common_contacts += 1
               break  # Break the inner loop if a common contact is found
   return common_contacts


def superimpose_proteins(hit_files, output_file):
   structures = []


   # Parse structures
   for hit_file in hit_files:
       parser = PDBParser()
       structure = parser.get_structure(hit_file.stem, hit_file)
       structures.append(structure)


   # Assign unique chain identifiers
   assign_chain_identifiers(structures)


   # Initialize Superimposer
   superimposer = Superimposer()


   # Initialize output structure with the first structure
   output_structure = structures[0].copy()


   # Track chains already added to output structure
   added_chains = set()


   # Calculate fraction of common contacts and RMSD
   fcc_values = []
   rmsd_values = []
   for i, ref_structure in enumerate(structures):
       for j, target_structure in enumerate(structures):
           if i >= j:
               continue  # Skip self-comparison and duplicate comparisons


           # Calculate RMSD
           ref_atoms = list(ref_structure.get_atoms())
           target_atoms = list(target_structure.get_atoms())
           num_atoms = min(len(ref_atoms), len(target_atoms))
           ref_atoms = ref_atoms[:num_atoms]
           target_atoms = target_atoms[:num_atoms]
           superimposer.set_atoms(ref_atoms, target_atoms)
           superimposer.apply(target_structure)
           rmsd = superimposer.rms
           rmsd_values.append(rmsd)
           print(f"RMSD between {ref_structure.id} and {target_structure.id}: {rmsd:.2f}")


           # Calculate fraction of common contacts
           common_contacts = count_common_contacts(ref_structure, target_structure)
           total_atoms = len(list(ref_structure.get_atoms())) + len(list(target_structure.get_atoms()))
           fcc = common_contacts / total_atoms
           fcc_values.append(fcc)
           print(f"FCC between {ref_structure.id} and {target_structure.id}: {fcc:.2f}")


           # Add chains to output structure
           for chain in target_structure.get_chains():
               chain_id = chain.id
               if chain_id not in added_chains:
                   output_structure[0].add(chain.copy())
                   added_chains.add(chain_id)


   # Save the superimposed structures into a single PDB file
   io = PDBIO()
   io.set_structure(output_structure)
   io.save(str(output_file))  # Convert Path object to string before saving


   # Perform hierarchical clustering
   fcc_matrix = squareform(fcc_values)
   rmsd_matrix = squareform(rmsd_values)
   linkage_matrix = linkage(fcc_matrix)


   # Choose the number of clusters (k) based on your requirements
   k = 5
   clusters = fcluster(linkage_matrix, k, criterion='maxclust')
   print("Cluster assignments:")
   for i, cluster_id in enumerate(clusters):
       print(f"{structures[i].id}: Cluster {cluster_id}")


if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Superimpose protein structures.")
   parser.add_argument("files", metavar="FILE", type=str, nargs="+",
                       help="Filenames of protein structures to superimpose")
   parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", type=str,
                       help="Filename to save the superimposed structure", required=True)
   args = parser.parse_args()


   # Convert filenames to Path objects
   hit_files = [Path(file) for file in args.files]
   output_file = Path(args.output)


   # Call the superimpose_proteins function
   superimpose_proteins(hit_files, output_file)


   print(f"Superimposed structure saved to: {output_file}")


