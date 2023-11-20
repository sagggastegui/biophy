from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
import numpy as np

MAXDIST =   # Define distance for a contact

parser = PDBParser(QUIET=True)

# Load structure from PDB file
st = parser.get_structure('6M0J', '6m0j_fixed.pdb')

# Select CA atoms
select = [atom for model in st for chain in model for residue in chain for atom in residue if atom.id == 'CA']

select_ats = []
for at in st.get_atoms():
   # Skip Hydrogens to reduce time
    if at.element != 'H':
        select_ats.append(at)

# Preparing search
nbsearch = NeighborSearch(select_ats)

print("Contacts between residues in different chains:")

# Searching for contacts under MAXDIST
for at1, at2 in nbsearch.search_all(MAXDIST):
    res1 = at1.get_parent()
    res2 = at2.get_parent()
    chain1 = res1.get_parent().id
    chain2 = res2.get_parent().id
    resname1 = res1.get_resname()
    resname2 = res2.get_resname()

    if chain1 != chain2:  # Ensure residues are from different chains and avoid duplicates
        distance = np.linalg.norm(at1.coord - at2.coord)
        print(resname1, res1.id[1], chain1,'-', resname2, res2.id[1], chain2, '-', distance)

