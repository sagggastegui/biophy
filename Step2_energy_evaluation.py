import argparse
import sys
import os
import math

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBIO import PDBIO, Select
#This are functions that you will need to import the parameters for VanderWaals or the residue library:
class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.residue_data[r.id] = r
        self.nres = len(self.residue_data)

    def get_params(self, resid, atid):
        atom_id = resid + ':' + atid
        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        else:
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    def __init__(self,data):
        self.id     = data[0]+':'+data[1]
        self.at_type = data[2]
        self.charge  = float(data[3])
        
class AtomType():
    def __init__(self, data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): #extracted from GELPI's github
    #parameters for the VdW
    def __init__ (self, file_name):
        self.at_types = {}
        try:
            fh = open(file_name, "r")
        except OSError:
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            self.at_types[data[0]] = AtomType(data)
        self.ntypes = len(self.at_types)
        fh.close()
# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib('/home/juliapm/biofi/aaLib.lib')
# loading VdW parameters
ff_params = VdwParamset('/home/juliapm/biofi/vdwprm')
# set the pdb_path and load the structure
pdb_path = "/home/juliapm/biofi/6m0j_fixed.pdb"
# Setting the Bio.PDB.Parser object
parser = PDBParser(PERMISSIVE=1)
# Loading structure
st = parser.get_structure('st', pdb_path)
#Possible Atom names that correspond to Ala atoms"
ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}
def residue_id(res):
    '''Returns readable residue id'''
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def atom_id(at):
    '''Returns readable atom id'''
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)
  def MH_diel(r):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    '''Electrostatic interaction energy between two atoms at r distance'''
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)

def calc_solvation(st, res):
    '''Solvation energy based on ASA'''
    solv = 0.
    solv_ala = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
        if at.id in ala_atoms:
            solv_ala += s
    return solv, solv_ala
def add_atom_parameters(st, residue_library, ff_params):
    ''' Adds parameters from libraries to atom .xtra field
        For not recognized atoms, issues a warning and put default parameters
    '''
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = residue_library.get_params(resname, at.id)
        if not params:
            print("WARNING: residue/atom pair not in library ("+atom_id(at) + ')')
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
add_atom_parameters(st, residue_library, ff_params)
def get_interface(st, dist):
    ''' Detects interface residues within a distance(dist)
        Assumes two chains, i.e. a unique interface set per chain.
    '''
    select_ats = []
    for at in st.get_atoms():
        # Skip Hydrogens to reduce time
        if at.element != 'H':
            select_ats.append(at)
    nbsearch = NeighborSearch(select_ats)
    interface = {}
    # Sets are more efficient than lists. Use sets when order is not relevant
    for ch in st[0]:
        interface[ch.id] = set()

    for at1, at2 in nbsearch.search_all(dist):
        #Only different chains
        res1 = at1.get_parent()
        ch1 = res1.get_parent()
        res2 = at2.get_parent()
        ch2 = res2.get_parent()
        if ch1 != ch2:
            interface[ch1.id].add(res1)
            interface[ch2.id].add(res2)
    return interface
get_interface(st, 4)
def calc_int_energies(st, res):
    '''Returns interaction energies (residue against other chains)
        for all atoms and for Ala atoms
    '''
    elec = 0.
    elec_ala = 0.
    vdw = 0.
    vdw_ala = 0.

    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent():
                r = at1 - at2
                e = elec_int(at1, at2, r)
                elec += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    elec_ala += e
                e = vdw_int(at1, at2, r)
                vdw += e
                if at1.id in ala_atoms: #GLY are included implicitly
                    vdw_ala += e
    return elec, elec_ala, vdw, vdw_ala


#Energy of interface residues
io = PDBIO()
st_chains = {}
NACCESS_BINARY = '/home/juliapm/biofi/soft/NACCESS/naccess'

add_atom_parameters(st, residue_library, ff_params)
srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

# State cut off distance between residues on the interface
MAXDIST = 4


# Iterate over structure to get chains
class SelectChain(Select):
    def __init__(self, chid):
        self.id = chid

    def accept_chain(self, chain):
        if chain.id == self.id:
            return 1
        else:
            return 0

for ch in st[0]:
    io.set_structure(st)
    io.save('tmp.pdb', SelectChain(ch.id))
    st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
    add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
    srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
    os.remove('tmp.pdb')

# If cut off distance = 0 we would take all residues on the structure,
# Therefore if its bigger call function and get residues on interface
if MAXDIST > 0.:
    interface = get_interface(st, MAXDIST)
    

    
# Create dictionaries to have count for every chain individually
elec = {}
elec_ala = {}

vdw = {}
vdw_ala = {}

solvAB = {}
solvAB_ala = {}

solvA = {}
solvA_ala = {}

# Create variables to have count of the total
total_elec = 0.
total_vdw = 0.
total_solv = 0.
total_solvMon = {}


# Get chains ID, because not always will be A and B
# and create a variable to have count of the total solvation in each chain.
chids = []
for ch in st[0]:
    chids.append(ch.id)
    total_solvMon[ch.id] = 0

# Total counter for total binding energy
total = 0.

# get interaction energy between chains A and E of all residues
# write result on a file 
with open("result_dG.txt", "a") as res_dG:
    for ch in st[0]:
        for res in ch.get_residues():
            # If cut off distance is bigger than 0 and not in the interface -> don't sum and call functions
            if MAXDIST > 0 and res not in interface[ch.id]:
                continue
            # Get every value with return value of functions if previous conditions are accomplished (in interface)
            elec[res], elec_ala[res], vdw[res], vdw_ala[res] = calc_int_energies(st[0], res)
            solvAB[res], solvAB_ala[res] = calc_solvation(st[0], res)
            solvA[res], solvA_ala[res] = calc_solvation(st_chains[ch.id], st_chains[ch.id][0][ch.id][res.id[1]])

            # Add every value of iteration
            total_elec += elec[res]
            total_vdw += vdw[res]
            total_solv += solvAB[res]
            total_solvMon[ch.id] += solvA[res]
            # This follows formula for binding energy between two chains
            # Count for every chain individually minus the solvent difference (solvAB[res])
            total += elec[res] + vdw[res] + solvAB[res] - solvA[res]

    print('Total Electrostatic: ', total_elec, file = res_dG)
    print('Total Vdw: ', total_vdw, file = res_dG)
    print('Total Solvation Complex: ', total_solv, file = res_dG)
    print('Total Solv ', chids[0], total_solvMon[chids[0]], file = res_dG)
    print('Total Solv ', chids[1], total_solvMon[chids[1]], file = res_dG)
    print('DG AB-A-B', total, file = res_dG)

# Same script as before but in order to calculate total binding energy, 
# not only in the interface, the only difference will be in the MAXDIST variable, now set off to 0

io = PDBIO()
st_chains = {}
NACCESS_BINARY = '/home/juliapm/biofi/soft/NACCESS/naccess'

add_atom_parameters(st, residue_library, ff_params)
srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)


MAXDIST = 0


# Iterate over structure to get chains
class SelectChain(Select):
    def __init__(self, chid):
        self.id = chid

    def accept_chain(self, chain):
        if chain.id == self.id:
            return 1
        else:
            return 0

for ch in st[0]:
    io.set_structure(st)
    io.save('tmp.pdb', SelectChain(ch.id))
    st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
    add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
    srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
os.remove('tmp.pdb')

if MAXDIST > 0.:
    interface = get_interface(st, MAXDIST)
    
# Create dictionaries to have count for every chain individually
elec = {}
elec_ala = {}

vdw = {}
vdw_ala = {}

solvAB = {}
solvAB_ala = {}

solvA = {}
solvA_ala = {}

# Create variables to have count of the total
total_elec = 0.
total_vdw = 0.
total_solv = 0.
total_solvMon = {}



# Get chains ID
# and create a variable to have count of the total solvation in each chain.
chids = []
for ch in st[0]:
    chids.append(ch.id)
    total_solvMon[ch.id] = 0

total = 0.

# Get interaction energy of all residues, the ones in the interface and the ones outside.
with open("res_total_e.txt", "a") as res_total:
    for ch in st[0]:
        for res in ch.get_residues():
            if MAXDIST > 0 and res not in interface[ch.id]:
                continue
            # Get every value with return value of functions
            elec[res], elec_ala[res], vdw[res], vdw_ala[res] = calc_int_energies(st[0], res)
            solvAB[res], solvAB_ala[res] = calc_solvation(st[0], res)
            solvA[res], solvA_ala[res] = calc_solvation(st_chains[ch.id], st_chains[ch.id][0][ch.id][res.id[1]])

            # Add every value of iteration
            total_elec += elec[res]
            total_vdw += vdw[res]
            total_solv += solvAB[res]
            total_solvMon[ch.id] += solvA[res]
            # Calculate the formula
            total += elec[res] + vdw[res] + solvAB[res] - solvA[res]


    print('Total Electrostatic: ', total_elec, file = res_total)
    print('Total Vdw: ', total_vdw, file = res_total)
    print('Total Solvation Complex: ', total_solv, file = res_total)
    print('Total Solv ', chids[0], total_solvMon[chids[0]], file = res_total)
    print('Total Solv ', chids[1], total_solvMon[chids[1]], file = res_total)
    print('DG AB-A-B', total, file = res_total)
