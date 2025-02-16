#YES
from Bio import PDB
import os
def split_pdb_by_residue(input_pdb_file, residue_id, output_residue_file, output_rest_file):
    """
    Splits a PDB file into two based on the given residue id.

    Parameters:
    - input_pdb_file (str): Path to the input PDB file.
    - residue_id (int): The residue ID to use for the split.
    - output_residue_file (str): Path to save the PDB file containing just the given residue.
    - output_rest_file (str): Path to save the PDB file containing all residues except the given one.
    """

    # Create a parser and load the PDB structure from the file
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb_file)

    # Create two IO objects for saving the PDB files
    io_residue = PDB.PDBIO()
    io_rest = PDB.PDBIO()

    # Select only the given residue
    class ResidueSelect(PDB.Select):
        def accept_residue(self, residue):
            return residue.id[1] == residue_id

    # Select everything except the given residue
    class RestSelect(PDB.Select):
        def accept_residue(self, residue):
            return residue.id[1] != residue_id

    # Save the two PDB files
    io_residue.set_structure(structure)
    io_residue.save(output_residue_file, ResidueSelect())

    io_rest.set_structure(structure)
    io_rest.save(output_rest_file, RestSelect())

#YES
import sys, os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, PandasTools
import pandas as pd
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

import matplotlib.pyplot as plt
import sklearn#scikt-learn
import rmsd #rmsd-1.4

from ResidueX.generate_conformer import GenerateConformer ## import custom function
def NCAA_sdf_generation(NCAA_sml,input_pdb,saved_path):



    # Create a molecule object from the SMILES string
    new_mol = Chem.MolFromSmiles(NCAA_sml)

    # Sanitize the molecule
    Chem.SanitizeMol(new_mol)
    
    

    
    
    
    # Define the SMARTS pattern for the common backbone of an amino acid
    single_aa_backbone_patt = Chem.MolFromSmarts('NCC(=O)')


    # Find the backbone atoms in the molecule
    backbone_atoms = new_mol.GetSubstructMatch(single_aa_backbone_patt)

    # Get all atom indices in the molecule
    all_atom_indices = set(range(new_mol.GetNumAtoms()))

    
    
    # Identifying the atoms connected to the carbonyl carbon (C=O) and filtering out backbone atoms
    connected_to_CO = []
    carbonyl_carbon_index = backbone_atoms[2]
    for neighbor in new_mol.GetAtomWithIdx(carbonyl_carbon_index).GetNeighbors():
        if neighbor.GetIdx() not in backbone_atoms and (neighbor.GetSymbol() != 'O' or neighbor.GetTotalValence() != 2):
            connected_to_CO.append(neighbor.GetIdx())

    # Identifying the atoms connected to the nitrogen (N) and filtering out backbone atoms
    connected_to_N = []
    nitrogen_index = backbone_atoms[0]
    for neighbor in new_mol.GetAtomWithIdx(nitrogen_index).GetNeighbors():
        if neighbor.GetIdx() not in backbone_atoms:
            connected_to_N.append(neighbor.GetIdx())

    connected_to_CO, connected_to_N

    
    
    
    
    
    
    
    
    # Subtract the set of backbone atom indices from the set of all atom indices
    # The result is the set of side-chain atom indices
    side_chain_atom_indices = all_atom_indices - set(backbone_atoms)-set(connected_to_CO)-set(connected_to_N)

    # Convert the indices back to a sorted list
    side_chain_atom_indices = sorted(list(side_chain_atom_indices))

    side_chain_atom_indices
    
    

    #refMol = Chem.MolFromMol2File('data/2LXS_lig.mol2', removeHs = True) ###reference peptide with available conformer
    refMol = Chem.MolFromPDBFile(input_pdb, removeHs = True)



    Chem.MolToSmiles(new_mol)

    refMol = Chem.MolFromPDBFile(input_pdb, removeHs = True)
    
    print(side_chain_atom_indices)
    GenerateConformer(new_mol, refMol, saved_path, side_chain_atom_indices, optimize=False)

    
    
def read_sdf_to_Mol(sdf_file):
    # Create an SDMolSupplier object
    supplier = Chem.SDMolSupplier(sdf_file, sanitize=False)
    # Loop through the supplier to access each molecule
    for mol in supplier:
        if mol is not None:  # Check if the molecule is successfully read
            # Do something with the molecule
            print(Chem.MolToSmiles(mol))
        else:
            print("A molecule in the SDF could not be read.")
    # If you just want the first molecule, you can directly access it
    first_mol = supplier[0] if supplier[0] is not None else None
    

    
    return first_mol




def filter_residue_atoms(input_pdb_path, output_pdb_path, residue_id, allowed_atoms):
    with open(input_pdb_path, 'r') as input_file, open(output_pdb_path, 'w') as output_file:
        for line in input_file:
            # Check if line is an atom record
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract the residue sequence number
                current_residue_id = int(line[22:26].strip())
                
                # If it's the specified residue, filter by allowed atoms
                if current_residue_id == residue_id:
                    atom_name = line[12:16].strip()
                    if atom_name in allowed_atoms:
                        output_file.write(line)
                else:
                    output_file.write(line)
            else:
                # Write non-atom lines (like headers) without changes
                output_file.write(line)


from rdkit import Chem
from rdkit.Chem import AllChem

def read_pdb_to_mol(pdb_path):
    """
    Reads a PDB file and returns it as an RDKit Mol object.
    """
    return Chem.MolFromPDBFile(pdb_path, removeHs=False)


#2023 12 11 
# *add carbon on CO and N variables
def get_pep_ready_carbon_alpha(pep_ready,residue_id):
    patt = Chem.MolFromSmiles('NCC(=O)')
    matches = pep_ready.GetSubstructMatches(patt)
    pep_carbon_alpha_to_link=matches[residue_id-1][1]#Accornign to 'NCC(=O)', the order should be (68, 70, 72, 73)
    
    pep_carbon_CO_to_link=matches[residue_id-1][2]
    pep_carbon_N_to_link=matches[residue_id-1][0]
    
    
    pep_ready_specific_residue_list=list(matches[residue_id-1])
    
    
        # Find hydrogen atoms connected to the atoms in the specific residue
    hydrogen_indices = []
    for idx in pep_ready_specific_residue_list:
        atom = pep_ready.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                hydrogen_indices.append(neighbor.GetIdx())

    return pep_carbon_alpha_to_link, pep_carbon_CO_to_link, pep_carbon_N_to_link, pep_ready_specific_residue_list, hydrogen_indices


def ncaa_atom_to_connect_and_delete(ncaa):
    from rdkit import Chem

    # Assuming ncaa is your molecule read from SDF
    # ncaa = read_sdf_to_Mol(each_ncaa_SDF_path)

    # Define the SMARTS pattern for the common backbone of an amino acid
    patt = Chem.MolFromSmiles('NCC(=O)')

    # Find matches of the pattern in the molecule
    matches = ncaa.GetSubstructMatches(patt)

    # Extract the alpha carbon index from the first match
    ncaa_carbon_alpha = matches[0][1]  # Alpha carbon is the second atom in 'NCC(=O)'

    # List of atom indices in the backbone
    bb_list = list(matches[0])

    # Lists to store connected atoms
    connected_to_N = []
    connected_to_CO = []

    # Find atoms connected to N and C=O
    for i, atom_idx in enumerate(bb_list):
        atom = ncaa.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'N' or (atom.GetSymbol() == 'C' and i == 2):  # Nitrogen or Carbonyl Carbon
            for neighbor in atom.GetNeighbors():
                # Exclude the backbone atoms and the alpha carbon
                if neighbor.GetIdx() not in bb_list and neighbor.GetIdx() != ncaa_carbon_alpha:
                    if atom.GetSymbol() == 'N':
                        connected_to_N.append(neighbor.GetIdx())
                    elif atom.GetSymbol() == 'C':
                        connected_to_CO.append(neighbor.GetIdx())

    # Find atoms connected to the atoms connected to N and C=O
    connected_to_connected_N = []
    connected_to_connected_CO = []

    for atom_idx in connected_to_N:
        for neighbor in ncaa.GetAtomWithIdx(atom_idx).GetNeighbors():
            if neighbor.GetIdx() not in bb_list and neighbor.GetIdx() not in connected_to_N:
                connected_to_connected_N.append(neighbor.GetIdx())

    for atom_idx in connected_to_CO:
        for neighbor in ncaa.GetAtomWithIdx(atom_idx).GetNeighbors():
            if neighbor.GetIdx() not in bb_list and neighbor.GetIdx() not in connected_to_CO:
                connected_to_connected_CO.append(neighbor.GetIdx())

    # All atom indices list
    all_atom_indices_list = [atom.GetIdx() for atom in ncaa.GetAtoms()]

    # Exclude the specified atoms
    excluded_indices = set(connected_to_N + connected_to_CO + connected_to_connected_N + connected_to_connected_CO + bb_list)-set([ncaa_carbon_alpha])
    remaining_indices = [idx for idx in all_atom_indices_list if idx not in excluded_indices]

    #remaining_indices
    
    atom_side_chain_list=remaining_indices
    atom_to_delete_list=list(excluded_indices)
    atom_side_chain_to_connect=ncaa_carbon_alpha
    
    
    ##Hydrogen atoms on C-alpha(this will be used to avoid being calculated in min distance)
    
        # Find hydrogen atoms connected to the alpha carbon
    connected_hydrogens_to_alpha_C = []

    # Get the alpha carbon atom
    alpha_C_atom = ncaa.GetAtomWithIdx(ncaa_carbon_alpha)

    # Iterate through neighbors of the alpha carbon
    for neighbor in alpha_C_atom.GetNeighbors():
        # Check if the neighbor is a hydrogen atom
        if neighbor.GetSymbol() == 'H':
            connected_hydrogens_to_alpha_C.append(neighbor.GetIdx())

    return atom_side_chain_list,atom_to_delete_list,atom_side_chain_to_connect,connected_hydrogens_to_alpha_C


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import rdchem

def integrate_NCAA_into_peptide(each_ncaa_SDF_path,residue_id, input_pep_pdb_path,output_pep_pdb_path,NCAA_peptide_SDF_path,distance_info_txt_path):

    ncaa=read_sdf_to_Mol(each_ncaa_SDF_path)
    allowed_atoms = ['N', 'H', 'CA', 'C', 'O']
#     allowed_atoms = ['N', 'H',  'C', 'O']
    filter_residue_atoms(input_pep_pdb_path, output_pep_pdb_path, residue_id, allowed_atoms)
    print(input_pep_pdb_path,output_pep_pdb_path,residue_id, allowed_atoms)
    pep_ready=read_pdb_to_mol(output_pep_pdb_path)

    pep_carbon_alpha_to_link,pep_carbon_CO_to_link, pep_carbon_N_to_link,pep_ready_specific_residue_list, pep_ready_specific_residue_hydrogen_indices=get_pep_ready_carbon_alpha(pep_ready,residue_id)
    
    print('pep_carbon_alpha_to_link,pep_carbon_CO_to_link, pep_carbon_N_to_link,pep_ready_specific_residue_list, pep_ready_specific_residue_hydrogen_indices',pep_carbon_alpha_to_link,pep_carbon_CO_to_link, pep_carbon_N_to_link,pep_ready_specific_residue_list, pep_ready_specific_residue_hydrogen_indices)

    atom_side_chain_list,atom_to_delete_list,atom_side_chain_to_connect,atom_connected_hydrogens_to_alpha_C=ncaa_atom_to_connect_and_delete(ncaa)

    print('atom_side_chain_list,atom_to_delete_list,atom_side_chain_to_connect,atom_connected_hydrogens_to_alpha_C',atom_side_chain_list,atom_to_delete_list,atom_side_chain_to_connect,atom_connected_hydrogens_to_alpha_C)







    combined_pep=Chem.CombineMols(pep_ready, ncaa)

    # The number of atoms in pepbb before combining
    num_atoms_pep_ready = pep_ready.GetNumAtoms()

    # Mapping for pepbb remains the same, since it's the first molecule
    mapping_pep_ready = {i: i for i in range(num_atoms_pep_ready)}

    # For sdfrag, the mapping will increment the atom indices by the number of atoms in pepbb
    mapping_ncaa = {i: i + num_atoms_pep_ready for i in range(ncaa.GetNumAtoms())}
    
    # To get the original atom number and current atom number you can print or store the mappings
    print("Mapping for pepbb:", mapping_pep_ready)
    print("Mapping for ncaa:", mapping_ncaa)





    # Update the atom_side_chain_list with the new indices
    atom_side_chain_list_updated = [mapping_ncaa[index] for index in atom_side_chain_list]

    # Update the atom_to_delete_list with the new indices
    # This assumes that the atom_to_delete_list contains indices from the ncaa molecule
    atom_to_delete_list_updated = [mapping_ncaa.get(index, index) for index in atom_to_delete_list]

    # Update the atom_side_chain_to_connect with the new index
    atom_side_chain_to_connect_updated = mapping_ncaa[atom_side_chain_to_connect] if atom_side_chain_to_connect is not None else None

    
    # Update the connected_hydrogens_to_alpha_C list with the new indices
    atom_connected_hydrogens_to_alpha_C_updated = [mapping_ncaa.get(h_atom_index, h_atom_index) for h_atom_index in atom_connected_hydrogens_to_alpha_C]

    # Now you have the updated lists with the correct indices after the merge
    print("Updated atom side-chain list:", atom_side_chain_list_updated)
    print("Updated atom to delete list:", atom_to_delete_list_updated)
    print("Updated atom side-chain to connect:", atom_side_chain_to_connect_updated)
    print("atom_connected_hydrogens_to_alpha_C_updated:",atom_connected_hydrogens_to_alpha_C_updated)









    ## Check whether clashes

    ncaa_side_chain_no_connection_atom_list=atom_side_chain_list_updated #including C-alpha
    pep_ready_all_atoms_list=list(mapping_pep_ready.keys())




    list1=ncaa_side_chain_no_connection_atom_list
    list1.remove(atom_side_chain_to_connect_updated)
    try:
        for i in atom_connected_hydrogens_to_alpha_C_updated:
            list1.remove(i)
    except:
        pass
    list2=pep_ready_all_atoms_list
    list2.remove(pep_carbon_alpha_to_link)
    print('list1',list1)
    print('list2',list2)
    # Calculate the distances and find the minimum
    min_distance = float('inf')
    min_distance_atoms = (None, None)  # Tuple to store indices of atoms with minimum distance

    # Get the conformer to access atom positions
    conf = combined_pep.GetConformer()

    # Calculate the distances and find the minimum
    min_distance = float('inf')
    for idx1 in list1:
        for idx2 in list2:
            # Calculate Euclidean distance between atoms idx1 and idx2
            pos1 = conf.GetAtomPosition(idx1)
            pos2 = conf.GetAtomPosition(idx2)
            distance = (pos1 - pos2).Length()  # Use rdGeometry method to get the distance
            if distance < min_distance:
                min_distance = distance
                min_distance_atoms = (idx1, idx2)
    # min_distance now contains the minimum distance between the atoms in the two lists
    print(f"The minimum distance between atoms in the two lists is: {min_distance:.2f} Å (between atoms {min_distance_atoms[0]} and {min_distance_atoms[1]})")


    # Open a text file in write mode
    with open(distance_info_txt_path, 'w') as file:
        # Write the formatted string to the file
            file.write(f"The minimum distance is: {min_distance:.2f} Å (between atoms {min_distance_atoms[0]} and {min_distance_atoms[1]})\n")


    editable_combined_pep = Chem.EditableMol(combined_pep)

    atom_index_1=pep_carbon_CO_to_link
    atom_index_2=atom_side_chain_to_connect_updated
    atom_index_3=pep_carbon_N_to_link
    
    #20240125
    fianl_atom_to_delete_list_updated=atom_to_delete_list_updated +[pep_carbon_alpha_to_link]
    print('fianl_atom_to_delete_list_updated',fianl_atom_to_delete_list_updated)
    # If you need to delete atoms, you can do so using the RemoveAtom method
    # Make sure you delete atoms in reverse order of their index (highest first) to avoid reindexing issues
    
    editable_combined_pep.AddBond(atom_index_1, atom_index_2, order=rdchem.BondType.SINGLE)
    editable_combined_pep.AddBond(atom_index_3, atom_index_2, order=rdchem.BondType.SINGLE)
    for atom_index in sorted(fianl_atom_to_delete_list_updated, reverse=True):#20230225
        editable_combined_pep.RemoveAtom(atom_index)



    print('atom_index_1,atom_index_2,atom_index_3:',atom_index_1,atom_index_2,atom_index_3)
    # Get the final molecule
    final_combined_pep = editable_combined_pep.GetMol()


    Chem.MolToMolFile(final_combined_pep, NCAA_peptide_SDF_path)

    
    return min_distance

def rename_UNL_pdb_file(file_path, target_residue_id):

    modified_lines = []

    with open(file_path, 'r') as file:
        for line in file:
            # Modify HETATM lines or lines with target residue ID
            if line.startswith("HETATM") or line[22:26].strip() == str(target_residue_id):
                # Extract fields
                record_type = 'HETATM' if line.startswith("HETATM") else 'ATOM  '
                atom_serial = line[6:11]
                atom_name = line[12:16]
                alt_loc = line[16]
                res_name = 'UNL'
                chain_id = 'B'
                res_seq = str(target_residue_id).rjust(4)
                icode = line[26]
                rest_of_line = line[27:76]
                element = line[76:78].strip().rjust(2)  # Right-align the element symbol

                # Reconstruct the line with proper spacing
                new_line = f"{record_type}{atom_serial} {atom_name}{alt_loc}{res_name} {chain_id}{res_seq}{icode}{rest_of_line}{element}\n"
                modified_lines.append(new_line)
            else:
                modified_lines.append(line)

    # Write the modified lines back to the file
    with open(file_path, 'w') as file:
        for line in modified_lines:
            file.write(line)
