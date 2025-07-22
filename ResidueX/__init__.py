"""
ResidueX: A Python package for non-canonical amino acid peptide generation.

This package provides tools for:
- Splitting PDB files by residue IDs
- Generating conformers for non-canonical amino acids (ncAAs)
- Integrating ncAA structures into peptide backbones
- Creating modified peptides with custom amino acid substitutions

Authors: Xuhang Dai, Dr. Rui Wang
License: MIT
"""

__version__ = "0.1.0"
__author__ = "Xuhang Dai, Dr. Rui Wang"
__email__ = "xd638@nyu.edu"

# Import core functions from residuex module
from .residuex import (
    split_pdb_by_residue,
    NCAA_sdf_generation,
    integrate_NCAA_into_peptide,
    read_sdf_to_Mol,
    filter_residue_atoms,
    read_pdb_to_mol,
    get_pep_ready_carbon_alpha,
    ncaa_atom_to_connect_and_delete,
    rename_UNL_pdb_file
)

# Import conformer generation functionality
from .generate_conformer import GenerateConformer

# Import utility functions
from .utils import pair_files_in_directory

# Import professional API
from .api import (
    ResidueXProcessor,
    ResidueXUtils,
    quick_ncaa_substitution,
    generate_ncaa_conformers
)

# Define what gets imported with "from ResidueX import *"
__all__ = [
    # Core functionality
    'split_pdb_by_residue',
    'NCAA_sdf_generation', 
    'integrate_NCAA_into_peptide',
    
    # Molecular handling
    'read_sdf_to_Mol',
    'read_pdb_to_mol',
    'filter_residue_atoms',
    
    # Peptide processing
    'get_pep_ready_carbon_alpha',
    'ncaa_atom_to_connect_and_delete',
    'rename_UNL_pdb_file',
    
    # Conformer generation
    'GenerateConformer',
    
    # Utilities
    'pair_files_in_directory',
    
    # Professional API
    'ResidueXProcessor',
    'ResidueXUtils',
    'quick_ncaa_substitution',
    'generate_ncaa_conformers'
]

# Package metadata
PACKAGE_INFO = {
    'name': 'ResidueX',
    'version': __version__,
    'description': 'Non-canonical amino acid peptide generation toolkit',
    'author': __author__,
    'email': __email__,
    'license': 'MIT',
    'url': 'https://github.com/XDaiNYU/ResidueX'
} 