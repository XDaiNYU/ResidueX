"""
ResidueX API: High-level interface for non-canonical amino acid peptide generation.

This module provides a clean, professional API for the ResidueX package,
wrapping the core functionality in easy-to-use classes and methods.
"""

import os
import glob
from typing import List, Dict, Optional, Tuple, Union
from pathlib import Path

from .residuex import (
    split_pdb_by_residue,
    NCAA_sdf_generation,
    integrate_NCAA_into_peptide,
    read_sdf_to_Mol,
    filter_residue_atoms
)
from .generate_conformer import GenerateConformer
from .utils import pair_files_in_directory


class ResidueXProcessor:
    """
    Main processor class for ResidueX operations.
    
    Provides a high-level interface for non-canonical amino acid peptide generation.
    """
    
    def __init__(self, working_directory: str = "."):
        """
        Initialize the ResidueX processor.
        
        Args:
            working_directory (str): Base directory for operations
        """
        self.working_directory = Path(working_directory)
        self.working_directory.mkdir(exist_ok=True)
        
    def split_peptide_by_residue(
        self, 
        input_pdb: str, 
        residue_id: int, 
        output_dir: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Split a peptide PDB file by residue ID.
        
        Args:
            input_pdb (str): Path to input PDB file
            residue_id (int): Residue ID to extract
            output_dir (str, optional): Output directory (defaults to working directory)
            
        Returns:
            Dict[str, str]: Paths to residue and rest PDB files
        """
        if output_dir is None:
            output_dir = self.working_directory / f"residue_{residue_id}"
        else:
            output_dir = Path(output_dir)
            
        output_dir.mkdir(exist_ok=True)
        
        residue_file = output_dir / "residue.pdb"
        rest_file = output_dir / "rest.pdb"
        
        split_pdb_by_residue(
            input_pdb_file=str(input_pdb),
            residue_id=residue_id,
            output_residue_file=str(residue_file),
            output_rest_file=str(rest_file)
        )
        
        return {
            "residue_file": str(residue_file),
            "rest_file": str(rest_file),
            "output_dir": str(output_dir)
        }
    
    def generate_ncaa_conformers(
        self, 
        ncaa_smiles: str, 
        reference_pdb: str, 
        output_dir: Optional[str] = None,
        num_conformers: int = 13
    ) -> List[str]:
        """
        Generate conformers for a non-canonical amino acid.
        
        Args:
            ncaa_smiles (str): SMILES string of the ncAA
            reference_pdb (str): Reference PDB file for alignment
            output_dir (str, optional): Output directory for conformers
            num_conformers (int): Number of conformers to generate
            
        Returns:
            List[str]: List of generated SDF file paths
        """
        if output_dir is None:
            output_dir = self.working_directory / "ncAA_conformers"
        else:
            output_dir = Path(output_dir)
            
        output_dir.mkdir(exist_ok=True)
        
        NCAA_sdf_generation(ncaa_smiles, reference_pdb, str(output_dir))
        
        # Return list of generated SDF files
        sdf_files = glob.glob(str(output_dir / "*.sdf"))
        return sorted(sdf_files)
    
    def integrate_ncaa_into_peptide(
        self,
        ncaa_sdf: str,
        peptide_pdb: str,
        residue_id: int,
        output_dir: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Integrate a non-canonical amino acid into a peptide.
        
        Args:
            ncaa_sdf (str): Path to ncAA SDF file
            peptide_pdb (str): Path to peptide PDB file
            residue_id (int): Residue ID for integration
            output_dir (str, optional): Output directory
            
        Returns:
            Dict[str, str]: Paths to output files
        """
        if output_dir is None:
            output_dir = self.working_directory / f"ncaa_integration_res{residue_id}"
        else:
            output_dir = Path(output_dir)
            
        output_dir.mkdir(exist_ok=True)
        
        # Prepare output paths
        peptide_ready = peptide_pdb.replace('.pdb', '_ready.pdb')
        sdf_output = output_dir / f"ncaa_res{residue_id}.sdf"
        distance_info = output_dir / f"distance_info_res{residue_id}.txt"
        
        # Perform integration
        min_distance = integrate_NCAA_into_peptide(
            each_ncaa_SDF_path=ncaa_sdf,
            residue_id=residue_id,
            input_pep_pdb_path=peptide_pdb,
            output_pep_pdb_path=peptide_ready,
            NCAA_peptide_SDF_path=str(sdf_output),
            distance_info_txt_path=str(distance_info)
        )
        
        return {
            "modified_peptide_sdf": str(sdf_output),
            "distance_info": str(distance_info),
            "peptide_ready": peptide_ready,
            "min_distance": min_distance,
            "output_dir": str(output_dir)
        }
    
    def process_ncaa_substitution(
        self,
        peptide_pdb: str,
        ncaa_smiles: str,
        residue_id: int,
        output_dir: Optional[str] = None
    ) -> Dict[str, Union[str, List[str], float]]:
        """
        Complete workflow for ncAA substitution in a peptide.
        
        Args:
            peptide_pdb (str): Path to peptide PDB file
            ncaa_smiles (str): SMILES string of the ncAA
            residue_id (int): Residue ID to substitute
            output_dir (str, optional): Output directory
            
        Returns:
            Dict: Complete results including all generated files
        """
        if output_dir is None:
            output_dir = self.working_directory / f"complete_substitution_res{residue_id}"
        else:
            output_dir = Path(output_dir)
            
        output_dir.mkdir(exist_ok=True)
        
        # Step 1: Split peptide by residue
        split_results = self.split_peptide_by_residue(
            input_pdb=peptide_pdb,
            residue_id=residue_id,
            output_dir=output_dir / "split"
        )
        
        # Step 2: Generate ncAA conformers
        conformers = self.generate_ncaa_conformers(
            ncaa_smiles=ncaa_smiles,
            reference_pdb=split_results["residue_file"],
            output_dir=output_dir / "conformers"
        )
        
        # Step 3: Integrate each conformer
        integration_results = []
        for conformer in conformers:
            result = self.integrate_ncaa_into_peptide(
                ncaa_sdf=conformer,
                peptide_pdb=peptide_pdb,
                residue_id=residue_id,
                output_dir=output_dir / "integrations" / Path(conformer).stem
            )
            integration_results.append(result)
        
        return {
            "split_results": split_results,
            "conformers": conformers,
            "integration_results": integration_results,
            "output_directory": str(output_dir)
        }


class ResidueXUtils:
    """Utility class for ResidueX operations."""
    
    @staticmethod
    def read_molecule_from_sdf(sdf_path: str):
        """Read a molecule from SDF file."""
        return read_sdf_to_Mol(sdf_path)
    
    @staticmethod
    def filter_residue_atoms(
        input_pdb: str, 
        output_pdb: str, 
        residue_id: int, 
        allowed_atoms: List[str]
    ):
        """Filter specific atoms from a residue."""
        filter_residue_atoms(input_pdb, output_pdb, residue_id, allowed_atoms)
    
    @staticmethod
    def pair_files_in_directory(directory_path: str):
        """Pair files in a directory."""
        return pair_files_in_directory(directory_path)


# Convenience functions for quick operations
def quick_ncaa_substitution(
    peptide_pdb: str,
    ncaa_smiles: str,
    residue_id: int,
    output_dir: str = "."
) -> Dict[str, Union[str, List[str], float]]:
    """
    Quick function for ncAA substitution.
    
    Args:
        peptide_pdb (str): Path to peptide PDB file
        ncaa_smiles (str): SMILES string of the ncAA
        residue_id (int): Residue ID to substitute
        output_dir (str): Output directory
        
    Returns:
        Dict: Results of the substitution process
    """
    processor = ResidueXProcessor(output_dir)
    return processor.process_ncaa_substitution(
        peptide_pdb=peptide_pdb,
        ncaa_smiles=ncaa_smiles,
        residue_id=residue_id
    )


def generate_ncaa_conformers(
    ncaa_smiles: str,
    reference_pdb: str,
    output_dir: str = "."
) -> List[str]:
    """
    Quick function for generating ncAA conformers.
    
    Args:
        ncaa_smiles (str): SMILES string of the ncAA
        reference_pdb (str): Reference PDB file
        output_dir (str): Output directory
        
    Returns:
        List[str]: List of generated SDF file paths
    """
    processor = ResidueXProcessor(output_dir)
    return processor.generate_ncaa_conformers(
        ncaa_smiles=ncaa_smiles,
        reference_pdb=reference_pdb
    ) 