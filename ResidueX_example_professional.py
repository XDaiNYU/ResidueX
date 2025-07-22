#!/usr/bin/env python3
"""
Professional Example: ResidueX ncAA Peptide Generation

This script demonstrates the professional API for ResidueX,
showing how to perform non-canonical amino acid substitutions
in peptides using the clean, object-oriented interface.

Author: Xuhang Dai, Dr. Rui Wang
License: MIT
"""

import os
import sys
from pathlib import Path

# Import the professional ResidueX API
from ResidueX.api import ResidueXProcessor, ResidueXUtils, quick_ncaa_substitution


def main():
    """Main function demonstrating ResidueX professional API."""
    
    print("=" * 60)
    print("ResidueX Professional API Example")
    print("Non-canonical Amino Acid Peptide Generation")
    print("=" * 60)
    
    # Set up paths
    example_dir = Path("./example/6ox2_Z")
    peptide_pdb = example_dir / "ranked_100_sp_pep.pdb"
    
    # Check if example data exists
    if not peptide_pdb.exists():
        print(f"❌ Example data not found: {peptide_pdb}")
        print("Please ensure the example data is available.")
        return
    
    print(f"✅ Using peptide: {peptide_pdb}")
    
    # Example 1: Using the ResidueXProcessor class
    print("\n" + "=" * 40)
    print("Example 1: Using ResidueXProcessor Class")
    print("=" * 40)
    
    # Initialize processor
    processor = ResidueXProcessor(working_directory="./output_professional")
    
    # Define ncAA parameters
    ncaa_smiles = 'CN[C@@H](CC1=CN(C)C=N1)C(=O)C'  # Example ncAA
    residue_id = 8
    
    print(f"📋 Parameters:")
    print(f"   - ncAA SMILES: {ncaa_smiles}")
    print(f"   - Residue ID: {residue_id}")
    print(f"   - Working Directory: {processor.working_directory}")
    
    try:
        # Perform complete ncAA substitution workflow
        print("\n🔄 Processing ncAA substitution...")
        results = processor.process_ncaa_substitution(
            peptide_pdb=str(peptide_pdb),
            ncaa_smiles=ncaa_smiles,
            residue_id=residue_id
        )
        
        print("✅ Processing completed successfully!")
        print(f"\n📁 Output Directory: {results['output_directory']}")
        print(f"📄 Generated {len(results['conformers'])} conformers")
        print(f"🔗 Created {len(results['integration_results'])} integrations")
        
        # Show some results
        if results['integration_results']:
            first_result = results['integration_results'][0]
            print(f"📏 Minimum distance: {first_result['min_distance']:.2f} Å")
        
    except Exception as e:
        print(f"❌ Error during processing: {e}")
        return
    
    # Example 2: Using convenience functions
    print("\n" + "=" * 40)
    print("Example 2: Using Convenience Functions")
    print("=" * 40)
    
    try:
        # Quick ncAA substitution
        print("🔄 Performing quick ncAA substitution...")
        quick_results = quick_ncaa_substitution(
            peptide_pdb=str(peptide_pdb),
            ncaa_smiles=ncaa_smiles,
            residue_id=residue_id,
            output_dir="./output_quick"
        )
        
        print("✅ Quick substitution completed!")
        print(f"📁 Output: {quick_results['output_directory']}")
        
    except Exception as e:
        print(f"❌ Error in quick substitution: {e}")
    
    # Example 3: Using utility functions
    print("\n" + "=" * 40)
    print("Example 3: Using Utility Functions")
    print("=" * 40)
    
    try:
        # Use utility functions
        utils = ResidueXUtils()
        
        # Pair files in directory
        print("📂 Pairing files in example directory...")
        paired_files = utils.pair_files_in_directory(str(example_dir))
        print(f"✅ Found {len(paired_files)} file pairs")
        
        # Show first few pairs
        for i, pair in enumerate(paired_files.head(3).itertuples()):
            print(f"   {i+1}. {pair.file1} ↔ {pair.file2}")
        
    except Exception as e:
        print(f"❌ Error in utility functions: {e}")
    
    # Example 4: Step-by-step workflow
    print("\n" + "=" * 40)
    print("Example 4: Step-by-Step Workflow")
    print("=" * 40)
    
    try:
        step_processor = ResidueXProcessor(working_directory="./output_stepwise")
        
        # Step 1: Split peptide
        print("1️⃣ Splitting peptide by residue...")
        split_results = step_processor.split_peptide_by_residue(
            input_pdb=str(peptide_pdb),
            residue_id=residue_id
        )
        print(f"   ✅ Split complete: {split_results['residue_file']}")
        
        # Step 2: Generate conformers
        print("2️⃣ Generating ncAA conformers...")
        conformers = step_processor.generate_ncaa_conformers(
            ncaa_smiles=ncaa_smiles,
            reference_pdb=split_results['residue_file']
        )
        print(f"   ✅ Generated {len(conformers)} conformers")
        
        # Step 3: Integrate first conformer
        if conformers:
            print("3️⃣ Integrating first conformer...")
            integration = step_processor.integrate_ncaa_into_peptide(
                ncaa_sdf=conformers[0],
                peptide_pdb=str(peptide_pdb),
                residue_id=residue_id
            )
            print(f"   ✅ Integration complete: {integration['modified_peptide_sdf']}")
            print(f"   📏 Minimum distance: {integration['min_distance']:.2f} Å")
        
    except Exception as e:
        print(f"❌ Error in step-by-step workflow: {e}")
    
    print("\n" + "=" * 60)
    print("🎉 ResidueX Professional API Example Completed!")
    print("=" * 60)
    print("\n📚 For more information, see:")
    print("   - README.md: Installation and usage instructions")
    print("   - tutorial_ncAA_peptide_generation.ipynb: Detailed tutorial")
    print("   - https://github.com/XDaiNYU/ResidueX")
    print("\n📧 Contact: xd638@nyu.edu")


if __name__ == "__main__":
    main() 