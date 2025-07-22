#!/usr/bin/env python3
"""
ResidueX Example: Non-canonical Amino Acid Peptide Generation

This script demonstrates the complete workflow for generating ncAA peptides
using the ResidueX package. It performs the following steps:
1. Split a peptide PDB file by specific residue ID
2. Generate conformers for a non-canonical amino acid (ncAA)
3. Integrate ncAA structures into the peptide backbone
4. Convert molecular formats and calculate distances

Author: Xuhang Dai, Dr. Rui Wang
License: MIT
"""

import glob
import os
import pandas as pd
from ResidueX.residuex import *
from ResidueX.generate_conformer import *
from ResidueX.utils import *


def main():
    """
    Main function demonstrating the ResidueX workflow.
    """
    print("=" * 60)
    print("ResidueX: Non-canonical Amino Acid Peptide Generation")
    print("=" * 60)
    
    # ============================================================================
    # Configuration Parameters
    # ============================================================================
    
    # File paths and directories
    EXAMPLE_DIR = './example/'
    CASE_NAME = '6ox2_Z'
    WORK_PATH = f"{EXAMPLE_DIR}/{CASE_NAME}"
    
    # ncAA parameters
    NCAA_SMILES = 'CN[C@@H](CC1=CN(C)C=N1)C(=O)C'  # Example ncAA structure
    NCAA_NAME = 'PDB_NCAA'
    RESIDUE_ID = 8  # Residue to substitute
    
    # Input/output file names
    CHOSEN_PEPTIDE = 'ranked_100_sp_pep.pdb'
    INPUT_PDB_PATH = f"{WORK_PATH}/{CHOSEN_PEPTIDE}"
    OUTPUT_TAG = CHOSEN_PEPTIDE.split('.')[0]
    
    print(f"📁 Working Directory: {WORK_PATH}")
    print(f"🧬 ncAA SMILES: {NCAA_SMILES}")
    print(f"📍 Residue ID: {RESIDUE_ID}")
    print(f"📄 Input Peptide: {CHOSEN_PEPTIDE}")
    
    # ============================================================================
    # Step 1: Pair files in directory
    # ============================================================================
    
    print("\n" + "=" * 40)
    print("Step 1: Pairing files in directory")
    print("=" * 40)
    
    df_all_work = pd.DataFrame()
    df_each_pdb_case = pair_files_in_directory(WORK_PATH)
    df_all_work = pd.concat([df_all_work, df_each_pdb_case], ignore_index=True)
    
    print(f"✅ Found {len(df_each_pdb_case)} file pairs")
    
    # ============================================================================
    # Step 2: Split peptide by residue
    # ============================================================================
    
    print("\n" + "=" * 40)
    print("Step 2: Splitting peptide by residue")
    print("=" * 40)
    
    # Create output directory
    OUTPUT_PDB_PATH = f"{WORK_PATH}/{OUTPUT_TAG}"
    os.makedirs(OUTPUT_PDB_PATH, exist_ok=True)
    
    # Define output file paths
    OUTPUT_RESIDUE_PATH = f"{OUTPUT_PDB_PATH}/residue.pdb"
    OUTPUT_REST_PATH = f"{OUTPUT_PDB_PATH}/rest.pdb"
    
    # Split the peptide
    split_pdb_by_residue(
        input_pdb_file=INPUT_PDB_PATH,
        residue_id=RESIDUE_ID,
        output_residue_file=OUTPUT_RESIDUE_PATH,
        output_rest_file=OUTPUT_REST_PATH
    )
    
    print(f"✅ Split complete:")
    print(f"   - Residue file: {OUTPUT_RESIDUE_PATH}")
    print(f"   - Rest of peptide: {OUTPUT_REST_PATH}")
    
    # ============================================================================
    # Step 3: Generate ncAA conformers
    # ============================================================================
    
    print("\n" + "=" * 40)
    print("Step 3: Generating ncAA conformers")
    print("=" * 40)
    
    # Create ncAA SDF output directory
    NCAA_SDF_SAVED_PATH = f"{OUTPUT_PDB_PATH}/NCAA_conf_SDF_{NCAA_NAME}"
    os.makedirs(NCAA_SDF_SAVED_PATH, exist_ok=True)
    
    # Generate conformers
    NCAA_sdf_generation(
        NCAA_sml=NCAA_SMILES,
        input_pdb=OUTPUT_RESIDUE_PATH,
        saved_path=NCAA_SDF_SAVED_PATH
    )
    
    print(f"✅ Conformers generated in: {NCAA_SDF_SAVED_PATH}")
    
    # ============================================================================
    # Step 4: Integrate ncAA into peptide
    # ============================================================================
    
    print("\n" + "=" * 40)
    print("Step 4: Integrating ncAA into peptide")
    print("=" * 40)
    
    # Get list of generated SDF files
    sdf_files_list = glob.glob(os.path.join(NCAA_SDF_SAVED_PATH, '*.sdf'))
    print(f"📊 Processing {len(sdf_files_list)} conformers...")
    
    # Initialize data collection
    csv_data = []
    list_fail = []
    
    # Process each conformer
    for each_ncaa_sdf_path in sdf_files_list:
        try:
            # Prepare file paths
            chosen_peptide_ready = f"{CHOSEN_PEPTIDE}_ready"
            output_pep_pdb_path = f"{WORK_PATH}/{chosen_peptide_ready}"
            
            each_ncaa_sdf_name = os.path.basename(each_ncaa_sdf_path).split('.')[0]
            output_ncaa_pep_path = f"{WORK_PATH}/{OUTPUT_TAG}/output_NCAA_peptides_{NCAA_NAME}"
            os.makedirs(output_ncaa_pep_path, exist_ok=True)
            
            # Define output file paths
            ncaa_peptide_sdf_path = f"{output_ncaa_pep_path}/ncaa_res{RESIDUE_ID}_{each_ncaa_sdf_name}.sdf"
            ncaa_peptide_pdb_path = f"{output_ncaa_pep_path}/ncaa_res{RESIDUE_ID}_{each_ncaa_sdf_name}.pdb"
            distance_info_txt_path = f"{output_ncaa_pep_path}/min_distance_ncaa_res{RESIDUE_ID}_{each_ncaa_sdf_name}.txt"
            
            # Integrate ncAA into peptide
            min_distance = integrate_NCAA_into_peptide(
                each_ncaa_SDF_path=each_ncaa_sdf_path,
                residue_id=RESIDUE_ID,
                input_pep_pdb_path=INPUT_PDB_PATH,
                output_pep_pdb_path=output_pep_pdb_path,
                NCAA_peptide_SDF_path=ncaa_peptide_sdf_path,
                distance_info_txt_path=distance_info_txt_path
            )
            
            # Convert SDF to PDB using Open Babel
            os.system(f'obabel -isdf {ncaa_peptide_sdf_path} -O {ncaa_peptide_pdb_path}')
            
            # Collect results
            csv_data.append({
                'conformer': each_ncaa_sdf_name,
                'min_distance': min_distance,
                'sdf_file': ncaa_peptide_sdf_path,
                'pdb_file': ncaa_peptide_pdb_path
            })
            
            print(f"   ✅ {each_ncaa_sdf_name}")
            
        except Exception as e:
            print(f"   ❌ {each_ncaa_sdf_name}: Failed - {e}")
            list_fail.append(each_ncaa_sdf_name)
    
    # ============================================================================
    # Summary
    # ============================================================================
    
    print("\n" + "=" * 40)
    print("Summary")
    print("=" * 40)
    
    print(f"🎯 Total conformers processed: {len(sdf_files_list)}")
    print(f"✅ Successful integrations: {len(csv_data)}")
    print(f"❌ Failed integrations: {len(list_fail)}")
    
    if csv_data:
        print(f"✅ All conformers processed successfully")
    
    print(f"\n📁 Output directory: {output_ncaa_pep_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()





