import glob
import os
import pandas as pd
from ResidueX.residuex import *
from ResidueX.generate_conformer import *
from ResidueX.utils import *





pathE = '/example/'
df_all_work = pd.DataFrame()
list_fail = []

# Initialize a list to collect details for CSV
csv_data = []


# pdb_id = '6ox2'
# pep_chain_ID = 'Z'


case_name='6ox2_Z'


# pdb_case_path = f"{pathE}/{pdb_id}_{pep_chain_ID}"

pdb_case_path = f"{pathE}/{case_name}"


df_each_pdb_case = pair_files_in_directory(pdb_case_path)
df_all_work = pd.concat([df_all_work, df_each_pdb_case], ignore_index=True)

NCAA_sml = 'CN[C@@H](CC1=CN(C)C=N1)C(=O)C'
NCAA_name = 'PDB_NCAA'

residue_id = 8
work_path = pdb_case_path


chosen_peptide = 'ranked_100_sp_pep.pdb'
input_pdb_path = f"{work_path}/{chosen_peptide}"
residue_to_extract = residue_id
output_tag = chosen_peptide.split('.')[0]
output_pdb_path = f"{work_path}/{output_tag}"
os.makedirs(output_pdb_path, exist_ok=True)

output_residue_path = f"{output_pdb_path}/residue.pdb"
output_rest_path = f"{output_pdb_path}/rest.pdb"
split_pdb_by_residue(input_pdb_path, residue_to_extract, output_residue_path, output_rest_path)

ncaa_SDF_saved_path = f"{output_pdb_path}/NCAA_conf_SDF_{NCAA_name}"
os.makedirs(ncaa_SDF_saved_path, exist_ok=True)
NCAA_sdf_generation(NCAA_sml, output_residue_path, ncaa_SDF_saved_path)

sdf_files_list = glob.glob(os.path.join(ncaa_SDF_saved_path, '*.sdf'))
for each_ncaa_SDF_path in sdf_files_list:
    input_pep_pdb_path = f"{work_path}/{chosen_peptide}"
    chosen_peptide_ready = f"{chosen_peptide}_ready"
    output_pep_pdb_path = f"{work_path}/{chosen_peptide_ready}"

    each_ncaa_SDF_name = os.path.basename(each_ncaa_SDF_path).split('.')[0]
    output_NCAA_pep_path = f"{work_path}/{output_tag}/output_NCAA_peptides_{NCAA_name}"
    os.makedirs(output_NCAA_pep_path, exist_ok=True)

    NCAA_peptide_SDF_path = f"{output_NCAA_pep_path}/ncaa_res{residue_id}_{each_ncaa_SDF_name}.sdf"
    NCAA_peptide_PDB_path = f"{output_NCAA_pep_path}/ncaa_res{residue_id}_{each_ncaa_SDF_name}.pdb"
    distance_info_txt_path = f"{output_NCAA_pep_path}/min_distance_ncaa_res{residue_id}_{each_ncaa_SDF_name}.txt"

    min_distance = integrate_NCAA_into_peptide(each_ncaa_SDF_path, residue_id, input_pep_pdb_path, output_pep_pdb_path, NCAA_peptide_SDF_path, distance_info_txt_path)

    os.system(f'obabel -isdf {NCAA_peptide_SDF_path} -O {NCAA_peptide_PDB_path}')





