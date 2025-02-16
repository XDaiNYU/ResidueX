import os
import os
import pandas as pd
import re

def pair_files_in_directory(directory_path):
    # Extract the directory name
    directory_name = os.path.basename(directory_path)

    # List all files in the directory
    all_files = [f for f in os.listdir(directory_path) if os.path.isfile(os.path.join(directory_path, f))]

    # Process the files to extract number and type
    file_pairs = {}
    for file in all_files:
        match = re.search(r'(\d+)_sp_(pep|pro)\.pdb', file)
        if match:
            number = int(match.group(1))  # Extract number
            file_type = match.group(2)    # Extract file type (pep/pro)

            # Group files by their number
            if number not in file_pairs:
                file_pairs[number] = {'pep': None, 'pro': None}
            file_pairs[number][file_type] = file

    # Create pairs and add them to a list
    data = []
    for number, files in file_pairs.items():
        if files['pep'] and files['pro']:  # Ensure both pep and pro files are present
            data.append([directory_name, number, files['pep'], files['pro']])

    # Create DataFrame
    df = pd.DataFrame(data, columns=['Directory', 'Identifier', 'Pep_File', 'Pro_File'])
    return df


def ensure_directory_exists(directory_path):
    """Creates a directory if it doesn't exist."""
    os.makedirs(directory_path, exist_ok=True)

