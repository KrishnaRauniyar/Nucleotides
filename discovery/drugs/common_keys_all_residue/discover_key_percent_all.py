import glob
import os
import argparse
import time
import pandas as pd
import csv

def get_keys_percents(files, chain):
    common_keys = {}
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                keyNo = line.split('\t')[0]
                if keyNo in common_keys:
                    common_keys[keyNo] += 1
                else:
                    common_keys[keyNo] = 1

    total_files = len(files)
    key_percentages = {key: round((count / total_files) * 100, 2) for key, count in common_keys.items()}
    df = pd.DataFrame(list(key_percentages.items()), columns=['key', f'{chain}%'])
    return df

## Arrange the filename as of the cross key filename
def rearrange_drug(drug):
    parts = drug.split('_')
    if len(parts) == 5:
        return f"{parts[0]}_{parts[1]}_{parts[3]}_{parts[2]}"
    return drug

## Extract from the cross key the protein residue name 
def extract_protein_type(protein_str):
    return protein_str.split('_')[2]

## create label based on NI and the others
def label_generator(row, filtered_cross_df):
    match = filtered_cross_df[filtered_cross_df['rearranged_drug'] == row['Protein']]
    if not match.empty:
        protein_type = extract_protein_type(match['protein'].values[0])
        return f'I_{protein_type}'
    return 'NI'

def create_label_to_proteins_dict(cross, drug, residue, result_dir):
    ## Reading the cross key file
    cross_df = pd.read_csv(cross)
    filtered_cross_df = cross_df[cross_df['drug'].str.split("_").str[2] == residue]
    filtered_cross_df.reset_index(drop=True, inplace=True)
    filtered_cross_df['rearranged_drug'] = filtered_cross_df['drug'].apply(rearrange_drug)
    
    # Reading the entire drug file
    drug_df = pd.read_csv(drug)
    filtered_drug_df = drug_df[drug_df['Protein'].str.split("_").str[-1] == residue]
    filtered_drug_df.reset_index(drop=True, inplace=True)
    
    # Generating labels
    filtered_drug_df['Label'] = filtered_drug_df.apply(label_generator, axis=1, args=(filtered_cross_df,))
    
    # Saving the labeled data to a CSV file
    label_drug_file_path = os.path.join(result_dir, f"label_drug_{residue}.csv")
    filtered_drug_df.to_csv(label_drug_file_path, index=None)

    ## genereating a dictionary with key as the labels and values as the drug list for common key calculations
    label_to_proteins = {}
    with open(label_drug_file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            label = row['Label']
            protein = row['Protein']
            # If the label is not already a key in the dictionary, add it with an empty list
            if label not in label_to_proteins:
                label_to_proteins[label] = []
            # Appending the protein to the list for this label
            label_to_proteins[label].append(protein)
    return label_to_proteins

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Common Keys calculation for discoveries')
    parser.add_argument('--path', '-p', help='Directory of the key generated.', required=True)
    parser.add_argument('--cross', '-c', help='CSV path to the cross keys', required=True)
    parser.add_argument('--drug', '-d', help='Path to the drug CSV file', required=True)
    parser.add_argument('--residue', '-r', help='Residue type', required=True)
    args = parser.parse_args()

    # Creating result directory if it does not exist
    result_dir = os.path.join(os.getcwd(), f'result_{args.residue}')
    os.makedirs(result_dir, exist_ok=True)

    label_to_proteins = create_label_to_proteins_dict(args.cross, args.drug, args.residue, result_dir)
    
    key_df = pd.DataFrame(columns=['key'])
    all_dfs = []
    for label, proteins in label_to_proteins.items():
        files = []
        for protein in proteins:
            search_pattern = os.path.join(args.path, f"{protein}*.keys_Freq_theta29_dist18*")
            files.extend(glob.glob(search_pattern))
        chain_df = get_keys_percents(files, label)
        all_dfs.append(chain_df)

        # Updating key_df with unique keys from current chain
        key_df = pd.merge(key_df, chain_df[['key']], how='outer')

    # Merging percentage columns for all chains
    final_df = key_df
    for df in all_dfs:
        final_df = pd.merge(final_df, df, on='key', how='left')

    # Filling NaN values with 0
    final_df.fillna(0, inplace=True) 
    file_path = os.path.join(result_dir, "common_keys_percent_one_vs_all.csv")
    final_df.to_csv(file_path, index=False)

