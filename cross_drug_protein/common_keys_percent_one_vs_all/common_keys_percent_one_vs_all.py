import glob
import os
import argparse
import time
import pandas as pd
from os.path import expanduser

parser = argparse.ArgumentParser(description='Common Keys')
parser.add_argument('--path', '-path', metavar='path', default=os.path.join(expanduser('~'), 'Research', 'Protien_Database', 'extracted_new_samples', 'testing'),
                    help='Directory of input sample and other files.')

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
    key_percentages = {key: round((count / total_files) * 100, 4) for key, count in common_keys.items()}
    df = pd.DataFrame(list(key_percentages.items()), columns=['key', f'{chain}_key%'])
    return df

if __name__ == '__main__':
    start = time.time()
    args = parser.parse_args()
    files = glob.glob(os.path.join(args.path, '*.keys_Freq_theta29_dist18*'))
    
    chains = set()
    for file in files:
        filename = os.path.basename(file)
        chain_part = filename.split('.')[0].split('_')[-1]
        chains.add(chain_part)
    chains = sorted(chains)
    print(f"The chains are: {chains}")

    key_df = pd.DataFrame(columns=['key'])
    all_dfs = []
    for chain in chains:
        chain_files = [file for file in files if os.path.basename(file).endswith(f'_{chain}.keys_Freq_theta29_dist18')]
        chain_df = get_keys_percents(chain_files, chain)
        all_dfs.append(chain_df)
        
        # Updating key_df with unique keys from current chain
        key_df = pd.merge(key_df, chain_df[['key']], how='outer')

    # Merging percentage columns for all chains
    final_df = key_df
    for df in all_dfs:
        final_df = pd.merge(final_df, df, on='key', how='left')

    # Filling NaN values with 0
    final_df.fillna(0, inplace=True)

    current_directory = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_directory, "common_keys_percent_one_vs_all.csv")
    final_df.to_csv(file_path, index=False)
    print(f'Time taken for Calculation: ', (time.time() - start) / 60)

