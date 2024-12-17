import glob, os, argparse, time
import argparse
from os.path import expanduser
import pandas as pd


parser = argparse.ArgumentParser(description='Commom Keys')
parser.add_argument('--path', '-path', metavar='path', default=os.path.join(expanduser('~'),'Research', 'Protien_Database','extracted_new_samples', 'testing'),
                    help='Directory of input sample and other files.')

def common_keys(files):
    common_keys = []
    start = time.time()
    total_keys = {}
    for file in files:
        keys = {}
        with open(file, 'r') as f:
            for line in f:
                key, freq = line.split('\t')
                keys[key] = int(freq)
        total_keys[file] = keys
        if common_keys:
            common_keys = list(set(common_keys) & set(keys.keys()))
        else:
            common_keys = list(set(keys.keys()))
    print('Time taken for Common Keys Calculation: ', (time.time() - start)/60)
    return total_keys, common_keys

if __name__ == "__main__":
    args = parser.parse_args()
    files = glob.glob(os.path.join(args.path,'*.keys_Freq_theta29_dist18*'))
    print(len(files))

    chains = ["G", "A", "C", "U", "DG", "DA", "DC", "DT"]
    for chain in chains:
        common_keys_list = []
        total_keys = {}
        
        chain_files = [file for file in files if os.path.basename(file).endswith(f'_{chain}.keys_Freq_theta29_dist18')]
        total_keys, common_keys_list = common_keys(chain_files)
        print(f"Common Keys of {chain}:", len(common_keys_list))		

        data = {'fileName': [], 'total_keys': [], 'common_keys': [], 'sum_common_keys_freq': []}
        for file, keys in total_keys.items():
            common_count = len(common_keys_list)
            data['fileName'].append(os.path.splitext(os.path.basename(file))[0])
            data['total_keys'].append(len(keys.keys()))
            data['common_keys'].append(common_count)
            common_keys_freq = sum([keys[key] for key in keys if key in common_keys_list])
            data['sum_common_keys_freq'].append(common_keys_freq)
            
        df = pd.DataFrame(data)
        df.to_csv(f'common_keys_{chain}.csv')



