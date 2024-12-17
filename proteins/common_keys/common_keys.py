import glob
import os
import argparse
import time
import pandas as pd
from os.path import expanduser

parser = argparse.ArgumentParser(description='Commom Keys')
parser.add_argument('--path', '-path', metavar='path', default=os.path.join(expanduser('~'), 'Research', 'Protien_Database', 'extracted_new_samples', 'testing'),
                    help='Directory of input sample and other files.')

def common_keys(files):
    common_keys = []
    start = time.time()
    total_keys = {}
    # Flag to continue the code to calculate the total keys even when the common_keys is empty
    common_keys_empty = False
    for file in files:
        keys = {}
        with open(file, 'r') as f:
            for line in f:
                key, freq = line.split('\t')
                keys[key] = int(freq)
        total_keys[file] = keys
        if not common_keys_empty:
            if len(common_keys) <= 0:
                common_keys = list(keys.keys())
            else:
                common_keys = list(set(common_keys) & set(keys.keys()))
            # if at any point common keys is empty, set the flag to stop further common keys calculation
            if not common_keys:
                print("Common keys became empty, stopping further common keys calculation.")
                common_keys_empty = True
    
    print('Time taken for Common Keys Calculation: ', (time.time() - start) / 60)
    return total_keys, common_keys



if __name__ == "__main__":
    args = parser.parse_args()
    files = glob.glob(os.path.join(args.path, '*.keys_Freq_theta29_dist18*'))
    print(len(files))
    total_keys, common_keys_list = common_keys(files)
    print("Common Keys:", common_keys_list)

    data = {'fileName': [], 'total_keys': [], 'common_keys':[], 'sum_common_keys_freq': []}
    for file, keys in total_keys.items():
        common_count = len(common_keys_list)
        data['fileName'].append(os.path.splitext(os.path.basename(file))[0])
        data['total_keys'].append(len(keys.keys()))
        data['common_keys'].append(common_count)
        common_key_freq = sum(keys[key] for key in keys if key in common_keys_list)
        data['sum_common_keys_freq'].append(common_key_freq)

    current_directory = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_directory, "common_keys.csv")
    df = pd.DataFrame(data)
    df.to_csv(file_path, index=False)
