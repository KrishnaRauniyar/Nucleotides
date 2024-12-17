import glob
import os
import argparse
import time
import pandas as pd
from os.path import expanduser

parser = argparse.ArgumentParser(description='Commom Keys')
parser.add_argument('--path', '-path', metavar='path', default=os.path.join(expanduser('~'), 'Research', 'Protien_Database', 'extracted_new_samples', 'testing'),
                    help='Directory of input sample and other files.')

def find_files(path, patterns):
    matched_files = []
    for pattern in patterns:
        files = glob.glob(os.path.join(path, pattern))
        matched_files.extend(files)
    return matched_files

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
            common_keys = list(keys.keys())
    print('Time taken for Common Keys Calculation: ', (time.time() - start) / 60)
    return total_keys, common_keys

if __name__ == "__main__":
    args = parser.parse_args()
    # files = glob.glob(os.path.join(args.path, '*.keys_Freq_theta29_dist18*'))
    # patterns = [
    # '*_A.keys_Freq_theta29_dist18*',
    # '*_G.keys_Freq_theta29_dist18*',
    # '*_DA.keys_Freq_theta29_dist18*',
    # '*_DG.keys_Freq_theta29_dist18*']

    patterns = [
    '*_C.keys_Freq_theta29_dist18*',
    '*_DC.keys_Freq_theta29_dist18*',
    '*_DT.keys_Freq_theta29_dist18*',
    '*_U.keys_Freq_theta29_dist18*']

    path = args.path 

    matched_files = find_files(path, patterns)

    print(len(matched_files))
    total_keys, common_keys_list = common_keys(matched_files)
    print("Common Keys:", common_keys_list)

    data = {'fileName': [], 'total_keys': [], 'common_keys':[], 'sum_common_keys_freq': []}
    for file, keys in total_keys.items():
        common_count = len(common_keys_list)
        data['fileName'].append(os.path.splitext(os.path.basename(file))[0])
        data['total_keys'].append(len(keys.keys()))
        data['common_keys'].append(common_count)
        common_key_freq = sum(keys[key] for key in keys if key in common_keys_list)
        data['sum_common_keys_freq'].append(common_key_freq)

    df = pd.DataFrame(data)
    # df.to_csv('common_keys_pu.csv')
    df.to_csv('common_keys_py.csv')

