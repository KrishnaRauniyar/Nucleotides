import glob
import os
import argparse
import time
import pandas as pd
from os.path import expanduser

parser = argparse.ArgumentParser(description='Commom Keys')
parser.add_argument('--path', '-path', metavar='path', default=os.path.join(expanduser('~'), 'Research', 'Protien_Database', 'extracted_new_samples', 'testing'),
                    help='Directory of input sample and other files.')

def get_keys_percents(files, chain):
	common_keys = {}
	start = time.time()
	for file in files:
		for line in open(file, 'r'):
			keyNo = line.split('\t')[0]
			if keyNo in common_keys:
				common_keys[keyNo] += 1
			else:
				common_keys[keyNo] = 1

	df = pd.DataFrame(common_keys.items(), columns = ['key', 'key_occurance'])
	df['total_files'] = len(files)
	df['key_percent'] = (100*df['key_occurance'])/(len(files))	

	current_directory = os.path.dirname(os.path.abspath(__file__))
	file_path = os.path.join(current_directory, f"common_keys_percent_{chain}.csv")
	df.to_csv(file_path, index=False)

	print(f'Time taken for Common Keys {chain} Percentage Calculation: ', (time.time() - start)/60)

if __name__ == '__main__':
	args = parser.parse_args()
	files = glob.glob(os.path.join(args.path, '*.keys_Freq_theta29_dist18*'))

	chains = set()
	for file in files:
		filename = os.path.basename(file)
		chain_part = filename.split('.')[0].split('_')[-1]
		chains.add(chain_part)
	chains = sorted(chains)
	print(f"The chains are: {chains}")

	for chain in chains:
		chain_files = [file for file in files if os.path.basename(file).endswith(f'_{chain}.keys_Freq_theta29_dist18')]
		get_keys_percents(chain_files, chain)




