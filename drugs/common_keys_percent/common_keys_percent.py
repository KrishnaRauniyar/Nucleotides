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
	df.to_csv(f'common_keys_percent_{chain}.csv', index = False)

	print(f'Time taken for Common Keys {chain} Percentage Calculation: ', (time.time() - start)/60)

if __name__ == '__main__':
	args = parser.parse_args()
	files = glob.glob(os.path.join(args.path, '*.keys_Freq_theta29_dist18*'))
	chains = ["G", "A", "C", "U", "DG", "DA", "DC", "DT"]
	for chain in chains:
		chain_files = [file for file in files if os.path.basename(file).endswith(f'_{chain}.keys_Freq_theta29_dist18')]
		get_keys_percents(chain_files, chain)




