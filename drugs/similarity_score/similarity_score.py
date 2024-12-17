import pandas as pd
import numpy as np
import argparse

def similarity_score(file_path, output_path):
    data = pd.read_csv(file_path, header=None, delimiter=';')
    keys = data[0]
    values = data[1].str.split(',', expand=True).astype(float)
    df = pd.DataFrame(values.values, index=keys, columns=keys)

    # Calculating the percentage difference from 100
    df_diff = 100 - (df * 100)

    residues = ['U', 'DT', 'A', 'DA', 'G', 'DG', 'C', 'DC']

    for residue in residues:
        # Filtering the DataFrame for the current residue
        filtered_df = df_diff[df_diff.columns[df_diff.columns.str.endswith(f'_{residue}')]]
        filtered_df = filtered_df.loc[filtered_df.index.str.endswith(f'_{residue}')]
        
        # Getting the upper triangle of the filtered DataFrame
        upper_triangle = np.triu(filtered_df, k=1)
        
        # Flattening and round the upper triangle, removing 0 values
        upper_triangle_list = upper_triangle.flatten().tolist()
        filtered_list = [round(val, 2) for val in upper_triangle_list if val != 0]
        
        # Saving the filtered list to a text file
        with open(f'{output_path}/{residue}_list.txt', 'w') as f:
            for item in filtered_list:
                f.write(f"{item}\n")

def main(input_path, output_path):
    similarity_score(input_path, output_path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculating the similarity score")
    parser.add_argument('-p', '--input_path', type=str, required=True, help="Generalized CSV input path")
    parser.add_argument('-o', '--output_path', type=str, required=True, help="results directory")
    args = parser.parse_args()
    main(args.input_path, args.output_path)
