import pandas as pd
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy.cluster import hierarchy
from sklearn.metrics import adjusted_rand_score
import argparse
import scipy.spatial as sp

print("Pandas version:", pd.__version__)
print("Seaborn version:", sns.__version__)
print("Matplotlib version:", matplotlib.__version__)
print("Scipy version:", scipy.__version__)

def load_and_process_data(csv_file):
    df = pd.read_csv(csv_file, header=None)
    df1 = df[0].str.split(';', expand=True)
    residue_split = df1[0].str.split('_', expand=True)
    residue_types = residue_split[3]
    df = df.drop(0, axis=1)
    merged_df = pd.concat([df1, df], axis=1)
    processed_df = merged_df.set_index(merged_df[0]).T.reset_index(drop=True).rename_axis(None, axis=1)
    processed_df = processed_df.drop(index=0).reset_index(drop=True)
    processed_df.index = processed_df.columns
    print("Processed DataFrame shape:", processed_df.shape)
    print("Residue Types:", residue_types.tolist())
    return processed_df, residue_types

def plot_clustermap(data, output_file):
    distance_matrix = 100 * data.astype(float) # in percentage

    # # Performing hierarchical clustering for rows and columns
    # row_linkage = hierarchy.linkage(data.astype(float).values, method='average')
    # col_linkage = hierarchy.linkage(data.astype(float).values.T, method='average')

    linkage = hierarchy.linkage(sp.distance.squareform(distance_matrix.values), method='average')

    # Creating a clustered heatmap
    plt.figure(figsize=(10, 8))
    clustered_heatmap = sns.clustermap(
        distance_matrix,
        row_cluster=True,
        col_cluster=True,
        row_linkage=linkage,
        col_linkage=linkage,
        cbar_kws={"shrink": 0.5},
        cbar_pos=(0.1, 0.83, 0.02, 0.18)
    )

    clustered_heatmap.savefig(output_file)
    plt.close()

    # Getting the order of rows after clustering
    row_order = clustered_heatmap.dendrogram_row.reordered_ind
    row_names = data.index[row_order]
    print("Row Order after clustering:", row_order)
    print("Row Names after clustering:", row_names.tolist())
    return row_names, row_order

def save_row_names(row_names, output_csv):
    # Writing the row names to a CSV file
    row_names_df = pd.DataFrame({'Drugs': row_names})
    row_names_df.to_csv(output_csv, index=False)

def calculate_ari(labels_true, labels_pred):
    ari = adjusted_rand_score(labels_true, labels_pred)
    return ari

def main(csv_file, plot_file, output_csv):
    data, residue_types = load_and_process_data(csv_file)
    row_names, row_order = plot_clustermap(data, plot_file)
    save_row_names(row_names, output_csv)

    # Aligning residue types with clustered data
    true_labels = residue_types.iloc[row_order].values
    print("True Labels:", true_labels)

    # Generating predicted labels from row order
    predicted_labels = row_order
    print("Predicted Labels:", predicted_labels)

    # Calculate ARI
    ari = calculate_ari(true_labels, predicted_labels)
    print(f'Adjusted Rand Index: {ari}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot heatmap clusters and calculate ARI")
    parser.add_argument('-p', '--input_path', type=str, required=True, help="Generalized CSV input path")
    args = parser.parse_args()
    main(args.input_path, 'clustermap.png', 'clustermap.csv')


