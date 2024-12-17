import pandas as pd
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.spatial as sp
import scipy.cluster.hierarchy as hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.metrics import adjusted_rand_score
import argparse

def process_data(input_path):
    df = pd.read_csv(input_path)
    headers = df.columns.to_list()
    headers.pop(0)
    df.insert(0, 'Protein', headers)
    df = df.drop('Unnamed: 0', axis=1)
    data = df.set_index('Protein')
    return data

def plot_clustermap(data, output_file):
    distance_matrix = 100 * data.astype(float) # in percentage
    linkage = hierarchy.linkage(sp.distance.squareform(distance_matrix.values), method='average')

    # Creating a clustered heatmap
    plt.figure(figsize=(10, 8))
    clustered_heatmap = sns.clustermap(
        distance_matrix,
        row_cluster=True,
        col_cluster=True,
        row_linkage=linkage,
        col_linkage=linkage,
        cbar_kws={"shrink": 0.5}
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
    data = process_data(csv_file)
    row_names, row_order = plot_clustermap(data, plot_file)
    save_row_names(row_names, output_csv)

    # # Aligning residue types with clustered data
    # true_labels = residue_types.iloc[row_order].values
    # print("True Labels:", true_labels)

    # # Generating predicted labels from row order
    # predicted_labels = row_order
    # print("Predicted Labels:", predicted_labels)

    # # Calculate ARI
    # ari = calculate_ari(true_labels, predicted_labels)
    # print(f'Adjusted Rand Index: {ari}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plot heatmap clusters and calculate ARI")
    parser.add_argument('-p', '--input_path', type=str, required=True, help="Generalized CSV input path")
    args = parser.parse_args()
    main(args.input_path, 'clustermap.png', 'clustermap.csv')