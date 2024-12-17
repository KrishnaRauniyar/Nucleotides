import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
import argparse

def mds(filePath):
    # Reading the CSV file
    df = pd.read_csv(filePath, header=None, sep=';')

    # Parsing the data
    labels = []
    distances = []

    for _, row in df.iterrows():
        label, dist = row[0], row[1]
        labels.append(label.split("_")[-1]) # Extract the class label
        distances.append(list(map(float, dist.split(","))))

    # Create the distance matrix
    distance_matrix = np.zeros((len(distances), len(distances)))

    for i in range(len(distances)):
        for j in range(len(distances[i])):
            distance_matrix[i][j] = distances[i][j]

    # Applying MDS
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
    X_mds = mds.fit_transform(distance_matrix)

    # Unique class labels
    unique_labels = list(set(labels))
    colors = plt.cm.get_cmap('tab10', len(unique_labels))

    # Plotting the results
    plt.figure(figsize=(10, 6))

    for i, label in enumerate(unique_labels):
        idx = [j for j, lbl in enumerate(labels) if lbl == label]
        plt.scatter(X_mds[idx, 0], X_mds[idx, 1], color=colors(i), label=label, s=50)

    plt.xlabel('MDS1')
    plt.ylabel('MDS2')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.title('MDS Plot')
    plt.savefig('mds_plot.png', bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser("MDS")
    parser.add_argument('-p','--input_path', type=str, required=True, help="CSV input path")
    args = parser.parse_args()
    mds(args.input_path)
