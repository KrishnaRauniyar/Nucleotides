import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def boxoverlay(input):
    groups = ['A', 'DA', 'G', 'DG', 'C', 'DC', 'U', 'DT']
    data = pd.DataFrame()

    for group in groups:
        filename = f"{input}/{group}_list.txt" 
        scores = pd.read_csv(filename, header=None, names=['Similarity'])
        scores['Group'] = group
        data = pd.concat([data, scores])

    # Setting the aesthetics for the plots
    sns.set(style="whitegrid")

    # Defining the palette
    palette = sns.color_palette("colorblind", len(groups))

    # Creating the box overlay plot
    plt.figure(figsize=(10, 6))
    ax = sns.stripplot(x='Group', y='Similarity', data=data, jitter=True, hue='Group', palette=palette, alpha=0.5, legend=False)
    sns.boxplot(x='Group', y='Similarity', data=data, whis=1.5, width=0.5, palette=palette, fliersize=0, ax=ax)

    # Annotating the mean values
    means = data.groupby('Group')['Similarity'].mean()
    for group in groups:
        mean = means[group]
        ax.text(groups.index(group), mean, f'{mean:.1f}', color='black', ha="center")

    # Setting labels and title
    ax.set_ylabel('Structural Similarity (%)')
    ax.set_title('Structural Similarity Scores by Group')
    plt.savefig('box_plot')


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Box Overlay Plot")
    parser.add_argument('-p','--input_path', type=str, required=True, help="CSV input path")
    args = parser.parse_args()
    boxoverlay(args.input_path)