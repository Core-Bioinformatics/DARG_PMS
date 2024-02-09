import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

cols = ['TSS1', 'TSS2', 'TSS3', 'RNA1', 'RNA2', 'RNA3', 'Distal1', 'Distal2', 'Distal3']

heatmaps_dir = '../../data/heatmaps'
output_dir = '../../figs/heatmaps'

for f in os.listdir(heatmaps_dir):
    if not f.endswith('.txt'):
        continue
    a = pd.read_csv(os.path.join(heatmaps_dir, f)).set_index('Genename')    
    cm = sns.clustermap(a[cols], cmap='coolwarm', col_cluster=False, metric='correlation')
    cm.ax_row_dendrogram.set_visible(False)
    fname = f.replace('.txt', '.pdf')
    plt.title(f)
    plt.savefig(os.path.join(output_dir, fname))
