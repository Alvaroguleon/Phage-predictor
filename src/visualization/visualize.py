import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import community

def plot_network(df, communities, path_to_save='reports/figures/phage_network.png'):
    print("-- Creating cluster map with NetworkX...")

    # Getting color for the nodes based on majority 'staining' in the community
    staining_map = df[["Accession", "staining"]].set_index('Accession')['staining'].to_dict()

    # Determine the majority class for each community
    community_staining = {}
    for community_number, community in enumerate(communities, start=1):
        staining_counts = {'positive': 0, 'negative': 0}
        for node in community:
            if node in staining_map:  # Count only if staining is known
                staining_counts[staining_map[node]] += 1
        majority_staining = max(staining_counts, key=staining_counts.get, default='negative')
        community_staining[community_number] = majority_staining

    # Calculate positions for all nodes to be used in plotting
    pos = nx.spring_layout(G, seed=42)  # Using a fixed seed for reproducibility

    # Setup figure
    plt.figure(figsize=(20, 20))

    # Draw all edges first to avoid redrawing
    nx.draw_networkx_edges(G, pos, alpha=0.1, edge_color='gray')

    # Draw nodes in batches based on community staining
    for community_index, community_nodes in enumerate(communities, start=1):
        color = 'blue' if community_staining[community_index] == 'positive' else 'red'
        nx.draw_networkx_nodes(G, pos, nodelist=list(community_nodes), node_color=color, node_size=50, alpha=0.5)

    # Draw community labels
    for community_index, community_nodes in enumerate(communities, start=1):
        representative_node = next(iter(community_nodes))
        plt.text(pos[representative_node][0], pos[representative_node][1], str(community_index),
                 horizontalalignment='center', verticalalignment='center',
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='black', boxstyle='round,pad=0.5'))

    plt.axis('off')
    plt.title('Network of Genomes by Gram Staining Classification')

    # Save the figure
    plt.savefig(path_to_save, dpi=300, bbox_inches='tight')  # Saves the figure to a file
    plt.close()  # Close the figure to free up memory
