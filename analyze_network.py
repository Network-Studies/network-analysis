import pandas as pd
import networkx as nx
from itertools import combinations

def run_analysis(filename):
    # 1. Load data and filter for significance
    df = pd.read_csv(filename)
    sig_df = df[df['Q_Value (FDR)'] < 0.05]
    
    # 2. Build the Network
    G = nx.Graph()
    for _, row in sig_df.iterrows():
        genes = [g.strip() for g in str(row['Genes']).split(',')]
        if len(genes) > 1:
            G.add_edges_from(combinations(genes, 2))
            
    # 3. Calculate the Metrics
    eff = nx.global_efficiency(G)
    dens = nx.density(G)
    assort = nx.degree_assortativity_coefficient(G)
    
    degs = [d for n, d in G.degree()]
    n = G.number_of_nodes()
    centralization = sum(max(degs) - d for d in degs) / ((n-1)*(n-2)) if n > 2 else 0

    print(f"--- Results for {filename} ---")
    print(f"Total Genes (Nodes): {n}")
    print(f"Global Efficiency: {eff:.3f}")
    print(f"Degree Centralization: {centralization:.3f}")
    print(f"Assortativity: {assort:.3f}")
    print(f"Network Density: {dens:.3f}\n")

# Run it on your two files
run_analysis('STRONTIUM_full_results.csv')
run_analysis('CALCIUM_full_results.csv')
