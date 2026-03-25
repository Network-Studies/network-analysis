import networkx as nx
import re

print("1. Cleaning the AmiGO data...")
gs_genes = set()
with open('calcium_genes.txt', 'r') as file:
    for line in file:
        gene = line.strip().upper()
        if not gene: continue
        # Strip any existing artifacts from AmiGO
        gene = re.sub(r'_HUMAN.*$', '', gene)
        gs_genes.add(gene)

print(f"   Found {len(gs_genes)} unique Gold Standard genes.")

print("2. Scanning HIPPIE database line by line... (this takes ~5 seconds)")
G = nx.Graph()

with open('hippie_current.txt', 'r') as file:
    for i, line in enumerate(file):
        parts = line.strip().upper().split('\t')
        
        # Print a sample of the first line so we can see the exact formatting
        if i == 0:
            print(f"   [DEBUG] First row data: {parts[:6]}")

        # Find the confidence score (the only float between 0 and 1 in the row)
        score = 0.0
        for p in parts:
            try:
                val = float(p)
                if 0.0 < val <= 1.0:
                    score = val
                    break
            except ValueError:
                pass
                
        # Immediately filter out low-confidence interactions
        if score <= 0.73:
            continue
            
        # Hunt for our genes anywhere in this row
        row_genes = set()
        for p in parts:
            # Strip the _HUMAN tag HIPPIE usually attaches
            clean_p = p.replace('_HUMAN', '')
            
            if clean_p in gs_genes:
                row_genes.add(clean_p)
                
        # If exactly 2 of our Gold Standard genes interact in this row, add the edge!
        row_list = list(row_genes)
        if len(row_list) == 2:
            G.add_edge(row_list[0], row_list[1])

print(f"\n3. Network built with {G.number_of_nodes()} connected nodes and {G.number_of_edges()} edges.")

if G.number_of_nodes() > 0:
    connected_components = list(nx.connected_components(G))
    lcc_nodes = max(connected_components, key=len)
    lcc_size = len(lcc_nodes)
    ci_gs = lcc_size / len(gs_genes)
    
    print("\n=== GOLD STANDARD NETWORK METRICS ===")
    print(f"LCC Size: {lcc_size}")
    print(f"Connectivity Index (CI_GS): {ci_gs:.4f}")
else:
    print("Error: Still 0 edges.")
