import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

def perform_q1_enrichment(gene_file, gmt_file, label):
    # 1. Load your extracted genes
    user_genes = set(open(gene_file).read().splitlines())
    n = len(user_genes)
    N = 20000  # Total Human Protein-Coding Genes

    # 2. Parse the GMT Database
    pathway_results = []
    with open(gmt_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            pathway_name = parts[0]
            pathway_genes = set(parts[2:])
            
            # Count Overlap (Hits)
            hits = user_genes.intersection(pathway_genes)
            k = len(hits)
            K = len(pathway_genes)
            
            if k >= 3: # Q1 Filter: Minimum 3 genes to be considered a 'pathway'
                # Hypergeometric Survival Function (1 - CDF)
                p_val = hypergeom.sf(k-1, N, K, n)
                pathway_results.append({
                    'Pathway': pathway_name,
                    'Hits': k,
                    'Path_Size': K,
                    'P_Value': p_val,
                    'Genes': ", ".join(list(hits))
                })

    # 3. Apply Benjamini-Hochberg (FDR) Correction
    df = pd.DataFrame(pathway_results)
    if df.empty:
        print(f"No significant pathways for {label}")
        return
        
    _, q_vals, _, _ = multipletests(df['P_Value'], method='fdr_bh')
    df['Q_Value (FDR)'] = q_vals
    
    # 4. Filter for Significance (Q < 0.05)
    df_sig = df[df['Q_Value (FDR)'] < 0.05].sort_values('Q_Value (FDR)')
    
    print(f"\n=== FINAL Q1 ENRICHMENT: {label} ===")
    print(df_sig[['Pathway', 'Hits', 'P_Value', 'Q_Value (FDR)']].head(10).to_string(index=False))
    df_sig.to_csv(f"{label}_full_results.csv", index=False)

# RUN THE REAL ANALYSIS
perform_q1_enrichment('strontium_genes.txt', 'c2.cp.reactome.v2023.2.Hs.symbols.gmt', 'STRONTIUM')
perform_q1_enrichment('calcium_genes.txt', 'c2.cp.reactome.v2023.2.Hs.symbols.gmt', 'CALCIUM')
