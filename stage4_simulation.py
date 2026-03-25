import networkx as nx
import random
import matplotlib.pyplot as plt

EXPERIMENTAL_CI = 0.3455
GENE_COUNT = 55
TRIALS = 1000

print("Loading HIPPIE network...")

G_full = nx.Graph()

with open("hippie_current.txt") as file:
    for line in file:
        parts = line.strip().upper().split("\t")

        try:
            score = float(parts[4])
        except:
            continue

        if score <= 0.73:
            continue

        p1 = parts[0].replace("_HUMAN","")
        p2 = parts[2].replace("_HUMAN","")

        G_full.add_edge(p1,p2)

all_pool = list(G_full.nodes())

print("HIPPIE nodes:",len(all_pool))

ci_values = []
better_than_strontium = 0

print("Starting Monte Carlo simulation...")

for i in range(TRIALS):

    random_genes = random.sample(all_pool,GENE_COUNT)

    G_rand = G_full.subgraph(random_genes)

    if G_rand.number_of_nodes() > 0:

        lcc = len(max(nx.connected_components(G_rand),key=len))
        ci_rand = lcc / GENE_COUNT

        ci_values.append(ci_rand)

        if ci_rand >= EXPERIMENTAL_CI:
            better_than_strontium += 1

    if (i+1) % 100 == 0:
        print("Completed",i+1,"trials")

p_value = better_than_strontium / TRIALS

print("P-value:",p_value)

print("Total CI values:", len(ci_values))

import matplotlib.pyplot as plt

plt.figure(figsize=(8,5))

plt.hist(ci_values, bins=10, edgecolor="black", linewidth=1.5)

plt.axvline(EXPERIMENTAL_CI, linestyle='--', linewidth=3)

plt.xlim(0,0.36)

plt.xlabel("Connectivity Index (CI)")
plt.ylabel("Frequency")
plt.title("Monte Carlo CI Distribution")

plt.tight_layout()

plt.savefig("monte_carlo_distribution.png", dpi=300)

print("Plot saved as monte_carlo_distribution.png")
