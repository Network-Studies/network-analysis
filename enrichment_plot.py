import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = {
"Pathway":[
"RTK Signaling",
"IL-4 & IL-13 Signaling",
"Signaling by Interleukins",
"Cytokine Signaling",
"Growth Factor Disorders",
"PI3K/AKT Regulation",
"IL-10 Signaling",
"RNA Pol II Transcription",
"VEGF Signaling",
"ECM Organization",
"MET Signaling",
"Hemostasis",
"MAPK Signaling",
"Estrogen Signaling"
],

"Sr_hits":[26,17,22,24,18,11,8,23,9,11,6,15,11,7],

"Sr_FDR":[
1.21e-24,
6.06e-24,
3.85e-20,
5.17e-18,
1.72e-15,
3.12e-13,
8.69e-12,
8.69e-12,
1.86e-10,
4.50e-09,
2.08e-09,
2.75e-09,
1.06e-08,
1.11e-08
],

"Ca_hits":[25,17,21,24,19,12,4,24,8,9,4,17,13,8],

"Ca_FDR":[
1.73e-21,
2.25e-22,
3.84e-17,
3.12e-16,
1.71e-15,
3.54e-14,
3.21e-05,
3.01e-11,
1.76e-08,
2.07e-06,
7.88e-06,
3.65e-10,
3.98e-10,
1.20e-09
]
}

df = pd.DataFrame(data)

df["Sr_score"] = -np.log10(df["Sr_FDR"])
df["Ca_score"] = -np.log10(df["Ca_FDR"])

y = np.arange(len(df))

plt.figure(figsize=(10,7))

plt.scatter(
df["Sr_score"],
y,
s=df["Sr_hits"]*25,
label="Strontium"
)

plt.scatter(
df["Ca_score"],
y,
s=df["Ca_hits"]*25,
marker="s",
label="Calcium"
)

plt.yticks(y,df["Pathway"])

plt.xlabel("-log10(FDR)")
plt.title("Reactome Pathway Enrichment Comparison (Sr vs Ca)")

plt.legend()

plt.subplots_adjust(left=0.40)

plt.savefig(
"reactome_comparison_final.png",
dpi=300,
bbox_inches="tight"
)

plt.show()
