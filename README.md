# mtb-pocket-analyzer
# M. tuberculosis Binding Pocket Analyser

A computational pipeline for identifying selective drug targets in *Mycobacterium tuberculosis* by comparing predicted protein structures against a human control protein.

Built as an independent portfolio project during my BSc in Biological Sciences at IISER Bhopal, alongside a research internship in organometallic chemistry at NISER Bhubaneswar (2026).

---

## Why I built this

Tuberculosis kills around 1.5 million people every year. The bacterium that causes it has become resistant to many existing drugs, which means we urgently need new ones. Finding a new drug starts with finding a good target — a bacterial protein you can disable without harming the patient.

The challenge is selectivity. A drug that sticks to a bacterial protein might also stick to a similar human protein, causing side effects or toxicity. So the question I tried to answer computationally was:

> *Which essential M. tuberculosis proteins have binding pockets shaped differently enough from human proteins that selective drug design is feasible?*

To answer that, I built a four-phase pipeline — structure download, quality filtering, pairwise structural alignment, and binding pocket analysis — ending in a ranked list of the most promising drug targets.

---

## Results

**Top 3 targets from the full analysis:**

| Rank | Gene | Pathway | Druggability | TM-score vs human | Selectivity index |
|------|------|---------|-------------|-------------------|------------------|
| #1 | kasB | FAS-II cell wall | 0.996 | 0.286 | 0.711 |
| #2 | fabH | FAS-II cell wall | 0.806 | 0.323 | 0.546 |
| #3 | fabD | FAS-II cell wall | 0.564 | 0.252 | 0.422 |

kasB came out on top — near-perfect druggability and a binding pocket structurally distant from the human control. All three top targets are from the FAS-II fatty acid synthesis pathway, which builds the mycolic acid layer unique to the M. tuberculosis cell wall. Humans don't have this pathway, which is exactly why it's such an attractive drug target space.

**Structural similarity heatmap (Phase 2):**

![TM-score heatmap](tm_score_heatmap%20(1).png)

The human DHFR control (red label) clusters entirely separately from all M. tuberculosis proteins in the dendrogram — confirming these bacterial enzymes have fundamentally different folds from the human protein. Also worth noting: hadA and hadB scored 0.60 with each other, the highest off-diagonal score in the matrix. They are the two subunits of the same enzyme, so the pipeline picking that up automatically is a nice sanity check.

**Pocket similarity clustering (Phase 3):**

![UMAP plot](pocket_umap%20(3).png)

DHFR sits isolated in pocket space (bottom right), far from most M. tuberculosis proteins. One exception: fabH clusters near DHFR in pocket space, suggesting its pocket chemistry is more human-like despite ranking #2 overall. This is a selectivity caveat worth noting — pocket-level analysis reveals risks that fold-level comparison alone can miss.

---

## Pipeline overview

```
UniProt IDs
     ↓
Phase 1 — Fetch AlphaFold predicted structures
          Filter low-confidence residues (pLDDT ≥ 70)
     ↓
Phase 2 — TM-align pairwise alignment (91 pairs)
          14×14 similarity matrix → clustered heatmap
     ↓
Phase 3 — fpocket binding site detection
          UMAP pocket clustering → selectivity ranking
     ↓
Phase 4 — Streamlit web app for interactive single-protein analysis
```

Each phase is a self-contained Google Colab notebook. They run in order — each reads the previous phase's output from Google Drive.

---

## Notebooks

| Notebook | What it does | Open in Colab |
|----------|-------------|---------------|
| `AlphaFold_Phase1_Colab.ipynb` | Downloads and cleans 14 M.tb + 1 human structure | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com) |
| `AlphaFold_Phase2_Colab.ipynb` | TM-align pairwise alignment, heatmap | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com) |
| `AlphaFold_Phase3_Colab.ipynb` | fpocket detection, UMAP, target ranking | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com) |

No local installation needed — everything runs on free Colab.

---

## Web app

The Streamlit app lets you analyse any protein by UniProt ID in the browser.

**Live demo:** [add your Streamlit URL here after deploying]

Enter any UniProt accession (e.g. `P9WJA1` for kasB) and the app:
- Fetches the AlphaFold predicted structure automatically
- Shows the per-residue pLDDT confidence profile
- Runs fpocket and returns druggability score, pocket volume, hydrophobicity, polarity
- Lists all detected pockets ranked by druggability

**To run locally:**
```bash
git clone https://github.com/yourusername/mtb-pocket-analyser
cd mtb-pocket-analyser
pip install -r requirements_streamlit.txt
streamlit run app.py
```

---

## Methods

**Structures:** AlphaFold Database v4 predicted structures fetched via REST API for 13 essential *M. tuberculosis* H37Rv enzymes (selected from the DEG database, Sassetti et al. 2003) and human DHFR (P00374) as a structural divergence control.

**Confidence filtering:** Per-residue pLDDT scores are stored in the B-factor column of AlphaFold PDB files. Low-confidence terminal regions were trimmed using a sliding window approach (window = 5, threshold = 70) before alignment to prevent disordered regions from distorting structural comparison scores.

**Structural alignment:** TM-align (Zhang & Skolnick, 2005) was run on all 91 pairwise combinations of the 14 filtered structures. TM-scores were averaged across both normalisation directions. The resulting 14×14 matrix was visualised as a clustered heatmap (seaborn clustermap, Ward linkage).

**Pocket detection:** fpocket v3 (Le Guilloux et al., 2009) was run on each filtered structure. The top-ranked pocket per protein (by druggability score) was selected as the representative binding site. Four features were extracted: druggability score, volume (Å³), hydrophobicity score, and polarity score. UMAP (McInnes et al., 2018) was applied to the standardised 4-feature matrix for 2D visualisation.

**Selectivity index:** Druggability score × (1 − TM-score vs DHFR). Rewards proteins with high pocket druggability and low structural similarity to the human control.

---

## Proteins studied

| UniProt | Gene | Pathway | Description |
|---------|------|---------|-------------|
| P9WIA3 | inhA | FAS-II | Enoyl-ACP reductase — primary isoniazid target |
| P9WJA3 | kasA | FAS-II | Beta-ketoacyl-ACP synthase I |
| P9WJA1 | kasB | FAS-II | Beta-ketoacyl-ACP synthase II |
| P9WIS3 | mabA | FAS-II | 3-ketoacyl reductase |
| P9WI77 | hadA | FAS-II | Hydroxyacyl-ACP dehydratase subunit A |
| P9WI75 | hadB | FAS-II | Hydroxyacyl-ACP dehydratase subunit B |
| P9WKG3 | fabH | FAS-II | 3-oxoacyl-ACP synthase III |
| P9WKR5 | fabD | FAS-II | Malonyl-CoA ACP transacylase |
| P9WNE9 | gyrA | DNA-rep | DNA gyrase subunit A |
| P9WNF1 | gyrB | DNA-rep | DNA gyrase subunit B |
| P0A5R3 | dnaN | DNA-rep | DNA polymerase III beta clamp |
| P9WPZ9 | echA6 | Lipid | Enoyl-CoA hydratase |
| P9WK71 | accA3 | Lipid | Acyl-CoA carboxylase alpha |
| P00374 | DHFR | Human | Human dihydrofolate reductase (control) |

---

## Repository structure

```
mtb-pocket-analyser/
├── app.py                          ← Streamlit web app
├── requirements_streamlit.txt      ← App dependencies
├── AlphaFold_Phase1_Colab.ipynb   ← Data acquisition
├── AlphaFold_Phase2_Colab.ipynb   ← Structural alignment
├── AlphaFold_Phase3_Colab.ipynb   ← Pocket analysis
└── results/
    ├── tm_score_heatmap.png        ← Phase 2 figure
    ├── pocket_umap.png             ← Phase 3 figure
    ├── pocket_features.png         ← Phase 3 bar charts
    ├── tm_score_matrix.tsv         ← 14×14 TM-score matrix
    ├── pairwise_alignments.tsv     ← All 91 alignment results
    ├── pocket_features.tsv         ← Pocket metrics per protein
    └── final_target_ranking.tsv    ← Ranked drug target list
```

---

## References

- Jumper et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596, 583–589.
- Zhang & Skolnick (2005). TM-align: a protein structure alignment algorithm based on the TM-score. *Nucleic Acids Research*, 33(7), 2302–2309.
- Le Guilloux et al. (2009). Fpocket: an open source platform for ligand pocket detection. *BMC Bioinformatics*, 10, 168.
- Sassetti et al. (2003). Genes required for mycobacterial growth defined by high density mutagenesis. *Molecular Microbiology*, 48(1), 77–84.
- McInnes et al. (2018). UMAP: Uniform Manifold Approximation and Projection for dimension reduction. *arXiv*, 1802.03426.



*BS-MS Biological Sciences major with Chemistry minor — IISER Bhopal*
*Independent computational biology project, 2026*
