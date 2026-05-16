# Papers

This directory can be used as a local cache for graph trend filtering papers.
PDF files are intentionally ignored by git in this folder to avoid repository
bloat and accidental redistribution of third-party material.

Track citations, links, and notes here. Store downloaded PDFs locally if useful,
but do not commit them unless we deliberately decide to vendor them or use Git
LFS.

## Core Papers

### Trend Filtering on Graphs

Yu-Xiang Wang, James Sharpnack, Alex Smola, Ryan J. Tibshirani.
`Trend Filtering on Graphs`.

- arXiv: <https://arxiv.org/abs/1410.7690>
- JMLR PDF: <https://jmlr.csail.mit.edu/papers/volume17/15-147/15-147.pdf>
- Relevance: main graph trend filtering construction. Defines
  \(\ell_1\)-penalized graph difference operators and contrasts graph trend
  filtering with Laplacian smoothing and graph wavelet smoothing.

### The DFS Fused Lasso

Oscar Hernan Madrid Padilla, James Sharpnack, James G. Scott,
Ryan J. Tibshirani. `The DFS Fused Lasso: Linear-Time Denoising over General
Graphs`.

- JMLR page: <https://www.jmlr.org/beta/papers/v18/16-532.html>
- Relevance: graph fused lasso / graph total variation denoising, the
  \(k=0\) graph trend filtering case, with scalable graph-to-chain ideas.

### Fused \(\ell_1\) Trend Filtering on Graphs

Vladimir Pastukhov. `Fused \(\ell_1\) Trend Filtering on Graphs`.

- arXiv: <https://arxiv.org/abs/2401.05250>
- Relevance: newer work combining fusion regularizers and graph trend
  filtering, including computationally feasible iterative solvers.

## Local PDF Cache

Recommended local filenames if you download PDFs:

- `wang_sharpnack_smola_tibshirani_2016_trend_filtering_on_graphs.pdf`
- `madrid_padilla_sharpnack_scott_tibshirani_2018_dfs_fused_lasso.pdf`
- `pastukhov_2024_fused_l1_trend_filtering_on_graphs.pdf`

Because `*.pdf` is ignored here, downloading them will not automatically add
them to git.
