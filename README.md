# pooled

##### By: Paul Doan ([Google Scholar](https://scholar.google.com/citations?view_op=list_works&hl=en&hl=en&user=O4q49qQAAAAJ))


- A straight-forward Python workflow for wet-lab biologists to analyze pooled, functional genetic perturbation screens (e.g., CRISPR, RNAi, ORF) coupled with either viability or FACS-based readouts.

- The pre-processing and QC steps of this workflow and the data demo below are duplicated from the Broad GPP's `poola` Python package. Please check it out [here](https://pypi.org/project/poola/). 

- The main difference in this `pooled` workflow is in the p-value calculation step. `pooled` determines the p-value for each element (sgRNA, shRNA, or ORF barcode) by comparing its normalized log2-fold change (LFC) gene-level rank against a randomly-permuted null distribution of gene-level ranks, whereas `poola`'s Method 1 scales the LFC to a Gaussian distribution. This approach was previously implemented in `R` by Dr. Mikolaj Slabicki from Dr. Benjamin Ebert's lab @ Dana-Farber Cancer Institute ([Slabicki and Kozicka et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32494016/)). Hence, the demo data below is from the paper's _genome-wide FACS-based CRISPR-Cas9 screen_ in HEK293T cell line engineered to express a GFP/mCherry reporter of *CCNK* a.k.a. cyclin K protein's stabiblity (Fig.  2G).  Special thanks to Dr. Slabicki whom generously made the original code and data were generously made available in the Method section.

- Special thanks to Dr. John Doench, Peter DeWeirdt, and Dr. Mudra Hegde (co-authors of `poola`) and the team at GPP for their commitment to open-source science. More extensive wet- and dry-lab GPP resources are available [here](https://portals.broadinstitute.org/gpp/public/).

### Installation
```
pip install pooled
```

### Tutorial
To see a worked example and all functionalities, check it this interactive Jupyter notebook here

 [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/kiddo18/pooled/master?urlpath=https%3A%2F%2Fgithub.com%2Fkiddo18%2Fpooled%2Fblob%2Fmaster%2Fnotebook%2Fpooled_implementation-Ebert-Official.ipynb)
