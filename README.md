# pooled

<img src="https://github.com/kiddo18/pooled/assets/43038912/da233ba4-ede4-409c-8df2-ba47d2cab7fb" width="300" height="300">

##### By: Paul Doan ([Google Scholar](https://scholar.google.com/citations?view_op=list_works&hl=en&hl=en&user=O4q49qQAAAAJ))


- A straight-forward Python package for wet-lab biologists to analyze pooled, functional genetic perturbation screens (e.g., CRISPR, RNAi, ORF) coupled with either viability or FACS-based readouts.

- The pre-processing and QC steps of this workflow below are duplicated from the Broad GPP's `poola` Python package. Please check it out [here](https://pypi.org/project/poola/). 

- The main difference in this `pooled` workflow is in the p-value calculation step. `pooled` determines the p-value for each element (sgRNA, shRNA, or ORF barcode) by comparing its normalized log2-fold change (LFC) gene-level rank against a randomly-permuted null distribution of gene-level ranks, whereas `poola`'s Method 1 scales the LFC to a Gaussian distribution. This approach was previously implemented in `R` by Dr. Mikolaj Slabicki from Dr. Benjamin Ebert's lab @ Dana-Farber Cancer Institute ([Slabicki and Kozicka et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32494016/)). Hence, the demo data below is from the paper's _genome-wide FACS-based CRISPR-Cas9 screen_ in HEK293T cell line engineered to express a GFP/mCherry reporter of *CCNK* a.k.a. cyclin K protein's stabiblity (Fig.  2G).  Special thanks to Dr. Slabicki whom generously made the original code and data available in the Method section.

- Special thanks to Dr. John Doench, Peter DeWeirdt, and Dr. Mudra Hegde (co-authors of `poola`) and the team at GPP for their commitment to open-source science. More extensive wet- and dry-lab GPP resources are available [here](https://portals.broadinstitute.org/gpp/public/).

### Installation
```
pip install pooled
```
Current version: 0.0.2

Dependencies:
- pandas
- numpy
- statsmodels
- matplotlib
- gpplot
- adjustText

### Tutorial
To see a worked example and all functionalities, please launch the Jupyter notebook using this button below:

[![nbviewer](https://img.shields.io/badge/render-nbviewer-orange.svg)](https://nbviewer.org/github/kiddo18/pooled/blob/master/notebook/pooled_implementation-Ebert-Official.ipynb)

### Future TODO!
- How to go from raw FASTQ to read count table
- Genetic Interactions (double Cas9 or Cas12a multi-gRNA)


### Acknowledgements
My knowledge of the design, execution and analyses of CRISPR screens came from my time working with my former post-doc mentor Dr. Sandor Spisak, PhD and my PI Dr. Nilay Sethi, MD, PhD at Dana-Farber Cancer Institute/HMS. This package is ultimately my trying to give back to the highly collaborative scientific and biomedical community I find myself in Boston.

### Questions?

Please leave them in the `Issues` section of this Github Repo. I'll try to compile some frequently asked questions here as they come. Highly encourage you to think about the data during the planning phase as that will greatly inform your screen design (library size, sequencing parameters, treatment time, number of treatment arms, minimum cell coverage, etc.)

### Alternatives

I would also recommend that you check out commonly used, well-cited, published tools such as MAGeCK, BAGEL, DrugZ, casTLE, CRISPhieRmix, CRISPRBetaBinomial, and many more at this [link](https://github.com/davidliwei/awesome-CRISPR).
