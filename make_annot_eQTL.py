#!/usr/bin/env python3
"""
Reads ldsc/Builds/baseline.*.annot.gz files
- containg chr, pos, rsID, cm an dappartenance to 42 categories used in Fullikan
et al 2020- and GTEx Specificity/Data/eQTL/*.v8..txt 
Returns a similar annotation files but with a new category (eQTL) cored 1 if the 
SNP is a significant (FDR<0.05) eSNP for any gene in any tissue.
"""


## eQTL import

# Parse eQTL files through tissues from .v8.full.txt (files with rsID)
import pandas as pd
   tissue_specific_files = glob.glob(os.path.join(
            INPUT_PATH, "eQTLs_per_tissue",
            "*.v8.signif_variant_gene_pairs.txt"
    ))
    first = True
    for i, file in enumerate(tissue_specific_files):
        tissue = os.path.basename(file).split(".")[0]
        brain = bool_2_brain[tissue.startswith("Brain")]
        #print(tissue)
        iter_input = pd.read_csv(file,
                                 iterator=True, chunksize=1000000,
                                 sep="\t") #, compression="gzip")
        for chunk in iter_input:
            chunk = chunk[INPUT_COLUMNS]
            chunk.rename(columns={
                "gene_id" : "Ensembl_ID",
                "pval_beta": "p_value",
                "slope": "effect_size",
                "maf":"MAF",
            }, inplace=True)

# Liftover position from 38 to 37 
from pyliftover import LiftOver
#lo = LiftOver("hg18", "Hg38")
# or use variant_id provided in GTEx lookup table = import it and retrieve each SNP


## read annot.gz files from 1000G_Phase1 (unless we can make Phase 3 work)
pd.read_csv

## add a column + parse eQTL
## or do a join ?