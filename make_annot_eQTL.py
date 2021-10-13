#!/usr/bin/env python3
"""
Reads ldsc/Builds/baseline.*.annot.gz files
- containg chr, pos, rsID, cm an dappartenance to 42 categories used in Fullikan
et al 2020- and GTEx Specificity/Data/eQTL/*.v8..txt 
Returns a similar annotation files but with a new category (eQTL) cored 1 if the 
SNP is a significant (FDR<0.05) eSNP for any gene in any tissue.
"""

import os.path
import pandas as pd
import time as time

## ----------- Functions
def get_esnps(chr, INPUT_PATH_EQTL):
    """Takes a chromosome chr (just an integer) and returns the corresponding list of 
    eSNPs in a panda datafrane contianing one column of rsID (called SNP) and a column
    of logical value for beeing a siginficant eQTL (all here are !)"""
    INPUT_COLUMNS = [
        "rsID" #, "ref", "alt"
    ]    
    input_eQTL_filename = os.path.join(INPUT_PATH_EQTL, "chr"+chr+".esnps.txt")
    esnps = pd.read_csv(input_eQTL_filename , sep="\t")[INPUT_COLUMNS]
    esnps["GTEx_v8_eSNP"] = 1 # logical annotation that these snps are eSNPs
    esnps.rename(columns={
        "rsID" : "SNP",
    }, inplace=True)
    return esnps

def get_baseline_annot(chr, INPUT_PATH_BASELINE):
    """Takes a chromosome chr (just an integer) and returns the corresponding list of 
    SNPs annotated in ldsc baseline model in a panda datafrane containing columns to
    identify the SNP (chr, bp, snp, cm) and then one column per category with logical 
    values (0 or 1) depending on the SNP appartenance to the given category"""
    input_baseline_filename = os.path.join(
            INPUT_PATH_BASELINE,
            "baseline."+chr+".annot.gz"
        )
    baseline_annot = pd.read_csv(input_baseline_filename,sep="\t",  compression="gzip")
    return baseline_annot

def build_new_annot(esnps, baseline_annot):
    """Takes list of eSNPs and existing baseline annotation tables and returns
    a annotation table with a new column corresponding to eSNPs. An eSNP significant
    in GTEx will have a 1, others a 0 (not using a list of tested SNPs in GTEx"""
    new_annot = pd.merge(
        baseline_annot, esnps, 
        on="SNP", # the column GTEx_v8_eSNP will be Nan if SNP not in esnps
        how="left" # only keep rows from original annotation files
    )
    new_annot["GTEx_v8_eSNP"] = new_annot["GTEx_v8_eSNP"].map(
        lambda x: str(0) if pd.isna(x) else str(1)
        )
    return new_annot

# ----------------------------------------------------------------------
# ------------------------------ Main ----------------------------------
# ----------------------------------------------------------------------

if __name__ ==  "__main__":
    # if the scrits is executed (python, ipython), __name__ == main
    # if the scripts is imported it has module name, only import functions
    startTime = time.time()

    ## --------------- Getting input file and output file PATHs
    # args = sys.argv

    # if len(args) <1:
    #     print("Error : missing argument. Try again with 1 argument \
    #         (input folder)")
    #     sys.exit(1)

    PATH = os.path.dirname(__file__)
    INPUT_PATH_EQTL= os.path.join(PATH, "Data/eQTL") # to input folder
    INPUT_PATH_BASELINE = os.path.join(PATH, "Builds/baseline")
    OUTPUT_PATH = os.path.join(PATH, "Builds/baseline_eQTL/")


    ## -------------- Importing input
    for chr in range (1,23):
        chr = str(chr)
        print("Annotating chromosome "+chr+"...\r")

        # Import esnps
        esnps = get_esnps(chr, INPUT_PATH_EQTL)
        
        # Import baseline model
        baseline_annot = get_baseline_annot(chr, INPUT_PATH_BASELINE)

        # Build new annotations
        new_annot = build_new_annot(esnps, baseline_annot)

        # Export annotations
        output_filename = os.path.join(OUTPUT_PATH, "baseline."+chr+".annot.gz")
        new_annot.to_csv(
                    output_filename,
                    header=True, mode='w', sep="\t", compression="gzip",
                    index=False
                )
                                
    executionTime = (time.time() - startTime)
    print('Done ! Execution time in: ' + str(executionTime) + " seconds")