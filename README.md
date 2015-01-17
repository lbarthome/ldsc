LDSC (LD SCore)
===============

v0.0001 (alpha)
Warning: under very active development

LDSC is a command line tool for estimating heritability and genetic correlation from GWAS
summary statistics. LDSC also computes LD Scores.

Support
-------

Before contacting us, please try the following:

1. Common issues are described in the [FAQ](docs/FAQ)
2. The methods are described in the papers (see citations below)
3. Search the [issue tracker](https://github.com/bulik/ldsc/issues)

Please report bugs on the [issue tracker](https://github.com/bulik/ldsc/issues). Also feel free to ask about statistical issues using the issue tracker. The issue tracker is searchable, and somebody else might have the same problem as you in the future. Better to archive the answer in a searchable location than lock it away in a private email. 

Citation
--------

For now, please cite

Bulik-Sullivan, et al. LD Score Regression Distinguishes Confounding from Polygenicity in Genome-Wide Association Studies.
In Press at Nature Genetics. ([bioRxiv Version](http://biorxiv.org/content/early/2014/02/21/002931))

We are currently preparing manuscripts describing the methods for estimating partitioned h2 and rg.

Requirements
------------

1. Python 2.7
2. argparse 1.2.1
3. bitarray 0.8.1
4. numpy 1.8.0
5. pandas 0.15.0
6. scipy 0.10.1

License
-------

This project is licensed under GNU GPL v3.


Authors
-------

Brendan Bulik-Sullivan (Broad Institute)

Hilary Finucane (MIT Math)
