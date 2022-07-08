# Network Inference

Gene co-expression network inference from gene expression data.

If you use this code, please cite the following paper: 

Identifying robust functional modules using three-body correlations in Escherichia coli,
T. Chen, P. Singh, K. E. Bassler, J. Phys. Complex. 3 015013 (2022)
https://doi.org/10.1088/2632-072X/ac5567 

The mutual information (MI) is calculated using B-Spline functions as described in Daub et al. Estimating mutual information using B-spline functions-an improved similarity measure for analysing gene expression data, BMC Bioinform. 5 118 (2004), and the relatedness scores are calculated using the context likelihood of relatedness (CLR) method described in Faith J J et al. Large-scale mapping and validation of E. coli transcriptional regulation from a compendium of expression profiles, PLoS Biol. 5 1–13 (2007).

To use the code, follow the steps below:

1. Expression data

	* Expression levels of each gene (rows) under different conditions (columns). First column is the gene name. See `example.txt`. 
	* Use bash script (work.sh) to obtain the number of genes `n` and number of experiments `N`. 
	* Example:

		`sh input-info.sh example.txt` 


2. Compile

	* Mutual information: 
		`gcc mi.c -lm -o mi.out`
	* Relatedness:
		`gcc relatedness.c -lm -o f.out`


3. Run
	* To obtain the mutual information matrix, run `mi.out` with 3 arguments, respectively as:
	* argument 1: input file name
	* argument 2: number of genes
	* argument 3: number of experiments

		`./mi.out example.txt 26 100`

	To obtain the relatedness matrix from mutual information, run `f.out` with 2 arguments, respectively as:
	* argument 1: mutual information file name
	* argument 2: number of genes

		`./f.out MI_example.txt 26`
