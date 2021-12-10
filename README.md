### The author of original papers offered the [R code](https://github.com/stegen/Stegen_etal_ISME_2013)

### [paper 1](https://www.nature.com/articles/ismej201393) and [paper 2](https://www.nature.com/articles/ismej201222)


### This R script may be more user-friendly.

```
Rscript assm_proportion.R -m otu.even.tab.txt  -t rep-seq.tree -o out -r 10

Rscript assm_proportion.R -h

	-m MAT, --mat=MAT
		Input the OTU table, with columns samples and rows OTUs, with last column is taxonomy

	-t TREE, --tree=TREE
		Input otu trees. newick format

	-o OUTDIR, --outdir=OUTDIR
		The output dirctory, default is ./

	-r REPS, --reps=REPS
		Number of randomizations, default is 999

	-p TOP, --top=TOP
		Choosing top number OTUs for calculation, default is 1000

	-h, --help
		Show this help message and exit
```

### 量化分析群落聚集过程.pdf
