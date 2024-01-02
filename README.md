


## About the remote

This is the remote in which we host all the workflows and utils needed to perform TWMR analyses that allows us to make causal inferences from GWAS and (sc)QTL-derived transcriptome wide data.  
Finally, we provide multiple approaches to infer which cell type has the highest burden of differentially expressed disease causing genes, based on GWAS-TWMR results.

The TWMR publication can be found [here](https://www.nature.com/articles/s41467-019-10936-0)


## Workflow

Below, we show a diagram summarising the entire workflow developed to perform cell type enrichment TWMR analyses.  

![workflow_img](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/cte_twmr_diagram.png)





## Getting Started


### Set up docker container to run the workflow

First step to get started performing TWMR and subsequent TWMR analyses is to get a working Docker container
having all the dependencieds needed to run the workflow


```bash
cd docker
docker build -t amc-celltyping .
bash dockerrun_twmrcontainer.sh
```


### Run each workflow step


Once docker is installed, you can run each workflow step following the R markdowns contained in the workflow directory.  
Please, note that some external data will need to be allocated on your working directory to be able to go through eachs step successfully.  
Below I provide a summary what each R markdown allows us to do.  




* [Do QC on eQTL data](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/1.eqtlfiles_qc.Rmd). Get eQTL data and QC it to match format needed for downstream analyses as well as meet the QC criteria summarised in the diagram.  

* [Do QC on GWAS data](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/2.gwasfiles_qc.Rmd). Get GWAS data and QC it to match format needed for downstream analyses as well as meet the QC criteria summarised in the diagram.  

* [Merge eQTL and GWAS data](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/3.eqtlgwasmerge.Rmd). Merge QCed eQTL and GWAS data to keep overlapping SNPs and get ready to run GCTA-COJO to find independently associated SNPs with an expression qt-trait

* [Get TWMR matrices](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/4.matrix_generate.Rmd). In this part of the workflow, we do some key steps. First, we run COJO to get a matrix of independent SNPs. Then, we prune results and also perform a matrix expansion job to allocate extra genes for the multiple instruments from COJO.  

* [Run TWMR](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/5.run_twmr.Rmd). In this part of the workflow, we run TWMR

* [Cell type enrichment](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/6.enrichment_analysis.Rmd). Here, we perform the enrichment analyses from TWMR results

* [Compare gene sets](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/7.compare_geneSets.Rmd). We compared results against state of the art methods to perform cell type enrichment analyses



### Get matrices and run TWMR through an Rscript

Based on our experience, the generation of gene matrices that we will finally use on inferences of causality based on TWMR is time consuming.  
For that, we have generated an Rscript stored in utils called **prepare_twmrinput.R**.  
In addition, we have creates an SBATCH script to further automate the process in case we need to repeat the process in multiple input GWASs at a time.  

To run:

```bash
sbatch utils/slurm_prepare_twmr.sh
```

Finally, we are ready to take the outputs the output of this job to then perform [Cell type enrichment](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/6.enrichment_analysis.Rmd) and [Compare gene sets](https://github.com/AMCalejandro/celltype_twmr/blob/master/img/7.compare_geneSets.Rmd).  