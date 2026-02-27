GO and DEGs
------------
  `GODEG` is a python and R script pipeline to annotate GO terms using UniProt database information, apply cutoffs for DEG identification, classify DEGS in the complete GO hierarchy and perform an enrichment analysis to identify overrepresented GO terms containing DEGs. This pipeline is a straightforward method to integrate DEGs and GO information.
  
  `GODEG` uses two main file inputs, a DEG table containing fold change and adjusted p values (resulted from software as DESEQ2 and EdgeR) and the genome GO annotation file. The DEG table is filtered to obtain DEGs based on adjusted p value and foldchange cutoffs, and these DEGs are assigned to GO categories (specific and parental categories). Typically, GO annotations contain the most specific GO terms. To obtain the ancestor GO categories we used the R package GO.db from bioconductor. After assigning the DEGs to GO terms, an enrichment analysis is performed to identify overrepresented GO terms containing DEGs.
  
  We also provide a python script to obtain GO annotation (`go_uniprot_annotator`) from a blast search against UniProt Database. This is achieved by parsing the UniProtIDs from a blast tabular file in the UniProt API to recover GO terms and protein names. `go_uniprot_annotator` automates the [Retrieve/ID mapping](https://www.uniprot.org/id-mapping) from UniProt, and provides GO annotation files for Biological Process, Molecular Function and Cellular Component, and also provides the protein name for each best hit from blast file. It is also possible (and optional) to inform a blast against SwissProt to provide a second annotation for the protein names. These annotation files are ready to use with `GODEG`.

How to install
--------------
Using a clean ubuntu 25.10, run the following commands using sudo. This pipeline is not resource intensive and can be run in most systems.
```
sudo apt update
sudo apt upgrade
```
Then, all dependencies can be installed in two commands:
```
sudo apt install python3 python3-pandas python3-scipy python3-statsmodels python3-xlsxwriter libcurl4-openssl-dev r-base python3-aiohttp
sudo Rscript -e 'if (!require("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install("GO.db")'
```
Alternatively, these dependencies can also be installed with `pip3`.

After installing the dependencies above, try:
```
python3 GODEG_v1.py -h
```
This will display the parameters descriptions if all dependencies were detected, and it will complain if there is a missing dependency. The R version must be >=4 and contains the latest GO.db package installed. The python version must be 3, and the following libraries must be installed manually and are not by default in python3: aiohttp, numpy, pandas, scipy, statsmodels and tqdm. The complete library list is: aiohttp, argparse, asyncio, concurrent.futures, functools, io, numpy, os, pandas, re, scipy.stats, shutil, statsmodels, subprocess, sys and tqdm.

Default usage
--------------
The default command for `GODEG` is:
```
python3 GODEG_v1.py -d deg_table.tsv -go genes_GO_biological_process.tsv -go_mode P --samples samples.txt -a uniprot_swissprot_protein_annotation.tsv -t 12
```
* The deg_table.tsv must be a table (tab separated), containing in the first column the seqIDs from the reference genome, followed by the remaining statistics as adjusted p value and foldChange, typically resulted from software as DESEq2 or EdgeR. This table can be pre-filtered to hide non-expressed genes, or unfiltered containing all the genome seqIDs. You can merge several comparisons in the same table (a common practice since they share the same seqIDs in the first column).
* genes_GO_biological_process.tsv (with header) is the GO annotation for the complete set of genes of the reference genome. It is composed of two columns: the first is the seqIDs from the reference genome, and the second column the corresponding GO terms comma separated.
* IMPORTANT: Be sure that the seqIDs in the first column of deg_table.tsv and genes_GO_biological_process.tsv matches. Be sure that both files seqIDs follow the same structure, and remove suffixes not shared by both seqIDs lists (i.e.: .1, .2, common to refer sequence version).
* samples.txt is a three column tsv file (without header) used to filter deg_table.tsv. samples.txt contains in the first column the foldchange column name from deg_table.tsv, in the second column the p value column name and the third column the sample name (the alias for the pvalue and foldchange pair treatment). It is necessary to inform those because deg_tables.tsv contains other statistics (i.e.: normalized reads counts) not used by this package. It is also vital for the pipeline properly identify foldChange and adjusted p value columns to identify DEGs. Below a samples.txt example for a deg_table.tsv with three sample comparisons:
```
sample1_log2FoldChange  sample1_padj  sample1
sample2_log2FoldChange  sample2_padj  sample2
sample3_log2FoldChange  sample3_padj  sample3
```
* go_mode is to choose among the three major go categories": P (biological process), F (molecular function) and C (cellular component).
* uniprot_swissprot_protein_annotation.tsv (with header) contains the protein names for the genome seqIDs. The first column contains the genome seqIDs, and the following column the annotation information. This parameter integrates the protein names in the seqIDs from the deg_table, and it is not mandatory. If you skip this parameter, it is highly advised to annotate the deg_table.tsv before.
* -t 12 is the number of threads (or processor cores) to run the pipeline.

If you wish to annotate the GO terms for the genome seqIDs using `go_uniprot_annotator`:
```
python3 go_uniprot_annotator_v1.py -bu blast_against_uniprot.tsv -bs blast_against_swissprot.tsv -t 200
```
* blast_against_uniprot.tsv is a tabular blast file (fmt 6) of the genomes seqIDs against UniProt. This blast table must be filtered by your annotation cutoffs (i.e.: e-value). We recommend using the most update Uniprot version as GO terms are constantly updated.
* blast_against_swissprot.tsv is a tabular blast file (fmt 6) of the genomes seqIDs against SwissProt. Provides a second annotation to the protein names, more accurate but less comprehensive than UniProt. This parameter is optional.
* -t 200 is the number of hits considered by each seqID. The higher the number, the more number of GO terms recovered for a seqID, but also significantly increases the run time of the annotation.
* This script needs internet access to recover the GO terms and protein names from UniProt server.

Result files
--------------

`go_uniprot_annotator` will output blast_input_GO_biological_process.tsv, blast_input_GO_cellular_component.tsv, blast_input_GO_molecular_function.tsv and uniprot_swissprot_protein_annotation.tsv (or uniprot_protein_annotation.tsv if -bs is not provided). The first three are the GO annotation files for each GO major category, and the fourth is the proteins names annotation.

`GODEG` will result in seven main outputs:
* deg_table_DEG.tsv (output 1): applies the FoldChange and adjusted p value cutoffs to filter DEGs. Next to each foldChange column there is a tag: UpRegulated for upregulated genes, DownRegulated for downregulated genes, NotDEG for genes that passed the adjusted p-value cutoff but not the foldChange cutoff, and null for genes that didn't pass the adjusted p-value cutoff. By default, the adjusted p-value cutoff is 0.05, foldChange upregulated > 1 and foldChange downregulated <-1.
* filter_deg_table_DEG.tsv (output 2): this output is a filtered version of output 1, containing only DEGs (UpRegulated and DownRegulated) foldChange values. seqIDs without at least one DEG are removed. This output is handy to build for heatmap figures.
* degs_count.tsv (output 3): contains the gene counts of DownRegulated, UpRegulated and total DEGs per sample.
* complete_GO_BP_classification_deg_table.xlsx and complete_GO_BP_classification_deg_table.tsv (output 4 and 5): these are the final GO classification for all genes (genome complete gene set and DEGs). The first column (GO Term) has the GO term number, the second (TERM) the GO term name, the third (Genome_seqid) the genome seqIDs belonging to that GO category and the forth column (Genome_seqid_count) the number of seqIDs from the genome in this GO category. Following those columns, for each treatment (third column of the samples.txt file) there are the seqIDs and seqIDs counts for UpRegulated and Downregulated Genes, the total of DEGs in that category (DEG_treatment_seqID_count), the enrichment p-value (DEG_sample_f_test_p_value) followed by the enrichment adjusted (p-value DEG_treatment_f_test_FDR). To browse this table it is advised to use the .xlsx file, as it hides the extensive seqIDs lists (which can be obtained in the .tsv) and has a pre-configured table view. In this table it is possible to identify GO categories containing most downregulated or upregulated genes and identify enriched GO categories.
* degs_by_go_term (output 6): inside this folder, there is a table file for each GO category present in output 4 and 5 containing the foldchange values. It uses output 2 table (containing only foldChange of DEGs), but filtered for every GO category. This is also useful to build heatmaps, as the user can straight away make a heatmap from a specific GO category annotation with the foldChange values.
* enriched_only.xlsx (output 7): contains only enriched terms with enrichment p-value >0.05 (default). GO terms without at least one enriched treatment are not present in this file.

About
--------------
`GODEG` and `go_uniprot_annotator` were developed by Dr. Gustavo Lima Rodrigues, supported by the Genomics Group (ESALQ, University of São Paulo) and Dr. Fei Lab (BTI, Cornell University). The academic paper about this tool is under development.
