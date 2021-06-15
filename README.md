# shotgun_microbiome_pipeline

## Shotgun metagenomic analytical pipeline

This pipeline was made to analyze shotgun microbial read data, including OTU assignment, functional profiling, diversity assessments, largely based off of tools from the biobakery suite. This is not meant to replace the information available there (cited and linked below); rather it provides a compilation of resources for the usage of those methods in upcoming publications (to be listed).

## Required software

* [biobakery](https://github.com/biobakery/biobakery)
  * [KneadData](https://github.com/biobakery/biobakery/wiki/kneaddata)
  * [MetaPhlAn](https://github.com/biobakery/MetaPhlAn)
  * [HUMAnN](https://github.com/biobakery/humann)
  * [PPANINI](https://github.com/biobakery/biobakery/wiki/ppanini)
* [R](https://www.r-project.org/)
  * [SourceTracker](https://github.com/danknights/sourcetracker)
  * [phyloseq](https://github.com/joey711/phyloseq)
  * [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)


## Pipeline

### Filtering

Reads are first filtered against a human-specific database ([KneadData](https://github.com/biobakery/biobakery/wiki/kneaddata) links to one option) as well as laboratory or environmental contaminants as appropriate. This can be done by creating a custom index as specified in [KneadData](https://github.com/biobakery/biobakery/wiki/kneaddata), and then filtering all samples against that index.

### Profiling

Profiling is done in [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) with "-t rel_ab_w_read_stats" specified. Check that this file has two header lines (#SampleID and #clade_name) or edit to ensure that this is the case. If using MetaPhlAn 3, for example, a script like the following may help convert the header:

```for f in *txt 
 do 
 sed -i '1d' $f
 sed -i '1d' $f 
 sed -i '2d' $f 
 sed -i '3d' $f
 done
```

After this, an OTU table is made using [this R script](merge.table.R), which makes this data easy to import/use in other contexts. Note that this R script can merge additional microbiome data into the total OTU table, so long as it fits 1 of 2 formats:

1. It is a profile with taxa names similar to the MetaPhlAn outputs
2. It is a profile with OTU's numbered instead of named, and has an associated table linking the OTU numbers to the taxonomic identifier

These varying files can be placed in the PerSample (1) or Paired (2) folders that the script references. If no samples are being used in the Paired folder, leave a set of paired empty files there. Further, using the $estimated_number_of_reads_from_the_clade statement will generate an OTU table with read counts (better for alpha diversity calculations) while using $relative_abundance will generate an OTU table with relative abundance. 

### Analyses

[SourceTracker](https://github.com/danknights/sourcetracker) can take the generated tables as an input. Some trimming may be necessary (the last column on merged_otu_lineage.txt is extraneous information, and if no paired samples are used, the empty sample column can be removed from both generated files). A metadata file is needed (as per the layout indicated in SourceTracker at minimum).

[Phyloseq](https://github.com/joey711/phyloseq) requires a taxonomy table (merged_otus_lineage.txt), an OTU table (merged_otus.txt), and a metadata file. The editing done above is also useful for this import. This allows for easy tests on alpha diversity via the `plot_richness()` function. Beta diversity can also be explored with `ordinate()` and `plot_ordination`.

[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) allows for a large number of operations - the following vignette describes many of these well: (https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). For the listed publications, variance stabilized transformations with the `vst()` function were useful for transforming high-throughput sequence count data. The differential abundance can be done in combination with `phyloseq` objects as described here (https://joey711.github.io/phyloseq-extensions/DESeq2.html). 

[HUMAnN](https://github.com/biobakery/biobakery/wiki/humann2) can be run on the original FASTQ (ideally trimmed and filtered) files. This generates multiple table outputs (gene family abundance, pathway abundance, and pathway coverage). There are multiple settings and options to edit here, including databases for alignment, thresholds, post-processing (how tables are sorted), that may all be altered depending on the study conditions. Normalized and stratified tables are similar in organization to the OTU tables produced earlier, and can be used as inputs for other analyses above if appropriate for the study - like being imported into phyloseq for ordination or differential abundance.

[PPANINI](https://github.com/biobakery/biobakery/wiki/ppanini) highlights functional novelty by plotting known and unknown components from samples. It utilizes the gene abundance table from [HUMAnN](https://github.com/biobakery/biobakery/wiki/humann2) with no modification, and the resulting output `.txt`. file can be imported into R for more flexible graphing. The package comes with an internal plotting mechanism through matplotlib, which is another option for visualization.

### Citations

Knights, D. et al. 2011. Bayesian community-wide culture-independent microbial source tracking. Nature Methods, 8, 9, 761-763. [Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3791591/)

Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8 [Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

McIver LJ, Abu-Ali G, Franzosa EA, Schwager R, Morgan XC, Waldron L, Segata N, Huttenhower C. bioBakery: a meta'omic analysis environment. Bioinformatics. 2018 Apr 1;34(7):1235-1237. [PMID: 29194469](https://pubmed.ncbi.nlm.nih.gov/29194469/)

McMurdie and Holmes (2013) phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data PLoS ONE 8(4):e61217 [Paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061217)
