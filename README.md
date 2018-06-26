# Generate a preliminary set of figures from tag sequencing data
The purpose of this is to create a preliminary set of figures (and a table!) in R. This will give you a first look of your new tag sequence data. As written it is specific to 18S tag-sequencing projects using the PR2 database. But, of course can be modified for other tag sequencing results.
Originally created for an in-lab tutorial to demonstrate the use of R, provie an example of initial R-based analyses to apply to tag sequence results, and generate figures to ward off an anxious PI.

## Requirements
* R v3.3.2 or higher
* (optional) Jupyter notebook (R)
* To follow directly, use products from [OTU clustering in qiime1 tutorial](https://github.com/shu251/V4_tagsequencing_18Sdiversity_q1) or [OTU clusering in qiime2](https://github.com/shu251/V4_tagsequencing_18Sdiversity_qiime2). Otherwise, ensure you have an OTU table, where each row is an OTU (with a unqiue ID) and each column is a sample (first row is a header/sample name) or follow along with provided text file: 'V4_OTUtable_test.txt'

## Jupyter R notebook OTUtable_to_PrelimFigs_R

Import test OTU table (same type of OTU table provided by QIIME). This includes 5 samples (columns) with a total of 6138 OTUs (each OTU is a row). Rscript imports OTU table and produces these three graphs which provide the user with a first look at the data!
Script then removed global singletons and give you stats on how many OTUs and sequences are in your dataset.
![headcount](https://github.com/shu251/figs/blob/master/headcount_output.png)

### Total number of sequences in each sample
![totalseqs](https://github.com/shu251/figs/blob/master/seq_stats_graphs.png)

### Total number of OTUs in each sample
![TotalOTUs](https://github.com/shu251/figs/blob/master/totalOTUs.png)

### Distribution of OTUs in each sample (i.e. how many are singleton OTUs?)
Note that the script removes global singletons (OTUs with only one sequence in the whole dataset), so the singleton OTUs shown in red are OTUs with 1 sequence in a sample, but that same OTU was detecte elsewhere.
![OTUdist](https://github.com/shu251/figs/blob/master/OTUdistribution.png)

### RData up to this point is saved as "Checkpoint1_PrelimFigs.RData", and can be loaded here.

Next the R script walks the user through removing unwanted samples (based on your preliminary look at the data). For instance, if one of your samples has too few sequences or too few OTUs, this is where you can remove that sample.

### Notes on normalization
The product of OTU or ASV clustering will be compositional data. There are better and better methods for dealing with this type of data, but at the end of the day, always keep this in mind with your interpretations.

Recommended references:
* Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V. & Egozcue, J. J. Microbiome Datasets Are Compositional: And This Is Not Optional. Front. Microbiol. 8, 57–6 (2017).
* Weiss, S. et al. Normalization and microbial differential abundance strategies depend upon data characteristics. Microbiome 5, 1–18 (2017).
* McMurdie, P. J. & Holmes, S. Waste Not, Want Not: Why Rarefying Microbiome Data Is Inadmissible. PLoS Comput Biol 10, e1003531 (2014).

For this example in R, we will stick with randomly subsampling our data. This normalization step serves to account for sample-to-sample differences in the raw sequence results. I'll continue updating the R script here as much as I can. 

## Randomly subsample data
It is required that you randomly subsample data so that all samples have the same number of sequences for any calculation of diversity estimates. Again, depending on your question, it may also be advisable to use subsampled data for all downstream analyses. Here, I'm using the vegan package to do this. 

## Calculate basic alpha diversity indices
Using the vegan package, I calculated Shannon and inverse Simpsons diversity metrics, and the total number of OTUs per sample (as a proxy for species richness).
https://cran.r-project.org/web/packages/vegan/vegan.pdf

Comparing diversity metrics between samples
![alpha](https://github.com/shu251/figs/blob/master/alpha_div.png)
In this test example, Sample 1 appears to be more diverse than the other samples.

## Generate rarefaction curve
This type of figure gives you an *idea* of how many OTUs you should *expect*. By continually subsampling your data (x-axis), how many OTUs do you get as a result (y-axis). *Ideally*, sampling/sequencing effort can be called good enough if your rarefaction curve reaches an asymptote. 

The example provided here does not even out because I've subsampled the data for run time.

![rare](https://github.com/shu251/figs/blob/master/rare.png)

## Simplifying PR2 taxonomy output and making a community composition plot 
The output from the PR2 database has a long list of taxa to the approximately species name. These next steps demonstrate how you can get a quick view of the whole community at the approximately phylum or class level. However, depending on your project/biome/advisor/personal preference, you are likely going to want to highlight specific groups or NOT highlight certain taxonomic groups. Since my lab focused on several of the major taxonomic groups found in marine ecosystems (like off the coast of SoCal!), I created a function to simplify the taxonomic names. This is a great way to get a 'birds eye' view of the data and give your advisor a general idea of what the whole community looks like.
This R function specifically separates out the PR2 output into Levels and then takes the information from the various levels to fill in a newly created "Taxa" column. Then I aggregate the data, or sum the data, by the taxonomic group designated by "Taxa". 

## Community composition - plotting

Several options for plotting community composition using ggplot. I've included barplots representing the relative abundance of reads (pictured below), barplots that demonstrate total reads for a given taxonomic group, and and area plot.  

![Community composition](https://github.com/shu251/figs/blob/master/CommunityComposition.png)

## Application of Bray-Curtis dissimilarity
This is typically used to evaluate the beta diversity of a data set. In this short example, the relative abundance of OTUs in each sample was computed. Then I demonstrate two ways to visualize the data, first to generate an average hierarchical cluster dendrogram. Then use the vegan package command "metaMDS()" to calculate a Bray-Curtis dissimilarity matrix and plot in multidimensional space (using either base R plot or extracting the points from the NMDS and plotting in ggplot). 

# Last updated June 25th, 2018 - S. Hu

