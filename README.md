# Creating a preliminary set of figures for 18S tag sequencing
Scripts provided follow QC_steps_V4_tagsequencing test files. Purpose of this is to generate a preliminary set of figures (and a table!) that give you a first look of your new tag sequence data. As written it is specific to 18S tag sequencing projects using the PR2 database. But, of course can be modified for other tag sequencing results.

Originally created for an in-lab tutorial for evaluating tag sequencing data.

## Prerequisites
To follow this step by step guide, use V4_OTUtable_test.txt
This is an OTU table from QIIME. Other OTU clustering programs will give you similarly formatted tables.

## Requirements
R v3.3.2 or higher
Jupyter notebook (R)

## How was V4_OTUtable_test.txt made?
See QC_steps_V4_tagsequencing repo.
### Quick review of how V4_OTUtable_test.txt was made:
Depending on the questions you want to ask your data, decide about what type of 
OTU clustering you want to do (both algorithm and percent similarity).
 
QIIME is a great resource for tutorials on OTU clustering - http://qiime.org/tutorials/otu_picking.html#running-the-otu-picking-workflows

In this example I will use open reference OTU picking (Rideout et al. 2014, Peer
J), via QIIME and the PR2 database.
```
#OTU clustering with uclust, will make a new directory (pick_open)

pick_open_reference_otus.py -i allseqs_test.fasta -o pick_open -m uclust -r /gal
adriel/sarah/PR2/pr2.qiime.fasta --suppress_step4 --suppress_taxonomy_assignment
cd pick_open

#assign taxonomy using PR2 database, creates a new directory (uclust_taxonomy)
assign_taxonomy.py -i rep_set.fna -t /galadriel/sarah/PR2/ids.names.2.txt -r /ga
ladriel/sarah/PR2/pr2.qiime.fasta -m uclust -o uclust_taxonomy

#make an OTU table, and convert to txt
make_otu_table.py -i final_otu_map.txt -o V4_tagseq_test.biom -t uclust_taxonomy/rep_set_tax_assignments.txt 
biom convert -i V4_tagseq_test.biom -o V4_OTUtable_test.txt --to-tsv --header-key=taxonomy
```
### Use a text editor to remove the comment lines from the QIIME formatted file.

## Jupyter R notebook OTUtable_to_PrelimFigs_R

Import test OTU table. This includes 5 samples (columns) with a total of 1,998 OTUs (each OTU is a row). Rscript imports OTU table and produces these three graphs which provide the user with a first look at the data!
![headcount](https://github.com/shu251/figs/blob/master/headcount_output.png)

### Total number of sequences in each sample

### Total number of OTUs in each sample

### Distribution of OTUs in each sample (i.e. how many are singleton OTUs?)
Note that the script removes global singletons (OTUs with only one sequence in the whole dataset), so the singleton OTUs shown in red are OTUs with 1 sequence in a sample, but that same OTU was detecte elsewhere.

### RData up to this point is saved as "Checkpoint1_PrelimFigs.RData", and can be loaded here.

Next the R script walks the user through removing unwanted samples (based on your preliminary look at the data). For instance, if one of your samples has too few sequences or too few OTUs, this is where you can remove that sample.

## Randomly subsample data
It is required that you randomly subsample data so that all samples have the same number of sequences for any calculation of diversity estimates. Again, depending on your question, it may also be advisable to use subsampled data for all downstream analyses. Here, I'm using the vegan package to do this. 

## Calculate alpha diversity
Using the vegan package, I calculated Shannon and inverse Simpsons diversity metrics, and the total number of OTUs per sample (as a proxy for species richness).
https://cran.r-project.org/web/packages/vegan/vegan.pdf

Comparing diversity metrics between samples

In this test example, Sample 1 appears to be more diverse than the other samples.

## Simplifying PR2 taxonomy output and making a community composition plot 
The output from the PR2 database has a long list of taxa to the approximately species name. These next steps demonstrate how you can get a quick view of the whole community at the approximately phylum or class level. However, depending on your project/biome/advisor/personal preference, you are likely going to want to highlight specific groups or NOT highlight certain taxonomic groups. Since my lab focused on several of the major taxonomic groups found in marine ecosystems (like off the coast of SoCal!), I created a function to simplify the taxonomic names. This is a great way to get a 'birds eye' view of the data and give your advisor a general idea of what the whole community looks like.
Function specifically separates out the PR2 output into Levels and then takes the information from the various levels to fill in a newly created "Taxa" column. Then I aggregate the data, or sum the data, by the taxonomic group designated by "Taxa". 

## Community composition - relative abundance of sequences





