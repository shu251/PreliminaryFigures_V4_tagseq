
#Load required libraries
library(reshape2)
library(ggplot2)
library(vegan)
library(plyr)

#Import OTU table
count<-read.table('V4_OTUtable_test.txt', sep="\t",header=TRUE, skip=1, comment.char = "")
# quick view:
colnames(count)[1]<-"OTU.ID"
head(count) #make sure samples are column names and OTU.IDs (from PR2) are row names
dim(count) #V4_OTUtable_test.txt should be 6138 rows and 7 columns


##Get quick stats of OTU results
length(count$OTU.ID) #Total number of OTUs generated
colsum<-apply(count[2:6],2,sum) #Colums 2:5 are my sample colums
colsum #number of sequences per sample

##Filter out OTUs with only 1 sequence, singletons
rowsum<-apply(count[2:6],1,sum) #remove global singletons.
count.no1 = count[ rowsum>1, ]  #count.no1 = OTU table without global singletons
dim(count)[1] - dim(count.no1)[1] #Outputs the number of OTUs (total) lost in this step

#Isolate only count colums
counts_only<-count.no1[2:6]
seq_total<-apply(counts_only,2,sum) #number of sequences per sample
OTU_count<-colSums(counts_only>0) #total number of OTUs
OTU_single<-colSums(counts_only==1) #Number of singleton OTUs
OTU_double<-colSums(counts_only==2) #Number of doubleton OTUs
OTU_else<-colSums(counts_only>2) #Number of OTUs with more than 2 sequences

sample_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_else)# Compile
head(sample_info) #dataframe with OTU stats per sample
# write.csv(sample_info, file="OTUstats.csv") #Option to write out to table

#Plot it
sample_info$samples<-row.names(sample_info)
allM<-melt(sample_info) #Melt to long format
allM$Figure<-"Sequences"
allM$Figure[allM$variable == "OTU_count"]<-"Total OTUs"
otudist<-c("OTU_single", "OTU_double", "OTU_else")
allM$Figure[allM$variable %in% otudist]<-"Breakdown of OTUs"
#head(allM)


#Basic bar plot
bar_stats<- ggplot(allM, aes(x=samples, y=value, fill=variable))+geom_bar(stat="identity",position="stack",color="black")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "top")

bar_stats %+% subset(allM, Figure %in% "Breakdown of OTUs")+labs(title="Distribution of OTUs",x="Samples", y="Total OTUs")+scale_fill_manual(values=c("#e41a1c","#fee08b","#4393c3"))
bar_stats %+% subset(allM, Figure %in% "Total OTUs")+labs(title="Total number of OTUs",x="Samples", y="Total OTUs")+scale_fill_manual(values=c("darkblue"))
bar_stats %+% subset(allM, Figure %in% "Sequences")+labs(title="Total number of sequences",x="Samples", y="Total sequences")+scale_fill_manual(values=c("purple"))

#save(counts_only, count.no1, allM, file="Checkpoint1_PrelimFigs.RData") #Optional save R objects
#load("Checkpoint1_PrelimFigs.RData",verbose=T) #Option to load R objects from previous

#Based on the total number of sequences in each sample.
count.no1$Sample_2<-NULL
head(count.no1)

#Optional code to subsample
#Assign row names as OTU IDs
row.names(count.no1)<-count.no1$OTU.ID
#Isolate only columns with data
keep<-count.no1[2:5];head(keep)
sub<-min(colSums(keep)); sub #sub =total number of sequences that is fewest among all samples

#Requires vegan library
rare <- rrarefy(t(keep), sub) #Randomly subsamples data so that all sample have the same number of sequences
subsampled<-as.data.frame(t(rare))
colSums(subsampled) #all should be equal to sub (in test data, n=99,530)

#Calculate alpha diversity

#diversity measurement that accounts for both abundance and evenness of the species present (both evenness and richness). proportion of species relative to the total multiplied by the ln of the proportion
shannon<-diversity(subsampled,index="shannon",2)
#evennes in the community
invsimp<-diversity(subsampled,index="invsimpson",2)

OTU_count<-colSums(subsampled>0) #to evaluate species richness
alpha<-data.frame(shannon,invsimp,OTU_count) #combine measurements
head(alpha)

alpha$samples<-row.names(alpha)
alpha.m<-melt(alpha)
head(alpha.m)

ggplot(alpha.m, aes(x=samples, y=value, fill=variable, shape=variable))+geom_point(size=4, aes(color=variable))+facet_grid(variable~.,scales="free")+theme_bw()+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "top")

#head(count.no1) #original data, has OTU and taxonomy information
#head(subsampled) #subsampled data for analysis
key<-count.no1[c(1,6)] #OTU ID to taxonomy key
subsampled$OTU.ID<-row.names(subsampled)

#Get taxonomic information (by OTU.ID) back on subsampled data
count.subsampled<-join(subsampled, count.no1[c(1,6)], by="OTU.ID", type="left", match="first")
head(count.subsampled); dim(count.subsampled)

#This is specific to the PR2 database. 
#Depending on your scientific question, you will need to review the taxonomic grouping
#Function, makes a simplifed "Taxa" column based on PR2 output

pr2_rename_taxa<-function(df){
  library(reshape2)
  split<-colsplit(df$taxonomy, "; ", c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8","Level9", "Level10", "Level11", "Level12"))
  split[ is.na(split) ] = "XXX"
  split[ split == "" ] = "XXX"
  split$Taxa<-"Other/unknown"
  split$Taxa[split$Level1 == "No blast hit"]="No blast hit"
  split$Taxa[split$Level1 == "Unassigned"]="Unassigned"
  split$Taxa[split$Level1 == "None"]="None"
  split$Taxa[split$Level2=="Amoebozoa"]="Amoebozoa"
  split$Taxa[split$Level2=="Apusozoa"]="Other/unknown"
  split$Taxa[split$Level2=="Eukaryota_X"]="Other/unknown"
  split$Taxa[split$Level2=="Eukaryota_Mikro"]="Other/unknown"
  split$Taxa[split$Level2=="Stramenopiles"]="Stramenopiles-Other"
  split$Taxa[split$Level2=="Alveolata"]="Alveolates-Other"
  split$Taxa[split$Level2=="Opisthokonta"]="Opisthokonts-Other"
  split$Taxa[split$Level2=="Archaeplastida"]="Archaeplastids-Other"
  split$Taxa[split$Level2=="Excavata"]="Excavates"
  split$Taxa[split$Level2=="Rhizaria"]="Rhizaria-Other"
  split$Taxa[split$Level2=="Hacrobia"]="Other/unknown"
  split$Taxa[split$Level3=="Haptophyta"]="Haptophytes"
  split$Taxa[split$Level3=="Fungi"]="Opisthokont-Fungi"
  split$Taxa[split$Level3=="Metazoa"]="Opisthokont-Metazoa"
  split$Taxa[split$Level3=="Foraminifera"]="Rhizaria-Foraminifera"
  split$Taxa[split$Level3=="Dinophyta"]="Alveolates-Dinoflagellates"
  split$Taxa[split$Level4=="Syndiniales"]="Alveolates-Syndiniales"
  split$Taxa[split$Level3=="Cryptophyta"]="Cryptophytes"
  split$Taxa[split$Level3=="Ciliophora"]="Alveolates-Ciliates"
  split$Taxa[split$Level3=="Chlorophyta"]="Archaeplastids-Chlorophytes"
  split$Taxa[split$Level3=="Cercozoa"]="Rhizaria-Cercozoa"
  split$Taxa[split$Level4=="Acantharea"]="Rhizaria-Acantharia"
  split$Taxa[split$Level4=="Chrysophyceae-Synurophyceae"]="Stramenopiles-Chrysophytes"
  split$Taxa[split$Level4=="Pelagophyceae"]="Stramenopiles-Pelagophytes"
  split$Taxa[split$Level4=="Bacillariophyta"]="Stramenopiles-Diatoms"
  split$Taxa[split$Level4=="MAST"]="Stramenopiles-MAST"
  split$Taxa[split$Level4=="Polycystinea"]="Rhizaria-Polycystines"
  split$Taxa[split$Level4=="RAD-C"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-B"]="Rhizaria-RAD (A,B,C)"
  split$Taxa[split$Level4=="RAD-A"]="Rhizaria-RAD (A,B,C)"
  return(split)
} 
newtax<-pr2_rename_taxa(count.subsampled)
unique(newtax$Taxa) #Simplified taxonomic group naming schematic
data_binned<-data.frame(count.subsampled, newtax) 

data.m<-melt(data_binned) #melt
#head(data.m)
data.agg<-aggregate(data.m$value, by=list(Taxa=data.m$Taxa,Samples=data.m$variable),sum) #sum sequences by taxonomic group
save(data.agg, data.m, data_binned, file="Checkpoint2_PrelimFigs.Rdata")

tax_order=c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Alveolates-Syndiniales","Alveolates-Other","Archaeplastids-Chlorophytes","Archaeplastids-Other","Cryptophytes","Excavates","Haptophytes","Rhizaria-Acantharia","Rhizaria-Cercozoa","Rhizaria-Polycystines","Rhizaria-RAD (A,B,C)","Stramenopiles-Diatoms","Stramenopiles-Pelagophytes","Stramenopiles-MAST","Stramenopiles-Chrysophytes","Stramenopiles-Other","Opisthokont-Fungi","Opisthokont-Metazoa","Opisthokonts-Other","Other/unknown","Unassigned")
tax_color=c("#800026","#cb181d","#e7298a","#df65b0","#fc4e2a","#fd8d3c","#fed976","#c7e9b4","#7fcdbb","#41ae76","#238b45","#006d2c","#00441b","#c6dbef","#6baed6","#1d91c0","#225ea8","#253494","#081d58","#54278f","#8c510a","#bf812d","#dfc27d")
names(tax_color)<-tax_order
data.agg$tax<-factor(data.agg$Taxa, levels=rev(tax_order)) #factoring
head(data.agg)

#Bar plot of community composition
ggplot(data.agg[order(data.agg$tax),], aes(y=x,fill=tax,x=Samples))+geom_bar(position = "fill", stat = "identity", color="black",aes(fill=tax))+scale_fill_manual(values=tax_color)+labs(title="", x="",y="Relative abundance of reads")+theme_bw()+theme(legend.position="right",axis.text.x = element_text(angle=45, hjust=1,vjust=1,color="black"))
#Change position = "fill" to position="stack" to look at total abundance
