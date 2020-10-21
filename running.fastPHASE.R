##### fastPHASE haplotype analysis and plots
##### Daniele Filiault 07 Dec 2012

########################################################################################
####
####    Please note: I am happy to share this code, but...
####    PLEASE don't use results for publication/presentation without talking to me first.
####    I spent a lot of time developing this method and my analyses are not yet published.
####    Many thanks in advance!  -d  filiaultd@gmail.com
####
########################################################################################

####  Before you start, read the documentation "Haplotype plots with fastPHASE"
####  for more info about what you need to make this work.

### before you begin, set your working directory appropriately.  This is mine for testing...
setwd("/Volumes/daniele.filiault/Documents/006.FLC.evolution/001.FLC.comparative/001.sequences/001.DNA/005.peijin.fall2011/fastPHASE.functions")

### and source functions
source("fastPHASE.functions.R")

### one analysis per directory, please!  Otherwise fastPHASE will overwrite your old results...


#####################
####  STEP 1: Making sure SNPs are in the correct input format and running fastPHASE with the default settings.
####  To play with these, read the fastPHASE manual and edit the runFastphase function.
######################

runFastphase("SNPs.Take.callmethod52.filtered.csv","snps.fastphase.format")
# first argument is your .csv SNP file, second is the name you want to give the fastPHASE-formatted SNP file.



#########################
#### STEP 2: Parse fastPHASE results.
#########################

haplotypes<-parseFastphase("snps.fastphase.format")
# argument is the name of the fastPHASE-formatted SNP file generated in STEP 1.
# this extracts the MOST likely estimated haplotype group for each SNP for each accession.
# in case of weirdness, it could be good to look at actual numbers here.



#########################
#### STEP 3: Make plots
#########################
plotFastphase(haplotypes, "Phenotype_Take_sdXsdYd5_meidan_NPA.csv", 21254911, 21257618, "call_method_75_info.csv", "Take.plot.out.pdf")

## arguments are as follows:
#  haplotypes (haplotypes) - parsed fastPHASE output from step 2
#  pheno.file ("Phenotype_Take_sdXsdYd5_meidan_NPA.csv") - two column .csv.  First column is ecotype_id, second is phenotype values.  Column name of second column should be the phenotype name.
#  gene.start (21254911) - coordinate where the gene of interest begins.  Should be a lower number than gene.end
#  gene.end (21257618) - coordinate where the GOI ends.
#  meta.file (call_method_75_info.csv") - this is the file that describes the accessions used for call method 75.  I am using call method 75 accessions ONLY in the analysis because the "bad" arrays have been filtered out in call method 75 and not in call method 52.  
#  plot.pdf.file ("Take.plot.out.pdf") - name of the pdf file for the plot



#Take's gene At5g52350
