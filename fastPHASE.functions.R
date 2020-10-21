#####################################
##### R Functions for fastPHASE haplotype analysis
##### Daniele Filiault 10 Dec 2012
####################################

########################################################################################
####
####    Please note: I am happy to share this code, but...
####    PLEASE don't use results for publication/presentation without talking to me first.
####    I spent a lot of time developing this method and my analyses are not yet published.
####    Many thanks in advance!  -d  filiaultd@gmail.com
####
########################################################################################

################
### get .csv SNP data to fastPHASE format and run fastPHASE
### DLF 07Dec2012
#################

runFastphase <- function(snpfilein,fileout){
  snps.up <- read.csv(snpfilein,stringsAsFactors=FALSE)
  
  ## SNPs should be biallelic - check this
  snp.levels <- apply(snps.up,1,function(x){
    bp.up<-x[-c(1:2)]
    bp.factor <- summary(as.factor(bp.up))
    bp.factor.no.na <- length(bp.factor[names(bp.factor)%in%"NA's"==FALSE])})
  not.bi<- which(snp.levels!=2)
  if (length(not.bi)==0){print("All SNPs biallelic")} else {
    print(paste(length(not.bi),"non-biallelic SNPs detected and removed"))
    snps.up<-snps.up[-not.bi,]}

  ## NAs should be "?"
  snps.up <- apply(snps.up,2,as.character)
  snps.up <- apply(snps.up,2,function(x) {
    x[is.na(x)]<-"?"
    return(x)})
  
  ## Write to fastPHASE format
  write.csv(snps.up, "snps.for.fastphase.csv",quote=FALSE,row.names=FALSE)
  csvToFastphase("snps.for.fastphase.csv",fileout)


  ## Run fastPHASE.  These are the default settings.  You may want to change them...
  system(paste("fastphase -oFastphaseOutput -Pzp ",fileout,sep=""))
}


###########################
### from a .csv SNP input file to a fastPHASE input file
### used inside the runFastphase function
### DLF 30Jan2012
#############################

csvToFastphase <- function(filein,fileout){
 snps <- read.csv(filein)

 # reformat data 
 snps.positions <- snps[,2]
 snps <- snps[,-(1:2)]
 snps <-t(snps)
 rownames(snps)<-gsub("X","",rownames(snps))

 snp.no <- ncol(snps)
 acc.no <- nrow(snps)

 #write in input format for fastPHASE - see fastPHASE
 snps.positions <- c("P",as.character(snps.positions))

 zz <- file(fileout,"w")
 cat(acc.no,file=zz,sep="\n")
 cat(snp.no,file=zz,sep="\n")
 cat(snps.positions,file=zz)
 cat("\n",file=zz)
 for(up in c(1:acc.no)){
   acc.out <- do.call("paste",c(as.list(snps[up,]),sep=''))
   cat(paste("ID",rownames(snps)[up],sep=" "),set="\n",file=zz)
   cat(paste(acc.out,"\n",sep=""),file=zz)
   cat(paste(acc.out,"\n",sep=""),file=zz)
   print(up)
                   }
 close(zz)
}


#################################
#### Parse output from fastPHASE
#### DLF 10Dec2012
#################################

parseFastphase <- function(fastphaseSNPs){
  input.data <- read.csv(fastphaseSNPs)
  input.accs <- input.data[grep("ID",input.data[,1]),]
  input.accs <- as.character(input.accs)
  input.accs <- gsub("ID ","",input.accs)
  input.accs <- gsub(" ","",input.accs)

  input.pos <- input.data[2,]
  input.pos <- gsub("P ","",input.pos)
  input.pos <- strsplit(input.pos," ")
  input.pos <- unlist(input.pos)

  data.files<-dir()
  data.files <- data.files[grep("E-Cluster", data.files)]
  acc.nos <- as.numeric(sapply(data.files,function(x) {strsplit(x,"_")[[1]][3]}))

  tester <- read.table(data.files[1])
  haplotypes <- matrix(NA, ncol=length(data.files),nrow=nrow(tester))
  rm(tester)

  for (up in c(1:length(data.files))){
    data.up <- read.table(data.files[up])
    data.up<-as.matrix(data.up)
    haplotypes[,up] <- apply(data.up,1,which.max)}

  colnames(haplotypes)<-input.accs[acc.nos]
  rownames(haplotypes)<-input.pos

  haplotypes.t <- t(haplotypes)
  haplotypes.t <- haplotypes.t[order(rownames(haplotypes.t)),]
  return(haplotypes.t)
}

#################################
#### Plot output from fastPHASE
#### DLF 11Dec2012
#################################
#test variables for troubleshooting...
#haplotypes <- haplotypes
#pheno.file <- "Phenotype_Take_sdXsdYd5_meidan_NPA.csv"
#gene.start <- 21254911
#gene.end <- 21257618
#meta.file <- "/Volumes/daniele.filiault/Documents/001.GWAS.reference.files/call_method_75/call_method_75/call_method_75_info.csv"
#plot.pdf.file <- "Take.plot.out.pdf"



plotFastphase <- function(haplotypes,pheno.file,gene.start,gene.end,meta.file,plot.pdf.file){
  
  ### step one: group the accessions by DNA distance for haplotype plotting
  ### and put haplotype information in that order.
  ### does this by distance across whole region.  Could subset for coding region or gene or any other sub-region...
  snps <- read.csv("snps.for.fastphase.csv",na.strings = c("?"))
  snps.pos <- snps[,2]
  snps <- snps[,-(1:2)]
  snps <- t(snps)
  rownames(snps)<-gsub("X","",rownames(snps))
  snps<-as.data.frame(snps)
  snps<-snps[order(rownames(snps)),]
  for (up in c(1:ncol(snps))){
    snps[,up]<-as.numeric(snps[,up])}
  snps <- snps-1
  # here is where you could subset snps to cluster using different SNP sets:
  snps.dist <- dist(snps)
  snps.clust <-hclust(snps.dist)
  snps.clust.dend <- as.dendrogram(snps.clust)
  haplotypes.snp.sort <- haplotypes[snps.clust$order,]

  ### step two: get latitude, longitude, phenotype info for accessions and put in appropriate order
  phenos <- read.csv(pheno.file,stringsAsFactors=FALSE)
  pheno.name <- colnames(phenos)[2]
  colnames(phenos)[1]<-"ecotype_id"
  meta <- read.csv(meta.file, stringsAsFactors=FALSE)
  all.info<-as.data.frame(as.matrix(rownames(snps)))
  colnames(all.info)<-"ecotype_id"
  all.info<- merge(all.info,meta,all.x=TRUE)
  all.info <- merge(all.info,phenos,all.x=TRUE)
  all.info.sort <-all.info[snps.clust$order,]

  ### step three: do actual plotting
  col.up.fig <- colors()[c(33,91,125,451,103,29,94,148,477,640,47,459,68,258,26,542,81,562,525,610)]
  # these colors chosen as maximally visible set of 20.  fastPHASE usually picks 20 groups...
  
  pdf(file=plot.pdf.file,width=23,height=14)
  
  # element 1. longitude
  par(mai=c(0.1,0.8,0.1,0.1),fig=c(0.0,1,0,0.08))
  plot(all.info.sort$longitude,col="purple",xaxt="n",ylab="longitude",xaxs="i")
  
  # element 2. latitude
  par(mai=c(0.1,0.8,0.1,0.1),fig=c(0.0,1,0.08,0.16),new=TRUE)
  #plot(hout.acc.info$LATITUDE) # one of these is wrong, at -37!
  #Need to ignore this, also ignoring two more with lat <30 to make plot more clear
  plot(all.info.sort$latitude,col="blue",xaxt="n",ylab="latitude",xaxs="i", ylim=c(30,69))
  
  # element 3. haplotypes
  par(mai=c(0.1,0.8,0.1,0.1),fig=c(0.0,1,0.16,0.6),new=TRUE)
  x <- 1:nrow(haplotypes.snp.sort)
  y <- 1:ncol(haplotypes.snp.sort)
  image(x,y,as.matrix(haplotypes.snp.sort),axes=FALSE,xlab="",ylab="region of interest",col=col.up.fig,xaxs="i")
  interest.gene <- c(gene.start:gene.end)
  interest.snps <- which(as.numeric(colnames(haplotypes.snp.sort))%in%interest.gene)
  if (length(interest.snps)>0){abline(h=c(min(interest.snps),max(interest.snps)))} else {print("Warning: No SNPs in gene of interest. No gene region marked on haplotype plot.")}

  #element 4. phenotypes
  par(mai=c(0.1,0.8,0.1,0.1),fig=c(0,1,0.6,0.80),new=TRUE)
  plot(all.info.sort[,colnames(all.info.sort)==pheno.name],xaxt="n",ylab=pheno.name,xaxs="i",yaxs="i",xlab="")
  
  #5.tree
  par(mai=c(0.1,0.8,0.1,0.1),fig=c(0,1,0.80,1),new=TRUE)
  plot(snps.clust.dend,leaflab="none",ann=FALSE,yaxt="n",xlim=c(0,max(x)),xaxs="i",yaxs="i")

  dev.off()
}
