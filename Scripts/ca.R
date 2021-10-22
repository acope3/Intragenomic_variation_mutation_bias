library(seqinr)
library(AnaCoDa)
library(ca)
library(cluster)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="FASTA file with protein-coding sequences",type="character")
parser$add_argument("-o","--output",help="Output directory",type="character",default="./")
parser$add_argument("-k","--num_clusters",help="Number of clusters to use in CLARA",type="integer",default=2)

args <- parser$parse_args()
input <- args$input
output <- args$output
k <- args$num_clusters


genome.file <- file.path("Genomes","Genomes","cds_cleaned",input)

genome <- initializeGenomeObject(genome.file)
size <- length(genome)

genes <- lapply(1:size,function(x){tmp<-genome$getGeneByIndex(x,F);tmp$seq})

genes <- lapply(genes,function(x){unlist(strsplit(x,split="",fixed=T))})

ids <- lapply(1:size,function(x){tmp<-genome$getGeneByIndex(x,F);tmp$id})
df <- data.frame(genes=unlist(ids))
df[words()] <- numeric(length=size)

## Use absolute frequencies
codon.counts<-lapply(genes,function(x){
			tmp <- uco(x,frame=0,index="eff");
			tmp <- tmp/(length(x)/3)
			})


for (i in 1:size)
{
	df[i,words()] <- codon.counts[[i]][words()]
}
df[is.na(df)] <-0

rownames(df) <- df$genes
df<-df[,2:ncol(df)]
df<-df[,which(!colnames(df) %in% c("taa","tga","tag"))]
fit <- ca(df)
print("Beginning clustering...")
clarax <- clara(x=fit$rowcoord[,1:4],k=k,samples=200,sampsize=length(genome)/2,rngR=T,pamLike=T)
print("Done...")

clusters <- data.frame(Gene=names(clarax$clustering),Cluster=clarax$clustering,stringsAsFactors=F)
write.table(clusters,file.path(output,input),sep="\t",col.names=T,row.names=F,quote=F)

pdf(file.path(output,paste0(input,"_",k,".pdf")))
plot(fit,labels=c(0,2))
plot(factoextra::fviz_cluster(clarax,labelsize=0))
dev.off()

