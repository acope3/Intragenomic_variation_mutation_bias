#!/usr/bin/env Rscript

library(profmem)
library(argparse)
rm(list=ls())

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="FASTA file with protein-coding sequences",type="character",default="./")
parser$add_argument("-o","--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./")
parser$add_argument("-d","--div",help="Number of steps to diverge from starting values. Will be applied at beginning of each run, with the exception of the last.",type="integer",default=0)
parser$add_argument("-s","--samp",help="Number of samples",type="integer",default=5000)
parser$add_argument("-a","--adapt",help="Adaptive Width By Samples, i.e. will adapt every i samples",type="integer",default=50)
parser$add_argument("-t","--thin",help="Thinning value. Total number of iterations will be samples * thinning",type="integer",default=20)
parser$add_argument("-n","--threads",help="Number of threads to use for MCMC",type="integer",default=1)
parser$add_argument("--dEta",help="Initial dEta values. Assumes csv format with columns AA,Codon,DEta. First line should be a header.",type="character")
parser$add_argument("--dM",help="Initial dM values. Assumes csv format with columns AA,Codon,DM. First line should be a header.",type="character")
parser$add_argument("--phi",help="Initial Phi values. Assumes csv format with Gene IDs in first column and Phi values in second column", type="character")
parser$add_argument("--est_csp",help="Use this flag to indicate estimation of Phi. Otherwise, Phi will not be estimated",action="store_true")
parser$add_argument("--est_phi",help="Use this flag to indicate estimation of CSP. Otherwise, CSP will not be estimated",action="store_true")
parser$add_argument("--est_hyp",help="Use this flag to indicate estimation of Hyperparameters. Otherwise, Hyperparameters will not be estimated",action="store_true")
parser$add_argument("--est_mix",help="Use this flag to indicate estimation of mixture assignment. Otherwise, assignment will not be estimated",action="store_true")
parser$add_argument("--max_num_runs",help="Max number of runs to do.",type="integer",default = 6)
parser$add_argument("--fix_dEta",help="Use this flag to fix dEta at starting value.",action="store_true")
parser$add_argument("--fix_dM",help="Use this flag to fix dM at starting value",action="store_true")
parser$add_argument("--mixture_assignment",type="character",default=NULL)
parser$add_argument("--codon_table",type="integer",default=1)
parser$add_argument("--with_phi",action="store_true")
parser$add_argument("--obs_phi",type="character",default=NULL)
parser$add_argument("--restart_file",type="character",default=NULL)
parser$add_argument("--prior_from_gc",type="double",default=NULL)
parser$add_argument("--development",help="Run a developmental version of AnaCoDa",action="store_true")


args <- parser$parse_args()
print(args)

div <- args$div
input <- args$input
directory <- args$output
thinning <- args$thin
adaptiveWidth <- args$adapt
samples <- args$samp
num_threads <- args$threads
phi.files <- args$phi
dEta.file <- args$dEta
dM.file <- args$dM
est.phi <- args$est_phi
est.csp <- args$est_csp
est.hyp <- args$est_hyp
est.mix <- args$est_mix
max_num_runs <- args$max_num_runs
fix_dEta <- args$fix_dEta
fix_dM <- args$fix_dM
with.phi <- args$with_phi
obs.phi <- args$obs_phi
dev <- args$development
prior.from.gc <- args$prior_from_gc
mix.assign <- args$mixture_assignment
restart.file <- args$restart_file
codon.table <- args$codon_table

if(dev)
{
	library(AnaCoDa,lib.loc="~/R_dev/")
} else{
	library(AnaCoDa)
}

calcDeltaMFromGCBias <- function(gcBias, include.ref=TRUE)
{
  #' Calculates mutation bias terms $\Delta M$ based on GC content.
  #'
  #' @description Calculates $\Delta M$ based on GC bias.  Calculations
  #  assume freq A = freq T and freq C = freq G.
  #' @param gcBias: frequency of GC within a region or genome.
  #' @param include.ref: indicates whether reference codons should be included in results.
  #' Default value is TRUE
  #'
  #' @return Vector of $\Delta M$ values for each codon. Codon strings are used for entry names.
  #'
  #' @details Reference codons in AnaCoDa are the last alphabetical codon for each
  #' amino acid. By convention, the mutation and selection bias parameters, $\Delta M$ and $\Delta \eta$,
  #' are scaled so that $\Delta M = 0$ and $\Delta \eta = 0$ for each reference codon.
  #' In general the single codon amino acids and the stop codons (M, W, and X) are ignored.
#' These reference codons values are not used by `initializeParameterObject()` when setting
#' values of the `mutation.prior.mean` and `mutation.prior.sd`

  atBias  <-  1- gcBias
  
  ## calculate nt frequencies
  fG <- gcBias/2
  fC <- fG
  fA <- atBias/2
  fT <- fA
  
  fVec = c(fA,fC,fG,fT)
  names(fVec) <- c("A", "C", "G", "T")
  ## calculate M values
  ## note that sum fi = 1
  mVec <- - log(fVec)
  
  aaList  <- AnaCoDa::aminoAcids()
  ## Drop stop codons, M, and W
  aaList <- aaList[grep("X|M|W", aaList, invert = TRUE)]
  
  ## Create a vector for mutation terms
  M <- vector()
  deltaM <- vector()
  refCodon.list <- vector()
  for(aa in aaList){
    ## Get all codons, including reference.
    ## I would expect focal = TRUE to be the correct syntax for this, but it's not
    codons <- AnaCoDa::AAToCodon(aa, focal = FALSE);
    nCodons = length(codons);
    refCodon <- codons[nCodons]
    ##print(c(aa, nCodons, codons))
    refCodon.list[refCodon] <- refCodon
    ## Calculate M for each codon
    for(codon in codons)
    {
      thirdNt <- substring(codon, 3,3)
      M[codon] <- mVec[thirdNt]
    }
    
    ## Scale relative to last codon
    for(codon in codons)
    {
      deltaM[codon] = M[codon]-M[refCodon]
    }
    
  }
  
  if(include.ref == FALSE)
  {
  
    deltaM <- deltaM[which(!names(deltaM) %in% refCodon.list)]
  }
  return(deltaM)
}


## Outputs CSP estimates
createParameterOutput <- function(parameter,numMixtures,samples,mixture.labels,samples.percent.keep=1,relative.to.optimal.codon=F,report.original.ref=T)
{
  for (i in 1:numMixtures)
  {
    getCSPEstimates(parameter,paste(dir_name,"Parameter_est",mixture.labels[i],sep="/"),i,samples*samples.percent.keep,relative.to.optimal.codon=relative.to.optimal.codon,report.original.ref = report.original.ref)#,thin=thin)
  }
}


## Outputs traces for CSPs and plots expected frequecies as function of log10(\phi) 
createTracePlots <- function(trace, model,genome,numMixtures,samples,mixture.labels,samples.percent.keep=1)
{
  for (i in 1:numMixtures)
  {
    plot(trace, what = "Mutation", mixture = i)
    plot(trace, what = "Selection", mixture = i)

    plot(model, genome, samples = samples*samples.percent.keep, mixture = i,main = mixture.labels[i])
  }
  plot(trace, what="AcceptanceRatio")
}


if (with.phi && !is.null(obs.phi))
{
  genome <- initializeGenomeObject(file=input,codon_table=codon.table,match.expression.by.id = FALSE,observed.expression.file=obs.phi)
} else{
  if (dev)
  {
    genome <- initializeGenomeObject(file=input,codon_table=codon.table)
  } else {
    genome <- initializeGenomeObject(file=input)
  }
}



mixDef <- "mutationShared" ## options are allUnique, mutationShared, selectionShared, but this only matters if have multiple categories

percent.to.keep <- 0.5 ## Script is set up to stop adapting after 50% of samples and calculate posterior estimates, convergence, mean acceptance rate, etc. from last 50% of samples.
size <- length(genome)
index <- c(1:size)

if (!is.null(mix.assign))
{
  tmp <- read.csv(mix.assign,sep="\t",header=T,stringsAsFactors=F)
  geneAssignment <- tmp[,2]
  numMixtures <- length(unique(tmp[,2]))
  mixture.labels <- c("1","2")
} else{
  geneAssignment <- rep(1,size)
  mixture.labels <- unlist(strsplit(input,"/"))
  mixture.labels <- mixture.labels[length(mixture.labels)]
  numMixtures <- 1
}

init_phi <- NULL
sphi_init <- rep(1.5,numMixtures)
if (!is.null(phi.files))
{
  segment_exp <- read.table(file=phi.files,sep=",",header=TRUE)
  init_phi <- c(init_phi,segment_exp[,2])
  sphi_init <- rep(sd(log(init_phi)),numMixtures)
  if(length(genome) != length(init_phi))
  {
    stop("length(genomeObj) != length(init_phi), but it should.")
  } else{
  print("Initial Phi values successfully files loaded:");
  }
} 

if (!is.null(obs.phi))
{
  obs.phi <- read.csv(obs.phi,header=T,row.names=1)
  n.obs.phi <- ncol(obs.phi)

  s_eps <- rep(0.1,n.obs.phi)
} else {
  s_eps <- 0.1
}

if (!is.null(prior.from.gc))
{
  mutation.prior.mean <- calcDeltaMFromGCBias(prior.from.gc,F)
  mutation.prior.mean <- rep(mutation.prior.mean,2)
} else {
  mutation.prior.mean <- 0
}

dir.create(directory)
done <- FALSE
done.adapt <- FALSE
run_number <- 1

while((!done) && (run_number <= max_num_runs))
{
  if (run_number == 1)
  {
    div_run = 0
  } else{
    div_run = div
  }
  percent.to.keep <- 0.75
  if (is.null(restart.file))
  {
      parameter <- initializeParameterObject(genome,sphi_init,numMixtures, geneAssignment,init.sepsilon = s_eps,split.serine = TRUE, mixture.definition = mixDef, initial.expression.values = init_phi,init.w.obs.phi=with.phi,mutation.prior.mean=mutation.prior.mean)
      if (length(dM.file) > 0)
      {
        parameter$initMutationCategories(dM.file,1,fix_dM)
        
      } 
      if (length(dEta.file) > 0)
      {
        parameter$initSelectionCategories(dEta.file,1,fix_dEta)
        
      }
      
  } else {
      previous <- stringr::str_extract(pattern="restart_[0-9]+",string=restart.file)
      run_number <- as.numeric(stringr::str_extract(pattern="[0-9]+",string=previous)) + 1
      parameter<-initializeParameterObject(init.with.restart.file = restart.file,model="ROC")
  }

  steps.to.adapt <- (samples*thinning)*(1-percent.to.keep)
  dir_name <- paste0(directory,"/run_",run_number)
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))

  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                               est.expression=est.phi, est.csp=est.csp, est.hyper=est.hyp,est.mix=est.mix)
 
  mcmc$setStepsToAdapt(steps.to.adapt)

 
  model <- initializeModelObject(parameter, "ROC", with.phi,fix.observation.noise=F)
  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
  

  sys.runtime <- system.time(
    runMCMC(mcmc, genome, model, num_threads,div=div_run)
  )

  ## Output runtime for mcmc
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)


  ## Creates R objects, which can be later loaded for re-analzying already completed runs
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
 

  ## Output CSP file
  createParameterOutput(parameter = parameter,numMixtures = numMixtures,samples = samples,mixture.labels = mixture.labels,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T)#,thin=thinning)

  ## Output phi file
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep,genome = genome)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)
  
  #plots different aspects of trace
  trace <- parameter$getTraceObject()
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"))
  plot(mcmc,what = "LogPosterior")
  plot(mcmc,what="LogLikelihood")
  plot(trace, what = "ExpectedPhi")
  if (est.hyp)
  {
    plot(trace,what="Sphi")
  }
  if (with.phi)
  {
    plot(trace,what="Sepsilon")
  }
  if (est.csp)
  {
    ## Calculate auto-correlation and convergence of CSP traces
    param.conv <- TRUE
    if (!fix_dEta)
    {
      acfCSP(parameter,csp="Selection",numMixtures = numMixtures,samples=samples*percent.to.keep)
      for (i in 1:numMixtures)
      {
        param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Selection",mixture=i,frac1=0.25,frac2=0.5)
        z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
        if (length(z.scores) > 5)
        {
          param.conv <- FALSE
        }
        write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_eta_",i,".txt"),ncolumns = 1)
      }
    }
    if (!fix_dM)
    {
      acfCSP(parameter,csp="Mutation",numMixtures = numMixtures,samples=samples*percent.to.keep)
      for (i in 1:numMixtures)
      {
        param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Mutation",mixture=i,frac1=0.25,frac2=0.5)
        z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
        if (length(z.scores) > 5)
        {
          param.conv <- FALSE
        }
        write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_delta_M_",i,".txt"),ncolumns = 1)
      }
    }
  }
  dev.off()
  
  pdf(paste(dir_name,"Graphs/CSP_traces_CUB_plot.pdf",sep="/"), width = 11, height = 12)
  createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numMixtures,samples=samples,samples.percent.keep = percent.to.keep,mixture.labels = mixture.labels)
  dev.off()
  
 
  diag <- convergence.test(mcmc,samples = samples*percent.to.keep,thin=thinning,frac1=0.1)
  z<-abs(diag$z)

  ## Can end if overall log(posterior) and CSP parameters have converged
  #done <- (z < 1.96) && param.conv
  done <- F
  rm(parameter)
  rm(trace)
  rm(model)
  run_number <- run_number + 1
}
