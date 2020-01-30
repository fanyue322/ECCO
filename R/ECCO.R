########################################################################################################################
# Package: ECCO
# Version: 1.0
# Date   : 2020-01-20
# Title  : Controlling for Confounding Effects in eQTL Mapping Studies through Joint Differential Expression and Mendelian Randomization Analyses
# Authors: Yue Fan and Xiang Zhou
# Contact: xafanyue@163.com
#          University of Michigan, Department of Biostatistics
########################################################################################################################

#' Controlling for Confounding Effects in eQTL Mapping Studies through Joint Differential Expression and Mendelian Randomization Analyses
#'
#' Determining the optimal number of PEER factors for eQTL mapping analysis through DE and MR analysis
#'
#' Instead of performing repetitive eQTL mapping, ECCO jointly applies differential expression analysis and Mendelian randomization (MR) analysis, leading to substantial computational savings.
#'
#' @param gene the gene expression data.
#' @param genename the name of the gene
#' @param gene_name a vector containg the names of all genes
#' @param geno a matrix containg all the cis-SNPs of the analyzed gene
#' @param ind the index of the gene
#' @return \item{iv_snp}{the snp that has the strongest association with the gene}
#' @author Yue Fan, Shiquan Sun, Xiang Zhou
#' @examples
#' data(exampledata)
#' attach(exampledata)
#' ind=1
#' genename=gene_name[ind]
#' gene=M_matrix[,ind]
#' geno=snp_raw[[ind]]
#' ivsnp=ecco0(gene,genename,gene_name,geno,ind)
#' closeAllConnections()
#' detach(exampledata)
ecco0 <- function(gene,genename,gene_name,geno,ind) {


  M <-gene


  snps = SlicedData$new();
  A=as.matrix(geno)
  snps$CreateFromMatrix(A)

  gene=SlicedData$new();
  M=as.matrix(M)
  M=t(M)
  gene$CreateFromMatrix(M)

  useModel = modelLINEAR;
  pvOutputThreshold = 9.99e-1;
  errorCovariance = numeric();
  output_file_name = tempfile();

  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = FALSE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);

  me=me$all$eqtls[1,]
  tmp=as.character(me$snps)
  onesnp=as.numeric(substr(tmp,4,nchar(tmp)))
  idx=which(gene_name==genename)
  iv_snp=c(A[onesnp,],idx)

  return(iv_snp)
}



#' Controlling for Confounding Effects in eQTL Mapping Studies through Joint Differential Expression and Mendelian Randomization Analyses
#'
#' Determining the optimal number of PEER factors for eQTL mapping analysis through DE and MR analysis
#'
#' Instead of performing repetitive eQTL mapping, ECCO jointly applies differential expression analysis and Mendelian randomization (MR) analysis, leading to substantial computational savings.
#'
#' @param pheno the phenotype data
#' @param gene a matrix containing the whole gene expression data after removing the PEER factors
#' @param gene_name a vector containg the names of all genes
#' @param iv_snp the instrumental variable, produced from ecco0
#' @param peer the number of peer factors to be examined
#' @return \item{Gene}{the name of the gene}
#' @return \item{PEER}{the number of PEER factors}
#' @return \item{p-value}{the p-value of alpha}
#' @return \item{beta_hat}{the estimation of beta from MR model}
#' @return \item{beta_tilde}{the estimation of beta from DE model}
#' @examples
#' data(exampledata)
#' attach(exampledata)
#' num_peer=1
#'  summary<-ecco(pheno,peer[[num_peer]],gene_name,iv_snp,num_peer)
#' closeAllConnections()
#' detach(exampledata)

ecco <- function(pheno,gene,gene_name,iv_snp,peer) {
  summary=c()
  N=nrow(iv_snp)
  samplesize=nrow(gene)



  for(ind in 1:N)
  {
    M=gene[,ind]
    ivsnp=iv_snp[ind,]
    genename=gene_name[ivsnp[samplesize+1]]

    ivsnp=ivsnp[1:samplesize]
    l1<-lm(M~ivsnp)
    l2 <- lm(Y ~ M)
    l3 <- lm(Y ~ ivsnp)
    summary <- rbind(summary, c(genename,peer, summary(l1)$coeff[2,4], summary(l3)$coeff[2,1]/summary(l1)$coeff[2,1], summary(l2)$coeff[2,1]))
    colnames(summary)=c('Gene','PEER','p-value','beta_hat','beta_tilde')
  }

  return(summary)
}
