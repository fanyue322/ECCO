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
#' @param Y the phenotype data
#' @return \item{iv_snp}{the snp that has the strongest association with the gene}
#' @return \item{summary]{the effect size estiamted when the number of PEER is 0}
#' @author Yue Fan, Shiquan Sun, Xiang Zhou
#' @examples
#' data(exampledata)
#' attach(exampledata)
#' ind=1
#' genename=gene_name[ind]
#' gene=M_matrix[,ind]
#' geno=snp_raw[[ind]]
#' result=ecco0(gene,genename,gene_name,geno,ind,Y)
#' closeAllConnections()
#' detach(exampledata)
ecco0 <- function(gene,genename,gene_name,geno,ind,Y) {


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

  M=t(M)

  l1<-lm(M~A[onesnp,])
  l2 <- lm(Y ~ M)
  l3 <- lm(Y ~ A[onesnp,])
  summary <- data.frame(Gene=genename,PEER=0, p_value=summary(l1)$coeff[2,4], beta_tilde=summary(l3)$coeff[2,1]/summary(l1)$coeff[2,1], beta_hat=summary(l2)$coeff[2,1])

  results=list()
  results[[1]]<-iv_snp
  results[[2]]<-summary
  names(results)=c('iv_snp','summary')
  return(results)
}

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
#' @param Y the phenotype data
#' @param r2 clumping r2 cutoff, default is 0.1
#' @param kb clumping window, default is 1
#' @return \item{iv_snp}{the snp that has the strongest association with the gene}
#' @return \item{summary]{the effect size estiamted when the number of PEER is 0}
#' @author Yue Fan, Xiang Zhou
#' @examples
#' data(exampledata)
#' attach(exampledata)
#' ind=1
#' genename=gene_name[ind]
#' gene=M_matrix[,ind]
#' geno=snp_raw_ivw[[ind]]
#' result=ecco0_ivw(gene,genename,gene_name,geno,ind,Y)
#' closeAllConnections()
#' detach(exampledata)
ecco0_ivw <- function(gene,genename,gene_name,geno,ind,Y,r2=0.1,kb=1) {
  ivsnp=c()
  M <-gene
  samplesize=length(M)
  snp<-geno[,-((ncol(geno)-3):(ncol(geno)))]
  rownames(snp)=NULL
  snp_infor<-geno[,((ncol(geno)-3):(ncol(geno)))]
  snps = SlicedData$new();
  A=as.matrix(snp)
  snps$CreateFromMatrix(A)

  gene=SlicedData$new();
  M=as.matrix(M)
  M=t(M)
  gene$CreateFromMatrix(M)

  useModel = modelLINEAR;
  pvOutputThreshold = 1;
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

  me=me$all$eqtls
  tmp=as.character(me$snps)
  onesnp=as.numeric(substr(tmp,4,nchar(tmp)))
  chr=snp_infor[1,1]
  sig_snp_set1<-data.frame("chr_name"=chr,"chrom_start"=snp_infor[onesnp,2],SNP=snp_infor[onesnp,4],"pval.exposure"=me[,4])
  ind_sig_snp <- clump_data(sig_snp_set1,clump_r2 =r2,clump_kb=kb)
  if(nrow(ind_sig_snp)==1)
  {
    idx=which(sig_snp_set1[,2]==as.numeric(ind_sig_snp[2]))
    tmp=as.character(me$snps[idx])
    onesnp=as.numeric(substr(tmp,4,nchar(tmp)))
    ind=which(gene_name==genename)
    tmp=c(A[onesnp,],ind)
    iv_snp=rbind(iv_snp,tmp)
    M=t(M)
    l1<-lm(M~A[onesnp,])
    l2 <- lm(Y ~ M)
    l3 <- lm(Y ~ A[onesnp,])
    pve_g=var(A[onesnp,])*summary(l1)$coeff[2,1]^2/var(M)
    Fstat=pve_g*(samplesize-2)/(1-pve_g)
    summary <- data.frame(Gene=genename,PEER=0, Fstat=Fstat, beta_tilde=summary(l3)$coeff[2,1]/summary(l1)$coeff[2,1], beta_hat=summary(l2)$coeff[2,1])
  }else{
    idx=match(ind_sig_snp[,2],sig_snp_set1[,2])
    tmp=as.character(me$snps[idx])
    onesnp=as.numeric(substr(tmp,4,nchar(tmp)))
    if(length(onesnp)>5)
    {
      onesnp=onesnp[1:5]
    }
    ind=which(gene_name==genename)
    tmp=cbind(A[onesnp,],ind)
    iv_snp=rbind(iv_snp,tmp)
    num_snp=nrow(tmp)
    summary_statistic=c()
    for(i in 1:num_snp)
    {
      l1 <- lm(t(M) ~ A[onesnp[i],])
      l3 <- lm(Y ~  A[onesnp[i],])
      summary_statistic=rbind(summary_statistic,c(summary(l1)$coeff[2,1],summary(l3)$coeff[2,1],summary(l1)$coeff[2,2],summary(l3)$coeff[2,2]))
    }
    IVW=mr_ivw(summary_statistic[,1],summary_statistic[,2],summary_statistic[,3],summary_statistic[,4])
    l2 <- lm(Y ~ t(M))
    lf<-lm(t(M)~t(A[onesnp,]))
    pve_g=sum(summary(lf)$coeff[2:6,1]^2*apply(A[onesnp,],1,var))/var(M)
    Fstat=pve_g*(samplesize-num_snp-1)/num_snp*(1-pve_g)
    summary <- data.frame(Gene=genename,PEER=0, Fstat=Fstat, beta_tilde=IVW$b, beta_hat=summary(l2)$coeff[2,1])
  }
  results=list()
  results[[1]]=iv_snp
  results[[2]]=summary
  names(results)=c('iv_snp','summary')
  return(results)
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
#' @param summary the effect size estiamted when the number of PEER is 0, produced from ecco0 or ecc0_ivw
#' @return \item{Gene}{the name of the gene}
#' @return \item{PEER}{the number of PEER factors}
#' @return \item{p-value}{the p-value of alpha}
#' @return \item{beta_hat}{the estimation of beta from MR model}
#' @return \item{beta_tilde}{the estimation of beta from DE model}
#' @examples
#' data(exampledata)
#' attach(exampledata)
#' num_peer=1
#' summary<-ecco(Y,peer[[num_peer]],gene_name,iv_snp,num_peer,summary)
#' closeAllConnections()
#' detach(exampledata)

ecco <- function(pheno,gene,gene_name,iv_snp,peer,summary) {
  N=nrow(iv_snp)
  samplesize=nrow(gene)
  p_value=summary[3]
  beta_tilde=summary[4]
  for(ind in 1:N)
  {
    M=gene[,ind]
    ivsnp=iv_snp[ind,]
    genename=gene_name[ivsnp[samplesize+1]]

    ivsnp=ivsnp[1:samplesize]
    l2 <- lm(Y ~ M)
    summary=rbind(summary,data.frame(Gene=genename,PEER=peer,p_value=p_value,beta_tilde=beta_tilde,beta_hat=summary(l2)$coeff[2,1]))
  }

  return(summary)
}
