#' Example dataset
#'
#' A simulated example dataset of ECCO.
#' The variables are as follows:
#'
#' @docType data
#' @keywords datasets
#' @name exampledata
#' @usage data(exampledata)
#' @format Contains the following objects:
#' \describe{
#'   \item{gene_name}{a vector containg the names of all genes.}
#'   \item{M_matrix}{a matrix containing the whole gene expression data.}
#'   \item{peer}{a data list containing the whole gene expression data after removing the PEER factors.}
#'   \item{snp_raw}{a data list containing the whole genotype data.}
#'   \item{Y}{the phenotype data.}
#'   \item{peerlist}{the number of peer factors to be examined.}
#' }
#' @examples
#' data(exampledata)
#' attach(exampledata)
#' N=length(gene_name)
#' iv_snp=c()
#' for(ind in 1:N)
#' {
#' tryCatch({
#'  gene=M_matrix[,ind]
#'  geno=snp_raw[[ind]]
#'  genename=gene_name[ind]
#'  ivsnp=codec0(gene,genename,gene_name,geno,ind)
#'  iv_snp=rbind(iv_snp,ivsnp)
#'  },
#'  error=function(e){})
#'  }
#'  summary_total=c()
#'  for(num_peer in 1:length(peerlist))
#'  {
#'   tryCatch({}
#'   pheno=Y
#'   gene=M_matrix
#'   geno=snp_raw
#'   gene_name=gene_name
#'   peerlist=c(1,2,5)
#'   summary<-ecco(pheno,peer[[num_peer]],gene_name,iv_snp,peerlist[num_peer])
#'   },
#'   error=function(e){})
#'   summary_total=rbind(summary_total,summary)
#'   }
#' closeAllConnections()
#' detach(exampledata)
"exampledata"
