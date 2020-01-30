# ECCO
ECCO is implemented as an open source R package for determining the optimal number of PEER factors in eQTL mapping studies. Instead of performing repetitive eQTL mapping, ECCO jointly applies differential expression analysis and Mendelian randomization (MR) analysis, leading to substantial computational savings. 

# Installation
ECCO is implemented as an R package, which can be installed from GitHub.

####  Install from GitHub
```
library(devtools)
install_github("3211895/ECCO")
```
# Usage
The main functions are ecco and ecco0. You can find the instructions and an example by '?ecco' and '?ecco0'.

##Example
```
data(exampledata)
attach(exampledata)
ind=1
genename=gene_name[ind]
gene=M_matrix[,ind]
geno=snp_raw[[ind]]
ivsnp=ecco0(gene,genename,gene_name,geno,ind)

num_peer=1
summary<-ecco(pheno,peer[[num_peer]],gene_name,iv_snp,num_peer)

A toy example for testing purposes only:
```
data(exampledata)
attach(exampledata)
N=length(gene_name)
iv_snp=c()
for(ind in 1:N)
{
tryCatch({
gene=M_matrix[,ind]
geno=snp_raw[[ind]]
genename=gene_name[ind]
ivsnp=ecco0(gene,genename,gene_name,geno,ind)
iv_snp=rbind(iv_snp,ivsnp)
},
error=function(e){})
}
res=c()
for(num_peer in 1:length(peerlist))
{
tryCatch({
pheno=Y
gene=M_matrix
geno=snp_raw
gene_name=gene_name
peerlist=c(1,2,5)
summary<-ecco(pheno,peer[[num_peer]],gene_name,iv_snp,peerlist[num_peer])
},
error=function(e){})
summary_total=rbind(summary_total,summary)
res=rbind(res,c(cor(as.numeric(summary[,4]),as.numeric(summary[,5])),peerlist[num_peer]))
}


## Our group

 <http://www.xzlab.org>.
