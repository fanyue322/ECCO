# ECCO
ECCO is implemented as an open source R package for determining the optimal number of PEER factors in eQTL mapping studies. Instead of performing repetitive eQTL mapping, ECCO jointly applies differential expression analysis and Mendelian randomization (MR) analysis, leading to substantial computational savings. 

# Installation
ECCO is implemented as an R package, which can be installed from GitHub.

####  Install from GitHub
```
library(devtools)
install_github("fanyue322/ECCO")
```
# Usage
The main functions are ecco and ecco0. You can find the instructions and an example by '?ecco' and '?ecco0'.

# Example
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
```
### Running the pipeline：

#### 1）Calculate the gene expression residuals with PEER package
```
# Input: the gene expression data dt, an N*P matrix; peer, the number of PEER factors to be removed
  model = PEER()
  PEER_setPhenoMean(model,as.matrix(dt)) 
  dim(PEER_getPhenoMean(model))
  PEER_setAdd_mean(model, TRUE)
  PEER_setNk(model,peer)   
  PEER_getNk(model)
  PEER_update(model)
  factors = PEER_getX(model)
  factors=factors[,-1]
  residuals = PEER_getResiduals(model)
  write.table(residuals, paste(tissue,'_peer', pc, ".txt", sep=""), quote=F, col.names=F, row.names=F)

```
#### 2）Select the cis-SNP with ecco0
```

```


#### 3）Select the cis-SNP with ecco0
```
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
res=data.frame(res)
optimal_num_peer=res[which(res[,1]==max(res[,1])),2]
```


# Results reproduced
All results from all methods used in the ECCO paper can be reproduced at 
 <https://github.com/fanyue322/ECCOreproduce>.

## Our group

 <http://www.xzlab.org>.
