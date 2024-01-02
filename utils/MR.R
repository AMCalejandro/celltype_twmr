#!/usr/bin/env Rscript


cmd_args=commandArgs(trailingOnly = TRUE)
# cmd_args = c("NULL", "ENSG00000213160", 
#                "/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/workdir/eqtlgen_height/chr2/ENSG00000213160/ENSG00000213160.matrix", 
#                "/home/rstudio/celltyping_wd/celltype_twmr/test_epqtl/workdir/eqtlgen_height/chr2/ENSG00000213160/ENSG00000213160.ld",
#              100000,
#              32000)

gene <- cmd_args[2]
Ngwas<-as.numeric(cmd_args[5])
N_eQTLs<-as.numeric(cmd_args[6])
out<-c("gene","alpha","SE","P","Nsnps","Ngene")


filecluster<-read.table(cmd_args[3],header=T,sep=",", dec = ".")


beta<-as.matrix(filecluster[,2:(length(filecluster[1,])-1)])
x<-colSums(abs(beta))
remove<-which(x==0)
if(length(remove)>0) {
  beta<-beta[,-remove]
}
beta<-as.matrix(beta)
gamma<-as.matrix(filecluster[,length(filecluster)])


C<-read.table(cmd_args[4],header=F,sep=" ",dec=".")
C<-as.matrix(C[,1:length(C[,1])])
 

S<-t(beta)%*%solve(C)%*%beta
H<-(1-1/sqrt(3781))*S+(1/sqrt(3781))*diag(length(S[,1]))
alpha<-solve(H)%*%(t(beta)%*%solve(C)%*%gamma)
alpha<-as.vector(alpha)
C_inv <- solve(C)
GCG_inv <- t(beta) %*% solve(C) %*% beta
GCG_inv<-(1-1/sqrt(3781))*GCG_inv+(1/sqrt(3781))*diag(length(GCG_inv[,1]))
GCG_inv<-solve(GCG_inv)


df_dg <- GCG_inv %*% t(beta) %*% C_inv
df_dG <- (GCG_inv %x% (t(gamma) %*% C_inv %*% ((beta %*% GCG_inv %*% t(beta)) %*% C_inv + diag(nrow(beta))))) + ((-t(gamma) %*% C_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% C_inv))
J <- cbind(df_dG, df_dg)


SEs<-c(rep(1/sqrt(32000),length(beta[1,])*length(beta[,1])),rep(1/sqrt(Ngwas),length(gamma[,1])))
R<-diag(length(beta[1,])+1)
Sigma <- (SEs %*% t(SEs)) * (C %x% R)
V <- J %*% Sigma %*% t(J)
se<- sqrt(V[1,1])


N=length(beta[,1])
Ngene=length(beta[1,])
Z<-alpha[1]/se
pval<-2*pnorm(abs(Z),lower.tail=FALSE)
line<-c(gene, alpha[1],se,pval,N,Ngene)
out<-rbind(out,line)


write.table(out,file=paste(tools::file_path_sans_ext(cmd_args[3]),".alpha",sep=""),quote=F,col.names=F,row.names=F)
