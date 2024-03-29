# install.packages("glasso")
# install.packages("igraph")
# install.packages("spacejam")
# install.packages("minet")#
# install.packages("WGCNA")
# install.packages("huge")
# install.packages("bnlearn")
# install.packages("pcalg")

library(glasso)
library(igraph)
library(spacejam)
library(minet)
library(WGCNA)
library(huge)

library(bnlearn)
library(pcalg)

#install.packages("glasso")

#install.packages("Bioconductor")
#library(BiocManager)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("WGCNA")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pcalg")

rm(list=ls())

##########################
## Undirected graphical models
##########################

##
## Read the Sachs et al data
##
#setwd('/Users/ashojaie/Dropbox/Teaching/shortcourse_NET/Ali/')

sachscov <- as.matrix(read.table("Network_omics/sachscov.txt"))
dim(sachscov)

sachscor <- cov2cor(sachscov)

sachsdat <- as.matrix(read.table("Network_omics/sachs.data.txt"))
dim(sachsdat)

p <- ncol(sachsdat)
n <- nrow(sachsdat)


##
## coexpression network
##
tau <- 0.5
A1 <- 1 * (round(abs(sachscor),3) > tau); diag(A1) <- 0
sum(A1)/2

##
## WGCNA
##
A2 <- round(adjacency(sachsdat), 3); diag(A2) <- 0

## co-expression network

## method 1: simple thresholding of correlations, at a cutoff chosen to give similar number of 
##			 edges to partial correlation methods

## a randomly chosen threshold, to get 50 nonzero values in the adjacency matrix (25 edges)
tau <- 0.22
A2 <- abs(sachscor) > tau; diag(A2) <- 0
sum(A2)

## method 2: testing for nonzero correlations

## testing for nonzero correlation, using Fisher Z-transform
fisherzs <- atanh(sachscor)
fisherps  <- 2*pnorm(abs(fisherzs), 0, 1/sqrt(n-3), lower.tail=FALSE)
A3 <- fisherps < (0.01/(p*(p-1)/2)); diag(A3) <- 0
sum(A3)


mim <- build.mim(sachsdat)
dim(mim)
A4<-round(minet::aracne(mim, eps=0.3), 3) # Need to specify threshold
A4<-1*(A4>0); diag(A4)<-0
## 
## plot the three networks
##
g1 <- graph.adjacency(A1, mode="undirected")
g2 <- graph.adjacency(A2, mode="undirected")
g3 <- graph.adjacency(A3, mode="undirected")
g4 <- graph.adjacency(A4, mode="undirected")


pdf('Network_omics/plot.pdf', width=9, height=3)
par(mfrow = c(1,4))
plot(g1,layout=layout.circle(g1), main='A1')
plot(g2,layout=layout.circle(g2), main='A2')
plot(g3,layout=layout.circle(g3), main='A3')
plot(g4,layout=layout.circle(g4), main='A4')

dev.off()
g0 <- g2

##
## ARACNE
## mutual information matrix
mim <- build.mim(sachsdat)
dim(mim)
A4<-round(minet::aracne(mim, eps=0.3), 3) # Need to specify threshold

##
## Partial correlation networks
##
abs(round(sachscor,3))
abs(round(solve(sachscor),3))

invcov <- abs(round(solve(sachscor),3))
1*(invcov > 0.5)

## calculate lambda, based on formula in the slides
alpha = 0.01
num = qt(p=alpha/(2*(p^2)),df=n-2, lower.tail=F)
lambda = num / sqrt(n-2 + num)

## glasso
glasso.est <- glasso(s=sachscor, rho=lambda, approx=FALSE, penalize.diagonal=FALSE)
A1 <- abs(glasso.est$wi) > 1E-16; diag(A1) <- 0
g1 <- graph.adjacency(A1, mode="undirected")

## neighborhood selection
ns.est <- glasso(s=sachscor, rho=lambda, approx=TRUE, penalize.diagonal=FALSE)
A2 <- abs(ns.est$wi) > 1E-16; diag(A2) <- 0
g2 <- graph.adjacency(A2, mode="undirected")

## nonparanormal
scor <- cor(sachsdat,method='spearman')
scor <- 2*sin(scor*pi/6)
npn.est <- glasso(s=scor, rho=lambda, approx=FALSE, penalize.diagonal=FALSE)
A3 <- abs(npn.est$wi) > 1E-16; diag(A3) <- 0
g3 <- graph.adjacency(A3, mode="undirected")

## nonparanormal -- alternative estiamtion 
library(huge)

npn.cor <- huge.npn(x=sachsdat, npn.func="skeptic", npn.thresh=NULL, verbose=FALSE)
npn.est <- glasso(s=npn.cor, rho=lambda, penalize.diagonal=FALSE)
A4 <- abs(npn.est$wi) > 1E-16; diag(A4) <- 0
g4 <- graph.adjacency(A4, mode="undirected")

## spacejam
library(spacejam)
spacejam.est <- SJ(sachsdat, lambda=0.5)
A5 <- 1*(spacejam.est$G)[,,1]
g5 <- graph.adjacency(A5, mode="undirected")


##
## binary network estimation
##

library(glmnet)

head(sachsdat)

sachsbin <- 1*(sachsdat > 0) + -1*(sachsdat <= 0)
head(sachsbin)

bin.est <- matrix(0,p,p)
## estiamte the neighborhood for each node 
for(j in 1:p){
	## this is the same method used in neighborhood selection, the only difference is 'family'
	nbr <- glmnet(x=sachsbin[,-j], y=sachsbin[,j], family='binomial', lambda=lambda) 
	bin.est[j,-j] <- 1*(abs(as(nbr$beta,"matrix")) > 0)	#store the estimates in jth row of matrix
}
A6 <- bin.est; diag(A6) <- 0
sum(A6)
g6 <- graph.adjacency(A6, mode="undirected")


## 
## plot the networks
##
pdf('plot.pdf', width=9, height=6)
par(mfrow = c(2,3), mar=c(1,1,4,1))
#plot(g0,layout=layout.circle(g0), main='co-expression')
plot(g1,layout=layout.circle(g1), main='glasso')
plot(g2,layout=layout.circle(g2), main='NS')
plot(g3,layout=layout.circle(g3), main='nonparanormal')
plot(g4,layout=layout.circle(g4), main='nonparanormal - v2')
plot(g5,layout=layout.circle(g5), main='spacejam')
plot(g6,layout=layout.circle(g3), main='Binary')
dev.off()


##########################
## Bayeisan networks
##########################
#install.packages("pcalg")
#biocLite("graph")
#biocLite("RBGL")
#biocLite("Rgraphviz")
library(pcalg)

dat <- read.table('sachs.data')
p <- ncol(dat)
n <- nrow(dat)

#inf <- read.table('sachs.info')
ps <- c("praf","pmek","plcg","PIP2","PIP3","P44","pakts","PKA","PKC","P38","pjnk")
colnames(dat) <- ps

X2 <- as.data.frame(scale(dat))


library(bnlearn)
## Grow-Shrink
A3 <- gs(x=X2, alpha=0.001)

## Hill climbing
A4 <- hc(X2)

compare(A3, A4)

## plot the graphs
pdf('plot.pdf', width=10, height=5)
par(mfrow=c(1,2), mar=c(1,1,4,1))
plot(A3, main='Grow-Shrink')
plot(A4, main='Hill Climbing')
dev.off()

## pcalg
indepTest <- gaussCItest

## define sufficient statistics
suffStat <- list(C=cor(dat), n=n)

## estimate CPDAG
pc.fit <- pc(suffStat=suffStat, indepTest=indepTest, labels=ps, alpha=0.1, verbose=FALSE)
#plot(pc.fit, main='PC Algorithm')

A5 <- pc.fit@graph		#get a graphNEL obj
A5 <- as(A5, "matrix")	#get the adjacency matrix for graphNEL obj
g5 <- graph.adjacency(A5, mode='directed')	#create an igraph obj
V(g5)$name <- ps		#assign vertex names

##get edge modes
##NOTE: the adjmat from pcalg is in "col"(biology) format, but
##igraph reads them in "row"(cs) format, need to transpose!!
getemode4CPDAG <- function(amat, est){    
    if(est == "pcalg") amat = t(amat)
    
    emode <- 2*amat
    emode[(amat==t(amat)) & (amat!=0)] <- -2
    emode <- emode[amat!=0]
    emode[emode==-2] <- 0
    
    return(emode)
}

par(mar=c(0,0,1,0))
plot(g5, layout=layout.circle(g5), vertex.size=25, vertex.color=NA,
	main='PC Algorithm', edge.arrow.mode=getemode4CPDAG(A5,"pcalg"))

##
##drawing all estimates together
##
A3 <- amat(A3)	#GS
g3 <- graph.adjacency(A3, mode='directed')	#create an igraph obj
V(g3)$name <- ps		#assign vertex names

A4 <- amat(A4)	#HC (scaled data)
g4 <- graph.adjacency(A4, mode='directed')	#create an igraph obj
V(g4)$name <- ps		#assign vertex names

pdf('plot.pdf', width=15, height=5)
par(mfrow=c(1,3), mar=c(0,0,1,0))
plot(g5, layout=layout.circle(g5), vertex.size=25, vertex.color=NA,
main='PC Algorithm', edge.arrow.mode=getemode4CPDAG(A5,"pcalg"))
plot(g3, layout=layout.circle(g3), vertex.size=25, vertex.color=NA,
main='Grow-Shrink', edge.arrow.mode=getemode4CPDAG(A3,"pcalg"))
plot(g4, layout=layout.circle(g4), vertex.size=25, vertex.color=NA,
main='Hill Climbing', edge.arrow.mode=getemode4CPDAG(A4,"pcalg"))
dev.off()

##
## Penalize likelihood estimation of DAGs
##
library(spacejam)
##NOTE: SJ.dag requires that columns of data matrix are ordered according to causal ordering.
##		Here, it is simply assumed that the variables in "sachsdat" are correctly ordered,
##		which is likely wrong!!
causalord <- c(1:ncol(sachsdat))
spacejam.DAG <- SJ.dag(sachsdat[,causalord], lambda=0.5)
A5 <- 1*(spacejam.est$G); diag(Ag) <- 0
g5 <- graph.adjacency(A5, mode="directed")



##########################
## Net based pathway analysis
##########################

##
## SPIA
##

rm(list=ls())

library(SPIA)

data(colorectalcancer)

res <- spia(de=DE_Colorectal, all=ALL_Colorectal, organism="hsa", beta=NULL,
nB=2000, plots=FALSE, verbose=FALSE, combine="fisher")

res$pG=combfunc(res$pNDE,res$pPERT,combine="norminv")
res$pGFdr=p.adjust(res$pG,"fdr")
res$pGFWER=p.adjust(res$pG,"bonferroni")
plotP(res,threshold=0.05)

points(I(-log(pPERT))~I(-log(pNDE)),data=res[res$ID=="05210",],col="green",
pch=19,cex=1.5)

install.packages("netgsa")

##
## NetGSA (see the NetGSA example)
##




##
## Additional code -- plots in slides etc
##

pdf('/Users/ashojaie/Dropbox/Teaching/shortcourse_NET/Ali/slides/Figures/cortest-1.pdf', 4, 3.5)
par(mar=c(2,2,1,1))
set.seed(1)
xseq <- seq(-2,2,.001)
n <- 10
plot(xseq, dnorm(x=xseq, 0, 1/sqrt(n-3)), type='l', main='', xlab='', ylab='', ylim=c(0,1.7))
abline(h=0)
segments(x0=0, x1=0, y0=0, y1=dnorm(x=0, 0, 1/sqrt(n-3)), lty=3)
lb <- round(-1.96 * (1/sqrt(n-3)), 2)
segments(x0=lb, x1=lb, y0=0, y1=dnorm(x=lb, 0, 1/sqrt(n-3)), lty=1, lwd=2, col=2)
ub <- round(1.96 * (1/sqrt(n-3)), 2)
segments(x0=ub, x1=ub, y0=0, y1=dnorm(x=ub, 0, 1/sqrt(n-3)), lty=1, lwd=2, col=2)
text(0.8,0.5,'n=10',col=1)
text(-1.3,1.5,expression(N(0,1/(n-3))),col=1)
dev.off()

pdf('/Users/ashojaie/Dropbox/Teaching/shortcourse_NET/Ali/slides/Figures/cortest-2.pdf', 4, 3.5)
par(mar=c(2,2,1,1))
set.seed(1)
xseq <- seq(-2,2,.001)
n <- 10
plot(xseq, dnorm(x=xseq, 0, 1/sqrt(n-3)), type='l', main='', xlab='', ylab='', ylim=c(0,1.7))
text(-1.3,1.5,expression(N(0,1/(n-3))),col=1)
abline(h=0)
segments(x0=0, x1=0, y0=0, y1=dnorm(x=0, 0, 1/sqrt(n-3)), lty=3)
lb <- round(-1.96 * (1/sqrt(n-3)), 2)
segments(x0=lb, x1=lb, y0=0, y1=dnorm(x=lb, 0, 1/sqrt(n-3)), lty=1, lwd=2, col=2)
ub <- round(1.96 * (1/sqrt(n-3)), 2)
segments(x0=ub, x1=ub, y0=0, y1=dnorm(x=ub, 0, 1/sqrt(n-3)), lty=1, lwd=2, col=2)
text(0.8,0.5,'n=10',col=1)

n <- 20
lines(xseq, dnorm(x=xseq, 0, 1/sqrt(n-3)), type='l', col='navy')
lb <- round(-1.96 * (1/sqrt(n-3)), 2)
segments(x0=lb, x1=lb, y0=0, y1=dnorm(x=lb, 0, 1/sqrt(n-3)), lty=1, lwd=2, col=4)
ub <- round(1.96 * (1/sqrt(n-3)), 2)
segments(x0=ub, x1=ub, y0=0, y1=dnorm(x=ub, 0, 1/sqrt(n-3)), lty=1, lwd=2, col=4)
text(0.6,1.2,'n=20',col='navy')
dev.off()



