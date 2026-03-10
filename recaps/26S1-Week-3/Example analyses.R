#Load R functions
source('LDBayes.R')
source('LDcluster.R')

###########Load data
x <- as.matrix(read.table('Genotype_data.txt',head=FALSE))
map <- read.table('Genetic_map.txt',head=TRUE)
Pheno <- read.table('Phenotype_data.txt',head=TRUE)

###########SNP clustering analysis using LDna
library('LDna')
xl <- x

xl[xl<0.5] <- 0
xl[xl>=0.5&xl<=1.5] <- 1
xl[xl>=1.5] <- 2
source('LDcluster.R')

mapl <- cbind(as.numeric(as.factor(map[,2])),map[,3])

data_0.5_0.7 <- LDnClustering2(snp=xl, map=mapl, nSNPs = dim(xl)[2], w1 = 10, w2 = 100,
LD_threshold1 = 0.5, LD_threshold2 = 0.7, PC_threshold = 0.8,
mc.cores = 1, plot.network = NULL, threshold_net = 0.9)


group <- rep(0,dim(mapl)[1])

for (i in 1:length(data_0.5_0.7$clusters)){

group[data_0.5_0.7$clusters[[i]]] <- i  

}

##############LD Bayes 
y <- Pheno$Seed_index
n <- dim(x)[1]
p <- dim(x)[2]
ID <- c(1:n,1:n,1:n)
X <- matrix(0,nrow=3*n,ncol=4*p)
X[1:n,1:p] <- x
X[(n+1):(2*n),1:p] <- x
X[(2*n+1):(3*n),1:p] <- x
X[1:n,(p+1):(2*p)] <- x 
X[(n+1):(2*n),(2*p+1):(3*p)] <- x  
X[(2*n+1):(3*n),(3*p+1):(4*p)] <- x

W <- matrix(0,nrow=3*n,ncol=4)
W[1:(3*n),1] <- 1
W[1:n,2] <- 1
W[(n+1):(2*n),3] <- 1
W[(2*n+1):(3*n),4] <- 1


Group <- c(group,group,group,group)

thin <- 5
burnin <- 200
K <- 500

result <- LDBayes(x=X, y=y, ID, W, K, group=Group, thin, burnin) 
