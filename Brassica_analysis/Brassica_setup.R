options(stringsAsFactors=F)
in_data = read.csv('rlog2_shade_brassica_shade.csv')
genotypes = read.csv('shade_gene_cis_genotypes_v1.csv')
cisGenotypes = t(as.matrix(genotypes[,-1]))#,dimnames=list(colnames(genotypes)[-1],genotypes[,1]))
rownames(cisGenotypes) = colnames(genotypes)[-1]
colnames(cisGenotypes) = genotypes[,1]

# replace NAs with 0.5
cisGenotypes[is.na(cisGenotypes)] = 0.5

data = data.frame(ID = colnames(in_data)[-c(1)])
data$Genotype = sapply(as.character(data$ID),function(x) paste(strsplit(x,'_')[[1]][1:2],collapse='_'))
data$TRT = sapply(as.character(data$ID),function(x) strsplit(x,'_')[[1]][3])

Y = t(in_data[,-1])
colnames(Y) = in_data[,1]
rownames(Y) = data$Genotype
Y = sweep(Y,1,rowMeans(Y),'-')
Y = sweep(Y,1,apply(Y,1,sd),'/')
# data = data[1:50,]
# Y = Y[1:50,1:50]

# drop genes with correlation > 0.95 (R2 > 0.90)
cor_Y = cor(Y)
hc = hclust(as.dist(1-abs(cor_Y)))
tips = cutree(hc,h=0.1)
Y = Y[,tapply(1:ncol(Y),tips,function(x) x[1])]

stopifnot(all(data$Genotype %in% rownames(cisGenotypes)))

cisGenotypes = cisGenotypes[data$Genotype,]
extra_genes = colnames(Y)[!colnames(Y) %in% colnames(cisGenotypes)]
extra_columns = matrix(0,nr = nrow(cisGenotypes),nc = length(extra_genes))
colnames(extra_columns) = extra_genes
cisGenotypes = cbind(cisGenotypes,extra_columns)[,colnames(Y)]

stopifnot(all(colnames(cisGenotypes) == colnames(Y)))
stopifnot(all(rownames(cisGenotypes) == rownames(Y)))

Bra_data = list()
Bra_data$Y            = Y
Bra_data$Z1           = model.matrix(~Genotype+0,data)
Bra_data$X            = matrix(model.matrix(~TRT,data)[,-1],nc=1)
Bra_data$A1           = diag(1,ncol(Bra_data$Z1))
Bra_data$cisGenotypes = cisGenotypes
Bra_data$n            = nrow(Bra_data$Y)
Bra_data$p            = ncol(Bra_data$Y)
Bra_data$b            = ncol(Bra_data$X)

sim_data = Bra_data
save(sim_data,file='Sim_data.RData')
