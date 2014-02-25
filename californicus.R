library(ggplot2)
library(RMySQL)
library(edgeR)
library(scales)
library(GOstats)
library(GSEABase)

	#     __________  ________________ ___      
	#    / ____/ __ \/ ____/ ____/ __ \__ \     
	#   / __/ / /_/ / /   / /   / /_/ /_/ /     
	#  / /___/ _, _/ /___/ /___ \__, / __/      
	# /_____/_/ |_|\____/\____//____/____/      
	                                          
mydb = dbConnect(MySQL(), user='californicus', password='pogo', dbname='californicus', host='ecoevo.unit.oist.jp')
ercc <- dbGetQuery(mydb,"
SELECT ercc.lib_id, ercc.gene, ercc.count, ercc.fpkm, ExFold_added.mix, ExFold_reference.`attamoles_ul` FROM (SELECT * FROM rnaseq WHERE gene REGEXP '^ERCC') AS ercc
JOIN ExFold_added
ON ercc.lib_id = ExFold_added.lib_id
JOIN ExFold_reference
ON  ExFold_added.mix = ExFold_reference.mix AND ercc.gene = ExFold_reference.id
ORDER BY ercc.lib_id, ercc.gene
")
dbDisconnect(mydb)
#figure out fpkm cutoff
ercc_min <- aggregate(ercc$fpkm[ercc$fpkm!=0],by=list(category=ercc$lib_id[ercc$fpkm!=0]),FUN=min)
cutoff_guess <- mean(ercc_min$x)

#take a guess cutoff and compute average correlation in exoected vs observed fold-change
corr_max <- function(cutoff) {
	mix1 <- levels(as.factor(subset(ercc, mix==1)$lib_id))
	mix2 <- levels(as.factor(subset(ercc, mix==2)$lib_id))
	cor_sum <-  0
	cor_count <-  0
	for (i in mix1)
		for (j in mix2)
		{
		cor_count <- cor_count + 1
		a <- subset(ercc,lib_id == i)
		b <- subset(ercc,lib_id == j)
		mdata <- na.omit(merge(a,b,by="gene"))
		keep <- apply(mdata[,c("fpkm.x","fpkm.y")],1,min) > cutoff & mdata$attamoles_ul.y > 0
		obs <- (mdata$fpkm.x/mdata$fpkm.y)[keep]
		exp <- (mdata$attamoles_ul.x/mdata$attamoles_ul.y)[keep]
		cor_sum <- cor_sum + cor(obs,exp,use="complete.obs")
		}
	return(cor_sum/cor_count)
}

#find best cutoff
cutoff <- optim(cutoff_guess,corr_max,method = "L-BFGS-B", lower = 0, upper = 1,control=c(fnscale=-10,trace=1))
corr_max(cutoff$value)

library(WGCNA)

# ERCC plots

#matrix plot fold-change, with kept points as solid circles, and eliminated points as hollow circles
op <- par(no.readonly = TRUE)
mix1 <- levels(as.factor(subset(ercc, mix==1)$lib_id))
mix2 <- levels(as.factor(subset(ercc, mix==2)$lib_id))
pdf('/Users/sasha/Dropbox/projects/californicus/plots/ercc_obs_exp.pdf',width=100, height=100, paper='special')
par(mfrow=c(length(mix1),length(mix2)),oma=c(0,0,0,0),mar=c(1,1,1,1))
for (i in mix1)
	for (j in mix2) {
		a = subset(ercc,lib_id == i)
		b = subset(ercc,lib_id == j)
		mdata <- na.omit(merge(a,b,by="gene"))
		mdata$keep <- 0
		mdata$keep[apply(mdata[,c("fpkm.x","fpkm.y")],1,min) > cutoff$value] <- 1
		mdata$obs <- mdata$fpkm.x/mdata$fpkm.y
		mdata$exp <- mdata$attamoles_ul.x/mdata$attamoles_ul.y
		mdata[sapply(mdata,is.infinite)] <- NA
		mdata <- mdata[complete.cases(mdata[,c("obs","exp")]),]
		with(subset(mdata, keep == 1),plot(obs, exp,log="xy",xlab="",ylab="",type = 'p', col = 'red', pch=16,cex=3))
		with(subset(mdata, keep == 0),points(obs, exp,log="xy",xlab="",ylab="",type = 'p', col = 'red',cex=3))
		legend("topleft", paste(i,j), bty="n",cex=2) 
	}
par(mfrow=c(1,1))
dev.off()
par(op)

#plot dymanic range and lower limit
pdf('/Users/sasha/Dropbox/projects/californicus/plots/dynamic range.pdf',width=10, height=8, paper='special')
fmt <- function(x) prettyNum(x,format="f",digits=2)
ercc$mix[ercc$mix == 1] <- "Mix 1"
ercc$mix[ercc$mix == 2] <- "Mix 2"
ggplot(data=ercc, aes(x=attamoles_ul,y=fpkm,color=lib_id))+geom_point()+scale_x_continuous(trans = log2_trans(),labels = fmt)+scale_y_continuous(trans = log2_trans())+theme_bw()+facet_grid(.~ mix)+xlab("concentration of spike-in")+theme(legend.position="none")+geom_hline(yintercept=cutoff$value,color="red")
dev.off()

	#        __                                __           _     
	#   ____/ /___ ____     ____ _____  ____ _/ /_  _______(_)____
	#  / __  / __ `/ _ \   / __ `/ __ \/ __ `/ / / / / ___/ / ___/
	# / /_/ / /_/ /  __/  / /_/ / / / / /_/ / / /_/ (__  ) (__  ) 
	# \__,_/\__, /\___/   \__,_/_/ /_/\__,_/_/\__, /____/_/____/  
	#      /____/                            /____/               

#get experimental factors
mydb = dbConnect(MySQL(), user='californicus', password='pogo', dbname='californicus', host='ecoevo.unit.oist.jp')
experimental_factors <- dbGetQuery(mydb,"SELECT * FROM factors")
rownames(experimental_factors) <- experimental_factors$id
experimental_factors <- experimental_factors[,-1]
group <- factor(paste(experimental_factors$phenotype,experimental_factors$context,sep=""))

#get transcripts and apply cutoff (this fetch takes a while)
counts <- dbGetQuery(mydb,
	"SELECT lib_id, gene,fpkm,count as coverage FROM rnaseq WHERE gene NOT REGEXP '^ERCC' " )
blast <- dbGetQuery(mydb,"SELECT * FROM blast ")
rownames(blast) <- blast$comp
go <-  dbGetQuery(mydb,"SELECT * FROM go")
go$evidence <- "ISS"
go <- go[,c("GO","evidence","gene")]
dbDisconnect(mydb)

fpkm <- unstack(counts,fpkm ~ lib_id) 
genes <- unstack(counts,coverage ~ lib_id)
rownames(genes) <- counts[counts$lib_id == "HA_153G",c("gene")]
rownames(fpkm) <- counts[counts$lib_id == "HA_153G",c("gene")]
# keep genes with at least half the libraries above the cutoff threshold
keep <- rowSums(fpkm > cutoff$value) > ncol(genes)/2

#we look only at transcripts with BLAST hits. Note blast hits were first filtered through the keep abundance filter above

keep <- rownames(genes) %in% blast$comp 

#calculate main model
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
y <- DGEList(counts=genes[keep,rownames(experimental_factors)],group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y,verbose=TRUE)
y <- estimateTagwiseDisp(y, verbose=TRUE)
fit <- glmFit(y, design)

my.contrasts <- makeContrasts(
# Haplos vs. pleos across all social contexts (singletons and pairs)
H_P = Hs + HH + HP - Ps - PP - PH,
# Haplos vs. pleos, singletons only
Hs_Ps = Hs - Ps,
# Haplos, singletons vs. in pairs
Hs_HH = Hs - HH,
# Pleos, singletons vs. in pairs
Ps_PP = Ps - PP,
# Haplos, in pure vs. mixed pairs (with other haplo vs. with pleo)
HH_HP = HH - HP,
# Pleos, in pure vs. mixed pairs (with other pleo vs. with haplo)
PP_PH = PP - PH,
# Aggressive vs. non-aggressive pairs
A_C = (HA + PA)/2 - (HH + PP + HP + PH)/4,
# Aggressive vs. singletons
A_s = HA + PA - Hs - Ps,
levels=design)

for (i in 1:ncol(my.contrasts)) {
	print(i)
	lrt <- glmLRT(fit, contrast=my.contrasts[,i])
	dt <- decideTestsDGE(lrt)
	print(colnames(my.contrasts)[i])
	print(summary(dt))
	z <- gzfile(paste0("/Users/sasha/Dropbox/projects/californicus/tests/",colnames(my.contrasts)[i],"_genes.csv.gz"),"w")
	write.csv(cbind(subset(lrt$table,select=-PValue),data.frame(significant=dt,blast=blast[rownames(lrt$table),"hit"])),z)
	close(z)
}

mds <- plotMDS(y)

ggplot(cbind(data.frame(mds1=mds$cmdscale.out[,1],mds2=mds$cmdscale.out[,2]),experimental_factors),aes(x=mds1,y=mds2,color=factor(phenotype),shape=factor(context)))+geom_point(size=5)+theme_bw()+scale_color_manual(values=c("red","blue"))+theme(legend.justification=c(1,0), legend.position=c(.3,.6))
ggsave("/Users/sasha/Dropbox/projects/californicus/plots/mds_genes.pdf")


#go term analysis

universe <- unique(go$gene)
goFrame=GOFrame(go[go$gene %in% universe,],organism="Pogonomyrmex californicus")
goAllFrame=GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# H_P
lrt <- glmLRT(fit, contrast=my.contrasts[,"H_P"])
dt <- decideTestsDGE(lrt)
H_upreg <- hyperGTest(GSEAGOHyperGParams(name = "haplo upregulated",
	geneSetCollection=gsc,geneIds = rownames(lrt$table[dt==1,]),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
summary(H_upreg)
write.table(summary(H_upreg),"/Users/sasha/Dropbox/projects/californicus/tests/go/H_upreg.txt")


P_upreg <- hyperGTest(GSEAGOHyperGParams(name = "pleo upregulated",
	geneSetCollection=gsc,geneIds = rownames(lrt$table[dt==-1,]),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
summary(P_upreg)
write.table(summary(P_upreg),"/Users/sasha/Dropbox/projects/californicus/tests/go/P_upreg.txt")


# A_s
lrt <- glmLRT(fit, contrast=my.contrasts[,"A_s"])
dt <- decideTestsDGE(lrt)
A_upreg <- hyperGTest(GSEAGOHyperGParams(name = "aggressive  upregulated vs singletons",
	geneSetCollection=gsc,geneIds = rownames(lrt$table[dt==1,]),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
summary(A_upreg)
write.table(summary(A_upreg),"/Users/sasha/Dropbox/projects/californicus/tests/go/A_upreg.txt")

s_upreg <- hyperGTest(GSEAGOHyperGParams(name = "singletons upregulated vs aggressive",
	geneSetCollection=gsc,geneIds = rownames(lrt$table[dt==-1,]),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
summary(s_upreg)
write.table(summary(s_upreg),"/Users/sasha/Dropbox/projects/californicus/tests/go/s_upreg.txt")

# A_C
lrt <- glmLRT(fit, contrast=my.contrasts[,"A_C"])
dt <- decideTestsDGE(lrt)
A_upreg_C <- hyperGTest(GSEAGOHyperGParams(name = "aggressive  upregulated vs non-aggressive",
	geneSetCollection=gsc,geneIds = rownames(lrt$table[dt==1,]),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
summary(A_upreg_C)
write.table(summary(A_upreg_C),"/Users/sasha/Dropbox/projects/californicus/tests/go/A_upreg_C.txt")

C_upreg <- hyperGTest(GSEAGOHyperGParams(name = "non-aggressive upregulated vs aggressive",
	geneSetCollection=gsc,geneIds = rownames(lrt$table[dt==-1,]),
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
summary(C_upreg)
write.table(summary(C_upreg),"/Users/sasha/Dropbox/projects/californicus/tests/go/C_upreg.txt")



#analysis of upreg modules

for (i in 1:1) {
	module <- as.character(read.table(paste0("/Users/sasha/Dropbox/projects/californicus/ACvsd_Modules/ACvsd_Module",i,".txt"),header=F)$V1)
	upreg <- hyperGTest(GSEAGOHyperGParams(name = i,
	geneSetCollection=gsc,geneIds = module,
	universeGeneIds=universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over"))
	print(summary(upreg))
}



