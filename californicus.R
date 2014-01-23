library(ggplot2)
library(RMySQL)
library(edgeR)
library(scales)

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

#read contrasts
contrasts <- dbGetQuery(mydb, "SELECT * FROM contrasts" )
rownames(contrasts) <- contrasts$name
contrasts <- contrasts[,-1]

#get transcripts and apply cutoff (this fetch takes a while)
counts <- dbGetQuery(mydb,
	"SELECT lib_id, gene,fpkm,count as coverage FROM rnaseq WHERE gene NOT REGEXP '^ERCC' " )
dbDisconnect(mydb)

fpkm <- unstack(counts,fpkm ~ lib_id) 
genes <- unstack(counts,coverage ~ lib_id)
rownames(genes) <- counts[counts$lib_id == "HA_153G",c("gene")]
rownames(fpkm) <- counts[counts$lib_id == "HA_153G",c("gene")]
# keep genes with at least half the libraries above the cutoff threshold
keep <- rowSums(fpkm > cutoff$value) > ncol(genes)/2

#calculate main model
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
y <- DGEList(counts=genes[keep,rownames(experimental_factors)],group=group)
y <- calcNormFactors(y)
y <- estimateCommonDisp(y,verbose=TRUE)
y <- estimateTagwiseDisp(y, verbose=TRUE)
fit <- glmFit(y, design)

my.contrasts <- makeContrasts(H_P = HA + HH + HP + Hs - PA - PH - PP - Ps, 
								Hs_Ps = Hs - Ps,
								Hs_HH = Hs - HH,
								Ps_PP = Ps - PP,
								HH_HP = HH - HP,
								PP_HP = PP - HP,
								levels=design)
lrt <- glmLRT(fit, contrast=my.contrasts[,"H_P"])
summary(dt <- decideTestsDGE(lrt))
lrt <- glmLRT(fit, contrast=my.contrasts[,"Hs_Ps"])
summary(dt <- decideTestsDGE(lrt))




#normalize contrasts, which are specified as ones and zeros (this is actually not necessary for the tests below)

#in logFC negative values are upregulated in positive contrasts
#note, this stage takes a while, so we will write the results to file at the end, so that this analysis can be resumed later
for(i in 1:nrow(contrasts)) {
	averaged_genes <- genes
	pair_id <- factor(sapply(colnames(averaged_genes),substr,4,6))  # codes correspond to experimental groups
	drop <-  c()
	# in all but the last contrast, average queens form the same treatment
	if (i < nrow(contrasts)) {
		for(j in levels(pair_id)) {
			if (sum(pair_id == j) > 1) {
				averaged_genes[,which(pair_id == j)[1]] <- rowMeans(genes[,pair_id == j])
				drop <- union(drop,c(which(pair_id == j)[-1]))
			}
		}
	}
	#drop columns that were averaged
	keep_columns <- setdiff(which(contrasts[i,]!=0),drop)
	# set up comparison vector for the test
	comparison <- c(rep(1,sum(keep_columns %in% which(contrasts[i,]>0))),rep(2,sum(keep_columns %in% which(contrasts[i,]<0))))
	# re-order columns so that they are matched up with the comparison vector
	keep_columns <- c(keep_columns[keep_columns %in% which(contrasts[i,]>0)],keep_columns[keep_columns %in% which(contrasts[i,]<0)])
	keep <- rowSums(averaged_genes[,keep_columns] > cutoff$value) > length(keep_columns)/2
	y <- DGEList(counts=averaged_genes[keep,keep_columns],group=comparison)
	y <- calcNormFactors(y)
	print("estimating common dispersion")
	y <- estimateCommonDisp(y,verbose=TRUE)
	print("estimating tagwise dispersion")
	y <- estimateTagwiseDisp(y, verbose=TRUE)
	print("significance testing")
	et <- exactTest(y)
	et$table$p_adj <- p.adjust(et$table$PValue,method="fdr")
	dge_list[[i]] <- list(talbe = et$table)
	if (nrow(et$table[et$table$p_adj<0.05,]) !=0) {
		out <- et$table[et$table$p_adj<0.05,]
		out$contrast <- rownames(contrasts[i,])
		print(out)
	}
#	write.csv(et$table,paste("/Volumes/mikheyev/Sasha/californicus/et/",rownames(contrasts[1,]),".csv",sep=""))
}

#write results to file