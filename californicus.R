library(ggplot2)
library(RMySQL)
library(edgeR)
library(scales)

mydb = dbConnect(MySQL(), user='californicus', password='pogo', dbname='californicus', host='ecoevo.unit.oist.jp')

	#     __________  ________________ ___      
	#    / ____/ __ \/ ____/ ____/ __ \__ \     
	#   / __/ / /_/ / /   / /   / /_/ /_/ /     
	#  / /___/ _, _/ /___/ /___ \__, / __/      
	# /_____/_/ |_|\____/\____//____/____/      
	                                          
ercc <- dbGetQuery(mydb,"
SELECT ercc.lib_id, ercc.transcript_id, ercc.descriptive_name,ercc.coverage, ercc.FPKM, ExFold_added.mix, ExFold_reference.`attamoles_ul` FROM (SELECT * FROM rnaseq WHERE descriptive_name REGEXP '^ERCC') AS ercc
JOIN ExFold_added
ON ercc.lib_id = ExFold_added.lib_id
JOIN ExFold_reference
ON  ExFold_added.mix = ExFold_reference.mix AND ercc.descriptive_name = ExFold_reference.id
ORDER BY ercc.lib_id, ercc.transcript_id
")

#figure out FPKM cutoff
ercc_min <- aggregate(ercc$FPKM[ercc$FPKM!=0],by=list(category=ercc$lib_id[ercc$FPKM!=0]),FUN=min)
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
		mdata <- na.omit(merge(a,b,by="transcript_id"))
		keep <- apply(mdata[,c("FPKM.x","FPKM.y")],1,min) > cutoff & mdata$attamoles_ul.y > 0
		obs <- (mdata$FPKM.x/mdata$FPKM.y)[keep]
		exp <- (mdata$attamoles_ul.x/mdata$attamoles_ul.y)[keep]
		cor_sum <- cor_sum + cor(obs,exp,use="complete.obs")
		}
	return(cor_sum/cor_count)
}

#find best cutoff
cutoff <- optim(cutoff_guess,corr_max,method = "L-BFGS-B", lower = 0, upper = 1)
corr_max(cutoff$value)

#matrix plot fold-change, with kept points as solid circles, and eliminated points as hollow circles
op <- par(no.readonly = TRUE)
mix1 <- levels(as.factor(subset(ercc, mix==1)$lib_id))
mix2 <- levels(as.factor(subset(ercc, mix==2)$lib_id))
pdf('/Volumes/mikheyev/Sasha/californicus/plots/ercc_obs_exp.pdf',width=100, height=100, paper='special')
par(mfrow=c(length(mix1),length(mix2)),oma=c(0,0,0,0),mar=c(1,1,1,1))
for (i in mix1)
	for (j in mix2) {
		a = subset(ercc,lib_id == i)
		b = subset(ercc,lib_id == j)
		mdata <- na.omit(merge(a,b,by="gene_id"))
		mdata$keep <- 0
		mdata$keep[apply(mdata[,c("FPKM.x","FPKM.y")],1,min) > cutoff$value] <- 1
		mdata$obs <- mdata$FPKM.x/mdata$FPKM.y
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
pdf('/Volumes/mikheyev/Sasha/californicus/plots/dynamic range.pdf',width=10, height=8, paper='special')
fmt <- function(x) prettyNum(x,format="f",digits=2)
ercc$mix[ercc$mix == 1] <- "Mix 1"
ercc$mix[ercc$mix == 2] <- "Mix 2"
ggplot(data=ercc, aes(x=attamoles_ul,y=FPKM,color=lib_id))+geom_point()+scale_x_continuous(trans = log2_trans(),labels = fmt)+scale_y_continuous(trans = log2_trans())+theme_bw()+facet_grid(.~ mix)+xlab("concentration of spike-in")+theme(legend.position="none")+geom_hline(yintercept=cutoff$value,color="red")
dev.off()

	#        __                                __           _     
	#   ____/ /___ ____     ____ _____  ____ _/ /_  _______(_)____
	#  / __  / __ `/ _ \   / __ `/ __ \/ __ `/ / / / / ___/ / ___/
	# / /_/ / /_/ /  __/  / /_/ / / / / /_/ / / /_/ (__  ) (__  ) 
	# \__,_/\__, /\___/   \__,_/_/ /_/\__,_/_/\__, /____/_/____/  
	#      /____/                            /____/               

#get regular transcripts, using the cutorr
counts <- dbGetQuery(mydb,
	"SELECT lib_id, transcript_id,locus_id,descriptive_name,coverage FROM rnaseq WHERE descriptive_name NOT REGEXP '^ERCC'" )

genes <- data.frame(gene_id = counts[counts$lib_id == "HA_153G",c("gene_id")])
rownames(genes) <- counts[counts$lib_id == "HA_153G",c("tracking_id")]
genes <- rbind(genes,unstack(counts,coverage ~ lib_id))
# keep genes with at least half the libraries above the threshols
keep <- rowSums(genes > cutoff$value) > ncol(genes)/2

contrasts <- dbGetQuery(mydb, "SELECT * FROM contrasts" )
rownames(contrasts) <- contrasts$name
contrasts <- contrasts[,-1]
#normalize contrasts, which are specified as ones and zeros (this is actually not necessary for the tests below)
for(i in 1:nrow(contrasts)) {
	contrasts[i,contrasts[i,] == 1 ] <- contrasts[i,contrasts[i,] == 1 ]/sum(contrasts[i,] == 1 )
	contrasts[i,contrasts[i,] == -1 ] <- contrasts[i,contrasts[i,] == -1 ]/abs(sum(contrasts[i,] == -1 ))
}

dge_list <- list()
pdf('/Volumes/mikheyev/Sasha/californicus/plots/nmds.pdf',width=10, height=8, paper='special')
for(i in 1:nrow(contrasts)) {
	pair <- factor(sapply(colnames(genes),substr,4,6))
	group <- factor(c(rep(1,sum(contrasts[i,]>0)),rep(2,sum(contrasts[i,]<0))))
	keep_columns <- c(colnames(contrasts)[contrasts[i,]>0],colnames(contrasts)[contrasts[i,]<0])
	keep <- rowSums(genes[,keep_columns] > cutoff$value) > length(keep_columns)/2
	y <- DGEList(counts=genes[keep,keep_columns],group=group)
	y <- calcNormFactors(y)
	y <- estimateCommonDisp(y)
	y <- estimateTagwiseDisp(y)
	et <- exactTest(y)
	de <- decideTestsDGE(et,p=0.05, adjust="BH")
	print(paste("contrast"),rownames(contrasts[i,]))
	print(table(de))
	dge_list[[rownames(contrasts[i,])]] <- list(et = et, de = de)
	mdsplot <- plotMDS(y,top=100)
	plot(mdsplot, main=rownames(contrasts[i,]),xlab="Dimension 1",ylab="Dimension 2")
	points(mdsplot$cmdscale.out[colnames(contrasts)[contrasts[i,]>0],1],mdsplot$cmdscale.out[colnames(contrasts)[contrasts[i,]>0],2],pch=16,col="red",cex=2)
	points(mdsplot$cmdscale.out[colnames(contrasts)[contrasts[i,]<0],1],mdsplot$cmdscale.out[colnames(contrasts)[contrasts[i,]<0],2],pch=16,col="blue",cex=2)
	# text(mdsplot$cmdscale.out[,1],mdsplot$cmdscale.out[,2],keep_columns)
}
dev.off()




# plot NMDS as a sanity check

# mdsplot <- plotMDS(dge,top=500)
# plot(mdsplot, main="",xlab="Dimension 1",ylab="Dimension 2")
# text(mdsplot$cmdscale.out[,1],mdsplot$cmdscale.out[,2],factors$factor)

#testing contrasts
dge_list <- list() # this list stores results of likelihood ratio tests, and a list of differentially expressed genes, after multiple comparison correction
for (i in 1:nrow(contrasts)) {
	lrt <- glmLRT(fit,contrast=t(contrasts[i,]))
	dge_list[[i]] <- list(lrt = lrt, lrt_et = decideTestsDGE(lrt,p=0.05,adjust="BH"))
}
