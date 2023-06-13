#--------------------------------------------------------------------
# Convert the log likelihoods to likelihood ratios
#--------------------------------------------------------------------
source('functions.R')

# choose an appropriate alpha level for the likelihood ratio test (LRT)
# there is 1 parameter difference between null (2 parameters) and alternative (6 parameters)
# so Chisquare using 4 degree freedom
# However, alpha threshold is affected by multiple testing, and needs correcting using Sidak:
# Allowing for 50 loci x 12 models is 600 multiple tests:

a <- 0.001
sidak <- ((1-a)^(1/600))
threshold <- qchisq(sidak, df=4)
#--------------------------------------------------------------------
outputs <- list.files('../results',full.names=T)
outputs <- outputs[grepl('parameter.search.output', outputs)]

for(output in outputs){
	x <- read.csv(output)
	domestic.major.LR <- round(exp(x[,4]-x[,2]),1)
	cluster.major.LR <- round(exp(x[,6]-x[,2]),1)
	popfluc.major.LR <- round(exp(x[,8]-x[,2]),1)

	domestic.minor.LR <- round(exp(x[,12]-x[,10]),1)
	cluster.minor.LR <- round(exp(x[,14]-x[,10]),1)
	popfluc.minor.LR <- round(exp(x[,16]-x[,10]),1)

	domestic.inv.major.LR <- round(exp(x[,18]-x[,2]),1)
	cluster.inv.major.LR <- round(exp(x[,20]-x[,2]),1)
	popfluc.inv.major.LR <- round(exp(x[,22]-x[,2]),1)

	domestic.inv.minor.LR <- round(exp(x[,24]-x[,10]),1)
	cluster.inv.minor.LR <- round(exp(x[,26]-x[,10]),1)
	popfluc.inv.minor.LR <- round(exp(x[,28]-x[,10]),1)

	res <- data.frame(locus=x$locus, 
		domestic.major.LR, cluster.major.LR, popfluc.major.LR, 
		domestic.minor.LR, cluster.minor.LR, popfluc.minor.LR,
		domestic.inv.major.LR, cluster.inv.major.LR, popfluc.inv.major.LR, 
		domestic.inv.minor.LR, cluster.inv.minor.LR, popfluc.inv.minor.LR)

	new.file <- gsub('parameter.search.output','likelihood ratio',output)
	write.csv(res, file=new.file,row.names=F)
	}
#--------------------------------------------------------------------
# plot likelihood ratios
#--------------------------------------------------------------------
outputs <- list.files('../results',full.names=T)
outputs <- outputs[grepl('likelihood ratio', outputs)]
max <- c()
for(n in 1:length(outputs))max <- c(max,max(read.csv(outputs[n])[,2:13]))
max <- log(max(max))
	
pdf('../plots/likelihood ratios.pdf',height=6, width=12)
par(mar=c(8,8,4,1))
for(c in 2:13){
	title <- gsub('.LR','',names(read.csv(outputs[1]))[c])
	plot(NULL,xlim=c(1,50),ylim=c(0,ceiling(max)),xlab='',xaxt='n',yaxt='n',ylab='likelihood ratio',pch=16, main=title)
	axis(1,at=1:50, labels=res[,1],las=2)
	axis(2,at=0:15, round(exp(0:15)),las=2)
	abline(v = 1:50, col='grey', lty=2)
	for(n in 1:length(outputs)){
		res <- read.csv(outputs[n])
		y <- res[,c]
		cols <- rep(1,50); cols[y>threshold] <- 2
		points(1:50,log(y),pch=16, col=scales::alpha(cols,0.5))
		}
	}
dev.off()
#--------------------------------------------------------------------
# plot allele models for anything over threshold
#--------------------------------------------------------------------
dom	<- read.csv('../model timeseries variables/domestic animal proportion.csv')
clust	<- read.csv('../model timeseries variables/cluster stat.csv')
popfluc <- read.csv('../model timeseries variables/pop fluctuations stat.csv')
dom.inv <- dom; dom.inv[,-c(1,2)] <- 1 - dom[,-c(1,2)] 
clust.inv <- clust; clust.inv[,-c(1,2)] <- 1 - clust[,-c(1,2)] 
popfluc.inv <- popfluc; popfluc.inv[,-c(1,2)] <- 1 - popfluc[,-c(1,2)] 

yearsBP <- c(dom$startBP[1],dom$endBP)

pdf('../plots/allele frequency models of interesting loci.pdf',height=3, width=12)
par(mfrow=c(1,4),mar=c(4,4,4,1))
files <- list.files('../results')
LR.files <- files[grepl('likelihood ratio', files)]
par.files <- files[grepl('parameter.search.output', files)]
N <- length(par.files)
for(c in 2:13){
	locus <- c()
	for(n in 1:N){
		LR <- read.csv(paste('../results/',LR.files[n],sep=''))
		res <- read.csv(paste('../results/',par.files[n],sep=''))
		y <- LR[,c]
		locus <- c(locus,LR$locus[y>threshold])
		model <- names(LR)[c]
		}
	locus <- unique(locus)

	title <- gsub('.LR','',model)
	short.name <- strsplit(model, split='.',fixed=T)[[1]]
	short.name <- paste(short.name[-c(length(short.name)-1,length(short.name))],collapse='.')
	if(short.name=='dom')vars <- dom
	if(short.name=='cluster')vars <- clust
	if(short.name=='popfluc')vars <- popfluc
	if(short.name=='dom.inv')vars <- dom.inv
	if(short.name=='clust.inv')vars <- clust.inv
	if(short.name=='popfluc.inv')vars <- popfluc.inv

	par <- c()
	for(loc in locus){
		for(region in 3:6){
			var <- vars[,region]
			plot(NULL, xlim=rev(range(yearsBP)), ylim=c(0,1), ylab='modelled allele frequency',xlab='yearsBP',main='')
			for(n in 1:N){
				res <- read.csv(paste('../results/',par.files[n],sep=''))
				pars <- res[res$locus==loc,paste(title,'.pars',sep='')]
				pars <- as.numeric(unlist(strsplit(pars,split=',')))
				par <- pars[c(1,2,region)]
				add.allele.frequency.curve(par, var, x=yearsBP)
				}
			title(loc,line=3)
			title(title,line=2,cex.main=0.8)
			title(names(vars)[region],line=1,cex.main=0.6)	
			}
		}
	}
dev.off()
#--------------------------------------------------------------------





