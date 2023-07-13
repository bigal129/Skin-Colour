#--------------------------------------------------------------------
# Convert the log likelihoods to likelihood ratios
#--------------------------------------------------------------------
source('functions.R')

# choose an appropriate alpha level for the likelihood ratio test (LRT)
# there is 1 parameter difference between null (5 parameters) and alternative (6 parameters)
# so Chisquare using 1 degree freedom
# However, alpha threshold is affected by multiple testing, and needs correcting using Sidak:
# Allowing for 15 loci x 16 models is 240 multiple tests:

a <- c(0.05,0.001) # 0.001 is highly conservative. 0.05 is more 'common'
sidak <- ((1-a)^(1/240))
threshold <- qchisq(sidak, df=1)
#--------------------------------------------------------------------
outputs <- list.files('../results',full.names=T)
outputs <- outputs[grepl('parameter.search.output', outputs)]

for(output in outputs){
	x <- read.csv(output)

	cluster.major.LR <- round(exp(x$cluster.major.loglik-x$null.major.loglik),1)
	popfluc.major.LR <- round(exp(x$popfluc.major.loglik-x$null.major.loglik),1)
	insol.major.LR <- round(exp(x$insol.major.loglik-x$null.major.loglik),1)
	milk.major.LR <- round(exp(x$milk.major.loglik-x$null.major.loglik),1)

	cluster.minor.LR <- round(exp(x$cluster.minor.loglik-x$null.minor.loglik),1)
	popfluc.minor.LR <- round(exp(x$popfluc.minor.loglik-x$null.minor.loglik),1)
	insol.minor.LR <- round(exp(x$insol.minor.loglik-x$null.minor.loglik),1)
	milk.minor.LR <- round(exp(x$milk.minor.loglik-x$null.minor.loglik),1)

	cluster.inv.major.LR <- round(exp(x$cluster.inv.major.loglik-x$null.major.loglik),1)
	popfluc.inv.major.LR <- round(exp(x$popfluc.inv.major.loglik-x$null.major.loglik),1)
	insol.inv.major.LR <- round(exp(x$insol.inv.major.loglik-x$null.major.loglik),1)
	milk.inv.major.LR <- round(exp(x$milk.inv.major.loglik-x$null.major.loglik),1)

	cluster.inv.minor.LR <- round(exp(x$cluster.inv.minor.loglik-x$null.minor.loglik),1)
	popfluc.inv.minor.LR <- round(exp(x$popfluc.inv.minor.loglik-x$null.minor.loglik),1)
	insol.inv.minor.LR <- round(exp(x$insol.inv.minor.loglik-x$null.minor.loglik),1)
	milk.inv.minor.LR <- round(exp(x$milk.inv.minor.loglik-x$null.minor.loglik),1)

	res <- data.frame(locus=x$locus, 
		cluster.major.LR, popfluc.major.LR, insol.major.LR, milk.major.LR,
		cluster.minor.LR, popfluc.minor.LR, insol.minor.LR, milk.minor.LR,
		cluster.inv.major.LR, popfluc.inv.major.LR, insol.inv.major.LR, milk.inv.major.LR,
		cluster.inv.minor.LR, popfluc.inv.minor.LR, insol.inv.minor.LR, milk.inv.minor.LR)

	new.file <- gsub('parameter.search.output','likelihood ratio',output)
	write.csv(res, file=new.file,row.names=F)
	}
#--------------------------------------------------------------------
# plot likelihood ratios
#--------------------------------------------------------------------
outputs <- list.files('../results',full.names=T)
outputs <- outputs[grepl('likelihood ratio', outputs)]
max <- c()
C <- ncol(read.csv(outputs[1]))
for(n in 1:length(outputs))max <- c(max,max(read.csv(outputs[n])[,2:C]))
max <- log(max(max))
	
pdf('../plots/likelihood ratios.pdf',height=6, width=12)
par(mar=c(8,8,4,1))
for(c in 2:C){
	R <- nrow(res)
	title <- gsub('.LR','',names(read.csv(outputs[1]))[c])
	plot(NULL,xlim=c(1,R),ylim=c(0,ceiling(max)),xlab='',xaxt='n',yaxt='n',ylab='likelihood ratio',pch=16, main=title)
	axis(1,at=1:R, labels=res[,1],las=2)
	axis(2,at=0:ceiling(max), round(exp(0:ceiling(max))),las=2)
	abline(v = 1:R, col='grey', lty=2)
	for(n in 1:length(outputs)){
		res <- read.csv(outputs[n])
		y <- res[,c]
		cols <- rep(1,R)
		points(jitter(1:R),log(y),pch=16, col=scales::alpha(cols,0.15),cex=1.6)
		}
	abline(h=log(threshold),col=2,lty=2)
	}
dev.off()
#--------------------------------------------------------------------
# plot allele models for anything really interesting (mean over threshold)
#--------------------------------------------------------------------
milk <- read.csv('../model timeseries variables/milk proportion.csv')
insol	<- read.csv('../model timeseries variables/midday insolation.csv')
clust	<- read.csv('../model timeseries variables/cluster stat.csv')
popfluc <- read.csv('../model timeseries variables/pop fluctuations stat.csv')

# inverse models 
milk.inv <- milk; milk.inv[,-c(1,2)] <- 1 - milk[,-c(1,2)]
insol.inv <- insol; insol.inv[,-c(1,2)] <- 1 - insol[,-c(1,2)]
clust.inv <- clust; clust.inv[,-c(1,2)] <- 1 - clust[,-c(1,2)] 
popfluc.inv <- popfluc; popfluc.inv[,-c(1,2)] <- 1 - popfluc[,-c(1,2)] 

yearsBP <- c(milk$startBP[1],milk$endBP)
regions <- names(milk)
regions <- gsub('.MAP','',regions)

files <- list.files('../results')
LR.files <- files[grepl('likelihood ratio', files)]
par.files <- files[grepl('parameter.search.output', files)]
N <- length(par.files)
C <- ncol(read.csv(outputs[1]))

# get mean LRs
loci <- read.csv(paste('../results/',LR.files[1],sep=''))[,1]
tmp <- read.csv(paste('../results/',LR.files[1],sep=''))[,-1]
for(n in 2:N) tmp <- tmp + read.csv(paste('../results/',LR.files[n],sep=''))[,-1]
tmp <- tmp/N

pdf('../plots/allele frequency models of interesting loci.pdf',height=3, width=12)
for(c in 1:(C-1)){
	for(l in 1:length(loci)){
		if(tmp[l,c] > max(threshold)){
			model <- names(tmp)[c]
			locus <- loci[l]
			short.name <- strsplit(model, split='.',fixed=T)[[1]]
			short.name <- paste(short.name[-c(length(short.name)-1,length(short.name))],collapse='.')
			if(short.name=='milk')vars <- milk
			if(short.name=='insol')vars <- insol
			if(short.name=='cluster')vars <- clust
			if(short.name=='popfluc')vars <- popfluc
			if(short.name=='milk.inv')vars <- milk.inv
			if(short.name=='insol.inv')vars <- insol.inv
			if(short.name=='cluster.inv')vars <- clust.inv
			if(short.name=='popfluc.inv')vars <- popfluc.inv

			par(mfrow=c(1,4),mar=c(4,4,4,1))
			for(region in 3:6){	
				plot(NULL, xlim=rev(range(yearsBP)), ylim=c(0,1), ylab='modelled allele frequency',xlab='yearsBP',main=locus)
				text(x=mean(yearsBP),y=0.88,labels=regions[region])
				text(x=mean(yearsBP),y=0.95,labels=gsub('.LR','',model), cex=1.3)	
				var <- vars[,region]
				for(n in 1:N){
					res <- read.csv(paste('../results/',par.files[n],sep=''))
					pars <- res[res$locus==locus,gsub('LR','pars',model)]
					pars <- as.numeric(unlist(strsplit(pars,split=',')))
					par <- pars[c(1,2,region)]
					add.allele.frequency.curve(par, var, x=yearsBP, alpha=0.05)
					}
				}
			}
		}
	}
dev.off()
#--------------------------------------------------------------------





