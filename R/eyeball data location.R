#-------------------------------------------------------------------------------------------------
# eyeball data
#-------------------------------------------------------------------------------------------------
require(maps)
require(mapdata)
require(maptools)
#-------------------------------------------------------------------------------------------------
x <- read.csv('../data/summary.csv',row.names=1)
meta <- read.csv('../data/meta.csv')
#-------------------------------------------------------------------------------------------------
# concentrate on the 15000 to 500 year range
i <- !is.na(meta$long) & !is.na(meta$lat) & meta$mean_date_BP>500 & meta$mean_date_BP<15000
meta <- meta[i,]
x <- x[i,]
#-------------------------------------------------------------------------------------------------
# separate into equal time quantiles and plot
kmls <- getKMLcoordinates('../kml/polygons.model.4.kml',TRUE)
xlim <- range(meta$long)
ylim <- range(meta$lat)

pdf(file='../plots/location of aDNA.pdf',width=20,height=8)
par(mfrow=c(2,3))
N <- 6
posts <- round(as.numeric(quantile(meta$mean_date_BP,probs=seq(0,1,length.out=N+1))))
for(n in 1:N){
	text <- paste(posts[n],'to',posts[n+1],'BP')
	sub <- subset(meta,mean_date_BP>=posts[n] & mean_date_BP<=posts[n+1])
	x <- jitter(sub$long,amount=1)
	y <- jitter(sub$lat,amount=0.5)
	plot(NULL,xlim=xlim,ylim=ylim,frame.plot=F,axes=F, xlab='',ylab='',main=text)
	map('world',xlim=xlim,ylim=ylim,col='#F1F8E9',add=T, fill=T, border='grey')
	points(x,y,pch=16,col=scales::alpha('firebrick',0.3),main=text,xlab='long',ylab='lat')
	for(k in 1:length(kmls))lines(kmls[[k]],col='blue')
	}
dev.off()
#-------------------------------------------------------------------------------------------------
