#--------------------------------------------------------------------------------------------------------------------------------
# Various functions
#--------------------------------------------------------------------------------------------------------------------------------
calculate.likelihoods.from.mpileups.single.locus.single.ID <- function(data, locus, ID){
 
	x <- data[data$master_ID==ID,grepl(locus,names(data))]
	names(x) <- gsub(paste(locus,'_',sep=''),'',names(x))
	alleles <- rev(strsplit(gsub("\\[|\\]|'| ","",x$alleles),split=',')[[1]])
	if(length(alleles)==0){
		summary <- data.frame(ID=ID, locus=locus, AA.AC.AG.AT.CC.CG.CT.GG.GT.TT.likelihoods=NA)
		return(summary)
		}
	liks <- rev(10^(-as.numeric(strsplit(gsub("\\[|\\]|'| ","",x$gls),split=',')[[1]])/10))
	N <- length(alleles)
	index <- gtools::combinations(N, 2, 1:N, repeats.allowed=T)
	geno <- as.data.frame(matrix(alleles[index],(N+1)*(N/2),2)); names(geno) <- c('chr1','chr2')
	geno <- cbind(geno, liks)

	# we interpret "<*>" as an unspecified alternate allele
	nucleotides <- c('A','C','T','G')
	specified <- unique(c(geno[,1],geno[,2]))
	specified <- specified[specified%in%nucleotides]
	unspecified <- nucleotides[!nucleotides%in%specified]
	geno[geno=='<*>'] <- paste(unspecified,collapse='|')

	genotypes <- as.data.frame(gtools::combinations(4, 2, nucleotides, repeats.allowed=T))
	for(n in 1:nrow(geno)){
		cond.1 <- grepl(geno[n,1],genotypes[,1]) & grepl(geno[n,2],genotypes[,2]) 
		cond.2 <- grepl(geno[n,1],genotypes[,2]) & grepl(geno[n,2],genotypes[,1]) 
		cond <- cond.1 | cond.2
		genotypes$liks[cond] <- geno$liks[n]
		}
	genotypes$liks[genotypes$liks=='1,1,1,1,1,1,1,1,1,1'] <- NA

	summary <- data.frame(ID=ID, locus=locus, AA.AC.AG.AT.CC.CG.CT.GG.GT.TT.likelihoods=paste(genotypes$liks,collapse=','))
return(summary)}
#--------------------------------------------------------------------------------------------------------------------------------
convert.mpileups <- function(data){
	loci <- gsub('_alleles','',names(data)[grepl('alleles',names(data))])
	meta <- data[,!grepl('alleles|gls|genotype|depth',names(data))]
	ids <- data$master_ID
	summary <- as.data.frame(matrix(,length(ids),length(loci))); row.names(summary) <- ids; names(summary) <- loci

	for(r in 1:nrow(summary)){
		print(paste(r,'completed of',nrow(summary)))
		for(c in 1:ncol(summary)){
			locus <- loci[c]
			ID <- ids[r]
			x <- calculate.likelihoods.from.mpileups.single.locus.single.ID(data, locus, ID )
			summary[r,c] <- x$AA.AC.AG.AT.CC.CG.CT.GG.GT.TT.likelihoods
			}
		}
	return(list(summary=summary,meta=meta))}
#--------------------------------------------------------------------------------------------------------------------------------
sample.genotypes <- function(likelihoods, individuals, locus, prior){

	# some basic housekeeping
	if(length(locus)>1)stop('Input only one locus')
	if(sum(prior)>1.001 | sum(prior)<0.999) stop('prior should sum to 1')
	if(length(prior)!=4)stop('prior should have length of 4, corresponding to A,C,G,T respectively')
	if(!locus%in%names(likelihoods))stop('locus is missing from likelihoods')
	missing.individuals <- individuals[!individuals%in%row.names(likelihoods)]
	if(length(missing.individuals)!=0)stop(paste('The following individuals are missing:',missing.individuals))

	# convert prior of alleles to prior of genotypes: AA AC AG AT CC CG CT GG GT TT
	tmp <- gtools::combinations(4, 2, 1:4, repeats.allowed=TRUE)
	genotype.priors <- prior[tmp[,1]]*prior[tmp[,2]] * (1+(tmp[,1]!=tmp[,2]))

	# extract genotype likelihodds
	x <- likelihoods[individuals,locus,drop=F]
	x <- x[!is.na(x),,drop=F]
	keep.individuals <- row.names(x)
	genotype.likelihoods <- matrix(as.numeric(unlist(strsplit(as.matrix(x),split=','))),10,nrow(x))

	# combine likelihoods and priors to get posterior probabilities
	genotype.tmp <- genotype.likelihoods *  genotype.priors
	genotype.probs <- t(genotype.tmp)/colSums(genotype.tmp)

	# sample alleles from genotype probabilities
	N <- nrow(genotype.probs)
	genotypes <- character(N)
	for(n in 1:N)genotypes[n] <- sample(c('AA','AC','AG','AT','CC','CG','CT','GG','GT','TT'),size=1,prob=genotype.probs[n,])


	# summary
	summary <- data.frame(individual=keep.individuals, genotype=genotypes)

return(summary)}
#--------------------------------------------------------------------------------------------------------------------------------
data.in.polygon <- function(data,kml.path,polygons=NULL,index=NULL){
	# returns the data in a particular polygon, use the full kml file path of the polygon
	# alternatively, if kml.path is NULL, polygons can be used: a list of matrixes of two columns (no names): long,lat
	# index: can be used to select particular polygons only in the kml
	require(maptools)
	require(splancs)
	require(rgdal)
	data.region <- NULL
	if(!is.null(kml.path)&!is.null(polygons))stop('Which one do you want? the kml.path or the polygons? Make one of them NULL')
	if(is.null(polygons))polygons <- getKMLcoordinates(kml.path,ignoreAltitude=T)
	if(is.null(index))index <- 1:length(polygons)
	for(n in index){
		polygon <- polygons[[n]]
		if('x'%in%names(data)){
			data <- subset(data,!is.na(x))
			georefs <- pip(pts=data.frame(x=data$x,y=data$y),poly=polygon,bound=T)
			index <- data$x%in%georefs$x & data$y%in%georefs$y
			}
		if('long'%in%names(data)){
			data <- subset(data,!is.na(long))
			georefs <- pip(pts=data.frame(x=data$long,y=data$lat),poly=polygon,bound=T)
			index <- data$long%in%georefs$x & data$lat%in%georefs$y
			}
		if('Longitude'%in%names(data)){
			data <- subset(data,!is.na(Longitude))
			georefs <- pip(pts=data.frame(x=data$Longitude,y=data$Latitude),poly=polygon,bound=T)
			index <- data$Longitude%in%georefs$x & data$Latitude%in%georefs$y
			}

		data.region <- rbind(data.region,data[index,])	
		}	
return(data.region)}
#--------------------------------------------------------------------------------------------------------------------------------
obj.fnc <- function(pars, vars, alleles.table, variant){
		
		times <- vars[,c(1,2)]
		scale <- pars[1]
		weight <- pars[2]
		f.start.brit <- pars[3]
		f.start.balt <- pars[4]
		f.start.rhine <- pars[5]
		f.start.med <- pars[6]

		loglik.brit <- get.overall.loglik.of.single.variable(vars$British.Isles.MAP,times,scale,f.start.brit,weight,subset(alleles.table,region=='brit'),variant)
		loglik.balt <- get.overall.loglik.of.single.variable(vars$Baltic.region.MAP,times,scale,f.start.balt,weight,subset(alleles.table,region=='balt'),variant)
		loglik.rhine <- get.overall.loglik.of.single.variable(vars$Rhine.Danube.axis.MAP,times,scale,f.start.rhine,weight,subset(alleles.table,region=='rhine'),variant)
		loglik.med <- get.overall.loglik.of.single.variable(vars$Mediterranean.Europe.MAP,times,scale,f.start.med,weight,subset(alleles.table,region=='med'),variant)

		loglik <- loglik.brit + loglik.balt + loglik.rhine  + loglik.med
return(-loglik)}
#--------------------------------------------------------------------------------------------------------------------------------
get.overall.loglik.of.single.variable <- function(variable,times,scale,f.start,weight,alleles.table,variant){
	
	# variable: 3 column dataframe comprising 'startBP' and 'endBP' and a variable(any name)
	# the variable must be between 0 and 1
	# scale: a coefficient to scale the variable to S
	# f.start: LP frequency at the start
	# weight: coefficient between 0 and 1 to give a skewing
	# alleles.table: regions specific subset of overall.alleles.table

	exponent <- (1/weight)-1
	# floating point bullshit
	if(exponent<0.00000001)exponent <- 0
	if(length(variable)!=nrow(times))stop('variables have to match the times')

	v. <- variable^exponent

	# generate a freq.table for the variable
	freq <- create.all.freqs(v.*scale, f.start)
	time <- c(times$startBP[1], times$endBP)
	freq.table <- data.frame(time=time,freq=freq)

	# calculate the overall loglik
	loglik <- calculate.loglik(freq.table,alleles.table,variant)

return(loglik)}
#--------------------------------------------------------------------------------------------------------------------------------
calculate.loglik <- function(freq.table,alleles.table,variant){
	if(!variant%in%c('A','C','G','T'))stop('variant should be A,C,G,T only')
	prob <- approx(freq.table, xout=alleles.table[,'mean_date_BP'])$y
	mat <- matrix(strsplit(paste(alleles.table$genotype,collapse=''),split='|')[[1]],2,nrow(alleles.table))
	present <- colSums(mat==variant)
	absent <- colSums(mat!=variant)
	loglik <- log(prob^present * (1-prob)^absent)	
return(sum(loglik, na.rm=T))}
#--------------------------------------------------------------------------------------------------------------------------------
create.all.freqs <- function(S.vec, f.start){
	# S.vec: vector of selection coefficients at each time point (not necessarily generations)
	# f.start: frequency at the start
	x <- -log(1/((1/f.start)-1)) + cumsum(log(1-c(0,S.vec)))
	freqs <- 1/(exp(x)+1)
return(freqs)}
#--------------------------------------------------------------------------------------------------------------------------------
add.allele.frequency.curve <- function(pars, var, x, alpha=1){

	require(scales)
	scale <- pars[1]
	weight <- pars[2]
	f.start <- pars[3]

	exponent <- (1/weight)-1
	# floating point bullshit
	if(exponent<0.00000001)exponent <- 0
	freq <- create.all.freqs(scale * var^exponent, f.start)
	lines(x,freq, col=alpha(1,alpha))
	}
#--------------------------------------------------------------------------------------------------------------------------------
