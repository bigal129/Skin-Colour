#------------------------------------------------------------------------------------------------------------------------------------------------
# Input a vector of individuals, a single locus, a prior for that locus, and the likelihoods csv, to randomly sample genotypes.
# NA data are excluded, rather than sampling from the prior, as they provide no value
# prior must be a vector of length 4, giving probabilities of A,C,G,T respectively.
#------------------------------------------------------------------------------------------------------------------------------------------------
source('functions.R')
require(DEoptimR)
#------------------------------------------------------------------------------------------------------------------------------------------------
# generate overheads
#------------------------------------------------------------------------------------------------------------------------------------------------
likelihoods <- read.csv('../data/summary.csv',row.names=1)
selected <- read.csv('../data/Selected genes.csv')
meta <- read.csv('../data/meta.csv')
meta <- meta[,c('master_ID','lat','long','mean_date_BP')]
#------------------------------------------------------------------------------------------------------------------------------------------------
# read variables
#------------------------------------------------------------------------------------------------------------------------------------------------
milk <- read.csv('../model timeseries variables/milk proportion.csv')
insol	<- read.csv('../model timeseries variables/midday insolation.csv')
clust	<- read.csv('../model timeseries variables/cluster stat.csv')
popfluc <- read.csv('../model timeseries variables/pop fluctuations stat.csv')
#------------------------------------------------------------------------------------------------------------------------------------------------
# inverse models 
#------------------------------------------------------------------------------------------------------------------------------------------------
milk.inv <- milk; milk.inv[,-c(1,2)] <- 1 - milk[,-c(1,2)]
insol.inv <- insol; insol.inv[,-c(1,2)] <- 1 - insol[,-c(1,2)]
clust.inv <- clust; clust.inv[,-c(1,2)] <- 1 - clust[,-c(1,2)] 
popfluc.inv <- popfluc; popfluc.inv[,-c(1,2)] <- 1 - popfluc[,-c(1,2)] 
#------------------------------------------------------------------------------------------------------------------------------------------------
# null model
#------------------------------------------------------------------------------------------------------------------------------------------------
null <- clust; null[,-c(1,2)] <- 1
#------------------------------------------------------------------------------------------------------------------------------------------------
# Main loop, once for each locus
#------------------------------------------------------------------------------------------------------------------------------------------------
summary <- NULL
N <- nrow(selected)
for(n in 1:N){
	locus <- selected$`RS.number`[n]

	# sample genotypes for a specific locus, given priors and likelihoods
	prior <- gtools::rdirichlet(1,as.numeric(selected[selected$`RS.number`==locus,c('A.count','C.count','G.count','T.count')])+1)
	geno <- sample.genotypes(likelihoods, individuals=row.names(likelihoods), locus, prior)
	geno <- merge(meta, geno, by.x='master_ID', by.y='individual') 

	# separate into the main four polygons of interest
	kml.path <- '../KML/polygons.model.4.kml'
	brit <- data.in.polygon(data=geno,kml.path=kml.path,index=1)
	balt <- data.in.polygon(data=geno,kml.path=kml.path,index=2)
	rhine <- data.in.polygon(data=geno,kml.path=kml.path,index=3)
	med <- data.in.polygon(data=geno,kml.path=kml.path,index=4)
	brit$region <- 'brit'
	balt$region <- 'balt'
	rhine$region <- 'rhine'
	med$region <- 'med'
	alleles.table <- rbind(brit,balt,rhine,med)

	# Parameter search and likelihoods
	variant.major <- selected$`Major.allele`[selected$`RS.number`==locus]
	variant.minor <- selected$`Minor.allele`[selected$`RS.number`==locus]


	upper <- c(1,1,1,1,1,1); lower <- c(0,0,0,0,0,0); NP <- 100; trace=T; tol <- 1e-08 # 1e-15

	null.best.major <- JDEoptim(lower=c(0,1,0,0,0,0), upper=c(1,1,1,1,1,1), fn=obj.fnc, vars=null, alleles.table=alleles.table, variant=variant.major, trace=trace, NP=NP)
	insol.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=insol, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)
	cluster.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=clust, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)
	popfluc.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=popfluc, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)

	null.best.minor <- JDEoptim(lower=c(0,1,0,0,0,0), upper=c(1,1,1,1,1,1), fn=obj.fnc, vars=null, alleles.table=alleles.table, variant=variant.minor, trace=trace, NP=NP)
	domestic.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=dom, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)
	cluster.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=clust, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)
	popfluc.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=popfluc, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)

	domestic.inv.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=dom.inv, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)
	cluster.inv.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=clust.inv, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)
	popfluc.inv.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=popfluc.inv, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)

	domestic.inv.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=dom.inv, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)
	cluster.inv.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=clust.inv, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)
	popfluc.inv.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=popfluc.inv, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)

	# Store results
	res <- data.frame(locus=locus,
			null.major.loglik=-null.best.major$value, 
			null.major.pars=paste(null.best.major$par,collapse=','),
			domestic.major.loglik=-domestic.best.major$value,
			domestic.major.pars=paste(domestic.best.major$par,collapse=','),
			cluster.major.loglik=-cluster.best.major$value,
			cluster.major.pars=paste(cluster.best.major$par,collapse=','),
			popfluc.major.loglik=-popfluc.best.major$value,
			popfluc.major.pars=paste(popfluc.best.major$par,collapse=','),

			null.minor.loglik=-null.best.minor$value, 
			null.minor.pars=paste(null.best.minor$par,collapse=','),
			domestic.minor.loglik=-domestic.best.minor$value,
			domestic.minor.pars=paste(domestic.best.minor$par,collapse=','),
			cluster.minor.loglik=-cluster.best.minor$value,
			cluster.minor.pars=paste(cluster.best.minor$par,collapse=','),
			popfluc.minor.loglik=-popfluc.best.minor$value,
			popfluc.minor.pars=paste(popfluc.best.minor$par,collapse=','),

			domestic.inv.major.loglik=-domestic.inv.best.major$value,
			domestic.inv.major.pars=paste(domestic.inv.best.major$par,collapse=','),
			cluster.inv.major.loglik=-cluster.inv.best.major$value,
			cluster.inv.major.pars=paste(cluster.inv.best.major$par,collapse=','),
			popfluc.inv.major.loglik=-popfluc.inv.best.major$value,
			popfluc.inv.major.pars=paste(popfluc.inv.best.major$par,collapse=','),

			domestic.inv.minor.loglik=-domestic.inv.best.minor$value,
			domestic.inv.minor.pars=paste(domestic.inv.best.minor$par,collapse=','),
			cluster.inv.minor.loglik=-cluster.inv.best.minor$value,
			cluster.inv.minor.pars=paste(cluster.inv.best.minor$par,collapse=','),
			popfluc.inv.minor.loglik=-popfluc.inv.best.minor$value,
			popfluc.inv.minor.pars=paste(popfluc.inv.best.minor$par,collapse=',')
			)	

	summary <- rbind(summary,res)
	print(paste(n,'of',N))
	}		
#------------------------------------------------------------------------------------------------------------------------------------------------
# save the overall parameter search results
#------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(summary,file=paste('../results/parameter.search.output.',runif(1),'.csv',sep=''),row.names=F)
#------------------------------------------------------------------------------------------------------------------------------------------------
install.packages("rgdal")