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
# read variables #1
cluster <- read.csv('../model timeseries variables/cluster stats.csv')
insol <- read.csv('../model timeseries variables/midday insolation.csv')
popfluc <- read.csv('../model timeseries variables/population fluctuation.csv')
milk <- read.csv('../model timeseries variables/milk proportion.csv')
#------------------------------------------------------------------------------------------------------------------------------------------------
# read variables #2
#------------------------------------------------------------------------------------------------------------------------------------------------
stress <- read.csv('../model timeseries variables/porotic hyperostosis proportion.csv')
trauma <- read.csv('../model timeseries variables/trauma deaths.csv')
#------------------------------------------------------------------------------------------------------------------------------------------------
# inverse models 
#------------------------------------------------------------------------------------------------------------------------------------------------
stress.inv <- stress; stress.inv[,-c(1,2)] <- 1 - stress[,-c(1,2)] 
trauma.inv <- trauma; trauma.inv[,-c(1,2)] <- 1 - trauma[,-c(1,2)] 
cluster.inv <- cluster; cluster.inv[,-c(1,2)] <- 1 - cluster[,-c(1,2)]
insol.inv <- insol; insol.inv[,-c(1,2)] <- 1 - insol[,-c(1,2)]
popfluc.inv <- popfluc; popfluc.inv[,-c(1,2)] <- 1 - popfluc[,-c(1,2)]
milk.inv <- milk; milk.inv[,-c(1,2)] <- 1 - milk[,-c(1,2)]
#------------------------------------------------------------------------------------------------------------------------------------------------
# null model
#------------------------------------------------------------------------------------------------------------------------------------------------
null <- stress; null[,-c(1,2)] <- 1
#------------------------------------------------------------------------------------------------------------------------------------------------
# Main loop, once for each locus
#------------------------------------------------------------------------------------------------------------------------------------------------
for (i in 1:10){
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
	  stress.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=stress, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)
	  trauma.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=trauma, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)
 
	  null.best.minor <- JDEoptim(lower=c(0,1,0,0,0,0), upper=c(1,1,1,1,1,1), fn=obj.fnc, vars=null, alleles.table=alleles.table, variant=variant.minor, trace=trace, NP=NP)
	  stress.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=stress, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)
	  trauma.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=trauma, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)

	  stress.inv.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=stress.inv, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)
	  trauma.inv.best.major <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=trauma.inv, alleles.table=alleles.table, variant=variant.major, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.major$par)

	  stress.inv.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=stress.inv, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)
	  trauma.inv.best.minor <- JDEoptim(lower=lower, upper=upper, fn=obj.fnc, vars=trauma.inv, alleles.table=alleles.table, variant=variant.minor, trace=trace, tol=tol, NP=NP, add_to_init_pop=null.best.minor$par)

	# Store results
	  res <- data.frame(locus=locus, 
			null.major.loglik=-null.best.major$value,
			null.major.pars=paste(null.best.major$par,collapse=','),
			stress.major.loglik=-stress.best.major$value,
			stress.major.pars=paste(stress.best.major$par,collapse=','),
			trauma.major.loglik=-trauma.best.major$value,
			trauma.major.pars=paste(trauma.best.major$par,collapse=','),
			
			null.minor.loglik=-null.best.minor$value, 
			null.minor.pars=paste(null.best.minor$par,collapse=','),
		  stress.minor.loglik=-stress.best.minor$value,
			stress.minor.pars=paste(stress.best.minor$par,collapse=','),
			trauma.minor.loglik=-trauma.best.minor$value,
			trauma.minor.pars=paste(trauma.best.minor$par,collapse=','),

			stress.inv.major.loglik=-stress.inv.best.major$value,
			stress.inv.major.pars=paste(stress.inv.best.major$par,collapse=','),
			trauma.inv.major.loglik=-trauma.inv.best.major$value,
			trauma.inv.major.pars=paste(trauma.inv.best.major$par,collapse=','),

			stress.inv.minor.loglik=-stress.inv.best.minor$value,
			stress.inv.minor.pars=paste(stress.inv.best.minor$par,collapse=','),
			trauma.inv.minor.loglik=-trauma.inv.best.minor$value,
			trauma.inv.minor.pars=paste(trauma.inv.best.minor$par,collapse=',')
			)	
	summary <- rbind(summary,res)
	print(paste(n,'of',N))
  }		
print(i)
#------------------------------------------------------------------------------------------------------------------------------------------------
# save the overall parameter search results
#------------------------------------------------------------------------------------------------------------------------------------------------
write.csv(summary,file=paste('../results/parameter.search.output.',runif(1),'.csv',sep=''),row.names=F)
}
#------------------------------------------------------------------------------------------------------------------------------------------------

