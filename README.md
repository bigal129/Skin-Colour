# Skin-Colour
Masters thesis - skin colour selection

Finished the main analysis coding with 4 different proxies: clustering, population fluctuation, milk usage and insolation.
Using major, minor and inverse of both 
This analysis is done over 12 genes in 15 loci
Compared agaisnt a major and minor null 

Outputted a results csv for this 
Outputted 100x for sensitivity testing

## Overview


15 loci of interest (some association with skin colour) are selected, and tested for signals of selection across the period 8000 BP to 2500 BP.
Specfically, we test if alleles at each locus (independently analysed) might be better explained by a model of constant selection (the null model),
or various alternative models where the strength of selection is somehow influenced (modulated) by some other hypothesised variable.
The null model therefore has two basic free parameters: the intial allele frequency (at 8000 BP) and the constant selection strength.
Meanwhile the alternative model as an extra parameter, that determines the extent to which the hypthesised variable modulates the selection strength.
In practice, this is complicated slightly by the fact that we are looking at selection in four different regions (polygons): 
British Isles; Baltic region; Rhine Danube axis; Mediterranean Europe.
The selection strength parameter (and the modulation parameter for the alternative models) are fixed across these four regions, whilst the initial allele frequency is free to vary in each region.
Therefore, the null model has five parameters (constant selection strength, intial allele frequency in each region), whilst the alternative models have six parameters (the additional modulation parameter).

Each variable that is hypothesised to influence the selection strength (termed the 'environmental variables') are time-series that vary between 0 and 1, across the period 8000 BP to 2500 BP.
For example, the enviromental variable 'residential clustering' comprises four time series (one for each region), that vary between 0 and 1, indicating periods where prehistoric settlements were more densely clustered , or more evenly spread out. 
Therefore, the parameters are simultaeously being optimised to find the maximum probability of the observed data (allele calls from ancient human DNA found in each respective polygon), across both space and time.
The only exception to this is the environmental variable 'midday insolation', which of course varies regionally on the surface of the earth (and therefore has different values in each polygon), but is a fixed value through time.

All models (both null and alternative) are constrained to explore positive selection. Therefore all make the key assumption that allele frequencies can only increase (in each polygon, during the study period). 
This of course includes the possibiloty of a zero increase, but crucially the model can never allow the allele frequency to decrease.

All 15 loci of interest happen to be bi-allelic (the analysis can be trivially extended to test tri-allelic loci if required), but make no assumptions about which of the two alleles might be under selection, nor which allele is ancestral.
Instead, we test both alleles (independently) at each locus (independently) for each of the four environmental variables. 
Furthermore, we make no assuption about in which direction the environmetal variable might effect selection - for example if an increase in settlement density causes in increase in selection, or if a decrease in settlement density causes an increase in selection.
Therefore, we also test four respective 'inverse models' ( 1 minus the values of each model).
This results in a total of (4 models + 4 inverse models) x 2 alleles x 15 loci, resulting in 240 independent tests.
Therefore our threshold of significance requires adapting our alpha value using the Sidak correction, and our Likelihood Ratio Test (LRT) comparing the null model to each alternative model also requires adjusting for one additional parameter (5 vs 6 respectively).
We choose two alpha values: 0.05, and a highly conservative value of 0.001.
Therefore our LRT thresholds are 13.7 and 21.2 calculated in R as follows:

a <- c(0.05,0.001) 
sidak <- ((1-a)^(1/240))
threshold <- qchisq(sidak, df=1)

## Repeated tests

A single LRT by comparing the model likelihood of selection modulated by a single environmental variable, compared to the null model likelihood of constant selection, requires calculating the probability of the observed data (the aDNA alleles) under the derived allele frequency curves.
However, the aDNA observed alleles ahve considerable uncertainty, largely due to the number of low reads (typical of aDNA). 
For example, the most likely genotype for an individual witha single read of 'A' is of course AA, but there is nevertheless a high chance the true genotype might be AC, or AG etc.
Therefore we take a Bayesian approach generating a probability of each of the 10 possible genotypes at each locus.

The first prior is that each of the four nucleotide (A,C,G,T) are equally probable. 
This is updated by the modern observed allele frequencies in the GNomAD database (with no exclusions). These four nucleotide probabilties are (trivially) coverted to 10 genotype probabilities.
This is updated by the PHRED scores provided by the aDNA pileup, that itself incorporates both read uncertainty and read depth.
Where there is a large read depth, these priors are washed out, and the data dominates.
These posterior probabilities for each genotype are used to randomly sample a single genotype for each ancient individual at each locus.
Therefore, this final sampling is performed many (100 + ) times, to build a distribution of ikelihoods, and subsequent LRTs.