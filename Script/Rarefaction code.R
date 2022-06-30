library(iNEXT) # This is the package developed by Chao et al. for rarefaction. 
# See : https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
# Reference papers :
# Why they use completeness rather than abundance for rarefaction :
# https://www.researchgate.net/publication/235713090_Coverage-based_rarefaction_and_extrapolation_Standardizing_samples_by_completeness_rather_than_size
# Rarefaction of all Hill numbers :
# https://www.researchgate.net/publication/273219459_Rarefaction_and_extrapolation_with_Hill_numbers_A_framework_for_sampling_and_estimation_in_species_diversity_studies


# The code :
# "communaute" is a dataframe filled with numbers only. It is written with columns as species and rows as stations. Be sure to save the name of the stations as rownames.
communaute <- BDD_sec[,-1]# remove the columns with data that are not numbers
rownames(communaute)=BDD_sec$TOPO # save the stations in rownames

#### iNEXT : computing Hill numbers and Hill curves.
# The iNEXT function needs the transposed table. q is the list of Hill numbers we need (0 = species richness, 1 is linked to shannon and 2 to simpson)
# We have to precise that the table contains abundance data (could be incidence).
iNEXT(t(communaute), q = c(0,1,2), datatype ="abundance")->test # test is a list of 3 elements (nbr de Hill)
test$DataInfo # summary of abundance, species richness, coverage and number of singletons, doubletons, tripletons, etc
test$iNextEst # estimation of coordinates for the rarefaction/accumulation curve. Completeness in X and hill number in Y (there is a list of coordinates for each hill order)
test$AsyEst # Asymptotic diversity estimates for each order. Sites are labeled with letters.

# estimateD is the rarefaction function. It works with Hill numbers 0, 1, 2.
# Transposed table is used. By default datatype ="abundance" (could be incidence)
# We have to precise which metric is used for completeness (coverage is better than abundance). Their definition of coverage is in reference papers.
# level = NULL by default (not written here). This way, the rarefaction use the smallest value of completeness among sites. Could be determined by user.
estimateD(t(communaute), base="coverage")->hill # m is the sample size for the reference level of completeness (we used coverage as completeness metric, thus sample size is not constant)
# SC = sample coverage (should be roughly equal to the minimum coverage among sites)
# the table says if the data is interpolated (rarefaction), observed (for sites with the lowest coverage) or extrapolated (should not happen as we chose to lower each station coverage to the lowest one)
#produce an ugly table, need some transformations (one station = 3 rows)
#values for hill numbers are given with upper and lower estimate (95% CI)
hillmatrix0= hill[hill[,"order"]==0,]
hillmatrix1= hill[hill[,"order"]==1,]
hillmatrix2= hill[hill[,"order"]==2,]
hillmatrix = hillmatrix0[,c(1,2,4)]
hillmatrix$SC0 = test$DataInfo$SC# retrieve the initial coverage from the table computed by iNEXT()
hillmatrix[,c("H0r","H0rU","H0rL")]=hillmatrix0[,c("qD","qD.UCL","qD.LCL")] #limite inférieur et supérieur de CI
hillmatrix[,c("H1r","H1rU","H1rL")]=hillmatrix1[,c("qD","qD.UCL","qD.LCL")]
hillmatrix[,c("H2r","H2rU","H2rL")]=hillmatrix2[,c("qD","qD.UCL","qD.LCL")]
rownames(hillmatrix)=rownames(communaute) # put the labels on the rows