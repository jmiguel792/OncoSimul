library(OncoSimulR)

###########################################################
r1fd <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("1",
                               "1.4 + 1*(f_2)",
                               "1.4 + 1*(f_1)",
                               "1.6 + f_1 + f_2"),
                   stringsAsFactors = FALSE)

afe4 <- allFitnessEffects(genotFitness = r1fd, 
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

mtfd <- allMutatorEffects(epistasis = c("A" = 0.1,"B" = 10))

############################################################

allg <- OncoSimulR:::generateAllGenotypes(fitnessEffects = afe4, order = FALSE, max = 256)
full2mutator_ <- OncoSimulR:::matchGeneIDs(mtfd, afe4)$Reduced

##ERROR## evalRGenotypeAndMut -> new_restrict.cpp
#############################################################
Rcpp::sourceCpp('OncoSimul-time/OncoSimulR/src/new_restrict.cpp')
allf <- t(vapply(allg$genotNums,
                 function(x) evalRGenotypeAndMut(x,
                                                 rFE = afe4,
                                                 muEF= mtfd,
                                                 full2mutator_ = full2mutator_,
                                                 verbose = FALSE,
                                                 prodNeg = FALSE,
                                                 currentTime = 0),
                 c(1.1, 2.2)))
##############################################################
