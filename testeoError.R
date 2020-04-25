library(OncoSimulR)

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

evalAllGenotypes(allFitnessEffects(genotFitness = r1fd, 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(10, 20, 30, 40)))

## Fitness is wrong
evalAllGenotypesFitAndMut(allFitnessEffects(genotFitness = r1fd, 
                                            frequencyDependentFitness = TRUE,
                                            frequencyType = "rel",
                                            spPopSizes = c(10, 20, 30, 40)),
                          mutatorEffects = mtfd)