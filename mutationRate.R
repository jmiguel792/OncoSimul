library(OncoSimulR)
# low mutationRate to high mutationRate

##############################################################
fl <- data.frame(
  Genotype = c("WT", "A", "B"),
  Fitness = c("1 + f_1*f_2", #WT 1 + f_1*f_2
              "1 + 1.5*f_*f_2", #A 1 + 1.5*f_*f_2
              "1 + 1.5*f_*f_1"), #B 1 + 1.5*f_*f_1
  stringsAsFactors = FALSE
)

fe <- allFitnessEffects(genotFitness = fl, 
                        frequencyDependentFitness = TRUE,
                        frequencyType = "rel")

set.seed(2)
sim1 <- oncoSimulIndiv(fe,
                      model = "McFL", 
                      onlyCancer = FALSE, 
                      finalTime = 1000,
                      mu = 1e-6,
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim1, show = "genotypes", col = c("green", "red", "yellow"))                    
##############################################################

set.seed(2)
sim2 <- oncoSimulIndiv(fe,
                      model = "McFL", 
                      onlyCancer = FALSE, 
                      finalTime = 1000,
                      mu = 1e-6,
                      muFactor = c(100,0),
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim2, show = "genotypes", col = c("green", "red", "yellow"))
###################################################################

set.seed(2)
sim3 <- oncoSimulIndiv(fe,
                      model = "McFL", 
                      onlyCancer = FALSE, 
                      finalTime = 1000,
                      mu = 1e-6,
                      muFactor = c(10000,400),
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim3, show = "genotypes", col = c("purple", "red", "yellow"))
#######################################################################