library(OncoSimulR)

# lower to higher mutationRate
##############################################################
fl <- data.frame(
  Genotype = c("WT", "A", "B"),
  Fitness = c("1 + f_1*f_2",
              "1 + 1.5*f_*f_2",
              "1 + 1.5*f_*f_1"),
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
                      mu = 1e-4,
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
                      muFactor = "100",
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim2, show = "genotypes", col = c("black", "red", "yellow"))
###################################################################

muexpression = "if( f_1 > 0.3 ) 10000; else 1;"
#muexpression = "if(T>400 and T<600) 10000; else if(T>600) 100; else 1;"
set.seed(2)
sim3 <- oncoSimulIndiv(fe,
                      model = "McFL", 
                      onlyCancer = FALSE, 
                      finalTime = 1000,
                      mu = 1e-6,
                      muFactor = muexpression,
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim3, show = "genotypes", col = c("green", "red", "yellow"))
#######################################################################

# higher to lower mutationRate
set.seed(2)
sim4 <- oncoSimulIndiv(fe,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 1000,
                       mu = 1e-2,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(sim4, show = "genotypes", col = c("green", "red", "yellow"))
##############################################################

set.seed(2)
sim5 <- oncoSimulIndiv(fe,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 1000,
                       mu = 1e-2,
                       muFactor = c((1/1000),0),
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(sim5, show = "genotypes", col = c("green", "red", "yellow"))
###################################################################

set.seed(2)
sim6 <- oncoSimulIndiv(fe,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 1000,
                       mu = 1e-2,
                       muFactor = c((1/10000),200),
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(sim6, show = "genotypes", col = c("purple", "red", "yellow"))
####################################################################
sim1$pops.by.time[nrow(sim1$pops.by.time), ]
sim2$pops.by.time[nrow(sim2$pops.by.time), ]
sim3$pops.by.time[nrow(sim3$pops.by.time), ]
sim4$pops.by.time[nrow(sim4$pops.by.time), ]
sim5$pops.by.time[nrow(sim5$pops.by.time), ]
sim6$pops.by.time[nrow(sim6$pops.by.time), ]
####################################################################
library(OncoSimulR)
rar3 <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                   Fitness = c("1",
                               "1.1 + .3*(n_2/N)", #(n_2/N)
                               "1.2 + .4*(n_1/N)", #(n_1/N)
                               "1.0 + .5 * ( n_1 > 20)"),
                   stringsAsFactors = FALSE)

afear3 <- allFitnessEffects(genotFitness = rar3, 
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")

set.seed(1)
tmp3 <- oncoSimulIndiv(afear3, 
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 100,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(tmp3, show = "genotypes")
#############################################################

muexpression = "if( n_1 > 10 ) 100; else 1;"
set.seed(1)
tmp4 <- oncoSimulIndiv(afear3, 
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 100,
                       mu = 1e-4,
                       muFactor = muexpression,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(tmp4, show = "genotypes")