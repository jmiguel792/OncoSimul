library(OncoSimulR)

##############################################################
## FIRST EXAMPLE ## -> I expect everything's normal
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
                      finalTime = 500,
                      mu = 1e-4,
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim1, show = "genotypes", col = c("green", "red", "yellow"))
##############################################################
## SECOND EXAMPLE ## -> just to see that exprtk get the value from muFactor expression. 
##                      I expect same as in the first example.
set.seed(2)
sim2 <- oncoSimulIndiv(fe,
                      model = "McFL", 
                      onlyCancer = FALSE, 
                      finalTime = 500,
                      mu = 1e-6,
                      muFactor = "100",
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim2, show = "genotypes", col = c("black", "red", "yellow"))
###################################################################
## THIRD EXAMPLE ## -> I expect that after 300 units of time appear all the genotypes
muexpression = "if(T>300) 10000; else 1;"
set.seed(2)
sim3 <- oncoSimulIndiv(fe,
                      model = "McFL", 
                      onlyCancer = FALSE, 
                      finalTime = 500,
                      mu = 1e-6,
                      muFactor = muexpression,
                      initSize = 5000, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = FALSE, 
                      errorHitWallTime = FALSE)

plot(sim3, show = "genotypes", col = c("green", "red", "yellow"))
#######################################################################
## FOURTH EXAMPLE ## -> this is more complex, but I actually expect something similar to the previous one
##                      I don't know if this specification alter significantly the resulting plot
muexpression = "if(T>300) 10000; else if(T>400) 100; else 1;"
set.seed(2)
sim4 <- oncoSimulIndiv(fe,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 500,
                       mu = 1e-6,
                       muFactor = muexpression,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(sim4, show = "genotypes", col = c("green", "red", "yellow"))
#######################################################################
## FIFTH EXAMPLE ## I will include an "and" just to show that we can use a mu expression following this configuration
muexpression = "if(T>300 and T<400) 10000; else if(T>400) 100; else 1;"
set.seed(2)
sim5 <- oncoSimulIndiv(fe,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 500,
                       mu = 1e-6,
                       muFactor = muexpression,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(sim5, show = "genotypes", col = c("green", "red", "yellow"))
####################################################################
## SIXTH EXAMPLE ## update after using relative frequencies
muexpression = "if(f_1>0.3) 100; else 1;"
set.seed(2)
sim6 <- oncoSimulIndiv(fe,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 500,
                       mu = 1e-6,
                       muFactor = muexpression,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

plot(sim6, show = "genotypes", col = c("green", "red", "yellow"))
####################################################################
####################################################################
sim1$pops.by.time[nrow(sim1$pops.by.time), ]
sim2$pops.by.time[nrow(sim2$pops.by.time), ]
sim3$pops.by.time[nrow(sim3$pops.by.time), ]
sim4$pops.by.time[nrow(sim4$pops.by.time), ]
sim5$pops.by.time[nrow(sim5$pops.by.time), ]
####################################################################
## SEVENTH EXAMPLE ## update after using absolute frequencies
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
##############################################################
library(OncoSimulR)

set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))

fe3 <- allFitnessEffects(noIntGenes = ni)

fm3 <- allMutatorEffects(epistasis = c("A" = 5,
                                       "B" = 10,
                                       "C" = 3,
                                       "A:C" = 70),
                         geneToModule = c("A" = "a1, a2",
                                          "B" = "b1, b2, b3",
                                          "C" = "c1, c2"))
set.seed(1) ## so that it is easy to reproduce
mue1 <- oncoSimulIndiv(fe3, muEF = fm3, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 500,
                       onlyCancer = FALSE)
plot(mue1)
################################################################
set.seed(1) ## for reproducibility
## 17 genes, 7 with no direct fitness effects
ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
               paste0("g", 1:10))

fe3 <- allFitnessEffects(noIntGenes = ni)

fm3 <- allMutatorEffects(epistasis = c("A" = 5,
                                       "B" = 10,
                                       "C" = 3,
                                       "A:C" = 70),
                         geneToModule = c("A" = "a1, a2",
                                          "B" = "b1, b2, b3",
                                          "C" = "c1, c2"))
set.seed(1) ## so that it is easy to reproduce
mue1 <- oncoSimulIndiv(fe3, muEF = fm3, 
                       mu = 1e-6,
                       muFactor = "if(T>200) 2; else 1;",
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 500,
                       onlyCancer = FALSE)
plot(mue1)
##############################################################
library(OncoSimulR)
sa <- 0.1
sb <- -0.2
sab <- 0.25
sac <- -0.1
sbc <- 0.25
sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A : C" = sac,
                                       "A:B" = sab,
                                       "-A:B:C" = sbc),
                         geneToModule = c(
                           "A" = "a1, a2",
                           "B" = "b",
                           "C" = "c"),
                         drvNames = c("a1", "a2", "b", "c"))
RNGkind("Mersenne-Twister")
set.seed(983)
ep1 <- oncoSimulIndiv(sv2, model = "McFL",
                      mu = 5e-6,
                      muFactor = "if(T>200) 100; else 1;",
                      sampleEvery = 0.025,
                      keepEvery = 0.5,
                      initSize = 2000,
                      finalTime = 3000,
                      onlyCancer = FALSE)
plot(ep1, show = "drivers", xlim = c(0, 1500),
     thinData = TRUE, thinData.keep = 0.5)
