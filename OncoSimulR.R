## ----setup, include=FALSE------------------------------------------------
## use collapse for bookdown, to collapse all the source and output
## blocks from one code chunk into a single block
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE)
options(width = 70)
require(BiocStyle)
require(pander)


## ----frelats, eval=TRUE,echo=FALSE, fig.cap="Relationships between the main functions in OncoSimulR."----
knitr::include_graphics("relfunct.png")


## ----firstload-----------------------------------------------------------
## Load the package
library(OncoSimulR) 


## ---- echo=FALSE---------------------------------------------------------
set.seed(2)
RNGkind("L'Ecuyer-CMRG")


## ----ex-dag-inf----------------------------------------------------------
## For reproducibility
set.seed(2)
RNGkind("L'Ecuyer-CMRG")

## Simulate a DAG
g1 <- simOGraph(4, out = "rT")

## Simulate 10 evolutionary trajectories
s1 <- oncoSimulPop(10, allFitnessEffects(g1, drvNames = 1:4),
                   mc.cores = 2, ## adapt to your hardware
                   seed = NULL) ## for reproducibility of vignette

## Sample those data uniformly, and add noise
d1 <- samplePop(s1, timeSample = "unif", propError = 0.1)

## You would now run the appropriate inferential method and
## compare observed and true. For example

## require(Oncotree)
## fit1 <- oncotree.fit(d1)

## Now, you'd compare fitted and original. This is well beyond 
## the scope of this document (and OncoSimulR itself).



## ----hidden-rng-exochs, echo = FALSE-------------------------------------
set.seed(NULL)


## ----hiddenochs, echo=FALSE----------------------------------------------
set.seed(2)
RNGkind("L'Ecuyer-CMRG")


## ----exochs--------------------------------------------------------------
## For reproducibility
set.seed(2)
RNGkind("L'Ecuyer-CMRG")

## Specify fitness effects. 

## Numeric values arbitrary, but set the intermediate genotype en
## route to ui as mildly deleterious so there is a valley.

## As in Ochs and Desai, the ui and uv genotypes
## can never appear. 

u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf

od <- allFitnessEffects(
    epistasis = c("u" = u,  "u:i" = ui,
                  "u:v" = uv, "i" = i,
                  "v:-i" = -Inf, "v:i" = vi))

## For the sake of extending this example, also turn i into a
## mutator gene

odm <- allMutatorEffects(noIntGenes = c("i" = 50))

## How do mutation and fitness look like for each genotype?
evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)


## ----exochsb-------------------------------------------------------------
## Set a small initSize, as o.w. unlikely to pass the valley
initS <- 10
## The number of replicates is tiny, 10, for the sake of speed
## of creation of the vignette
od_sim <- oncoSimulPop(10, od, muEF = odm,
                       fixation = c("u", "i, v"), initSize = initS,
                       model = "McFL",
                       mu = 1e-4, detectionDrivers = NA, 
					   finalTime = NA,
                       detectionSize = NA, detectionProb = NA,
                       onlyCancer = TRUE, 
					   mc.cores = 2, ## adapt to your hardware
                       seed = NULL) ## for reproducibility
## What is the frequency of each final genotype?
sampledGenotypes(samplePop(od_sim))


## ----hidden-rng-exochs33, echo = FALSE-----------------------------------
set.seed(NULL)


## ----hiddenrng0szen, echo=FALSE------------------------------------------
set.seed(7)
RNGkind("L'Ecuyer-CMRG")


## ----exszendro-----------------------------------------------------------
## For reproducibility
set.seed(7)
RNGkind("L'Ecuyer-CMRG")


## Repeat the following loop for different combinations of whatever
## interests you, such as number of genes, or distribution of the
## c and sd (which affect how rugged the landscape is), or 
## reference genotype, or evolutionary model, or stopping criterion, 
## or sampling procedure, or ...

##  Generate a random fitness landscape, from the Rough Mount
##  Fuji model, with g genes, and c ("slope" constant) and
##  reference chosen randomly (reference is random by default and 
##  thus not specified below). Require a minimal number of  
##  accessible genotypes

g <- 6
c <- runif(1, 1/5, 5)
rl <- rfitness(g, c = c, min_accessible_genotypes = g)

## Plot it if you want; commented here as it takes long for a
## vignette

## plot(rl)

## Obtain landscape measures from MAGELLAN. Export to MAGELLAN and
## call your own copy of MAGELLAN's binary
to_Magellan(rl, file = "rl1.txt")

## or use the binary copy provided with OncoSimulR
## see also below.
Magellan_stats(rl)

## Simulate evolution in that landscape many times (here just 10)
simulrl <- oncoSimulPop(10, allFitnessEffects(genotFitness = rl),
                        keepPhylog = TRUE, keepEvery = 1,
                        initSize = 4000,
                        seed = NULL, ## for reproducibility
                        mc.cores = 2) ## adapt to your hardware

## Obtain measures of evolutionary predictability
diversityLOD(LOD(simulrl))
diversityPOM(POM(simulrl))
sampledGenotypes(samplePop(simulrl, typeSample = "whole"))


## ----hidden-rng-exszend, echo = FALSE------------------------------------
set.seed(NULL)


## ----ex-tomlin1----------------------------------------------------------
sd <- 0.1 ## fitness effect of drivers
sm <- 0 ## fitness effect of mutator
nd <- 20 ## number of drivers
nm <- 5  ## number of mutators
mut <- 10 ## mutator effect

fitnessGenesVector <- c(rep(sd, nd), rep(sm, nm))
names(fitnessGenesVector) <- 1:(nd + nm)
mutatorGenesVector <- rep(mut, nm)
names(mutatorGenesVector) <- (nd + 1):(nd + nm)

ft <- allFitnessEffects(noIntGenes = fitnessGenesVector,
                        drvNames = 1:nd)
mt <- allMutatorEffects(noIntGenes = mutatorGenesVector)



## ----hiddentom, echo=FALSE-----------------------------------------------
set.seed(2)
RNGkind("L'Ecuyer-CMRG")


## ----ex-tomlin2----------------------------------------------------------
## For reproducibility
set.seed(2)
RNGkind("L'Ecuyer-CMRG")


ddr <- 4
st <- oncoSimulPop(4, ft, muEF = mt,
                   detectionDrivers = ddr,
                   finalTime = NA,
                   detectionSize = NA,
                   detectionProb = NA,
                   onlyCancer = TRUE,
                   keepEvery = NA, 
                   mc.cores = 2, ## adapt to your hardware
                   seed = NULL) ## for reproducibility

## How long did it take to reach cancer?
unlist(lapply(st, function(x) x$FinalTime))


## ----hidden-rng-tom, echo = FALSE----------------------------------------
set.seed(NULL)


## ----exusagebau----------------------------------------------------------
K <- 4
sp <- 1e-5
sdp <- 0.015
sdplus <- 0.05
sdminus <- 0.1

cnt <- (1 + sdplus)/(1 + sdminus)
prod_cnt <- cnt - 1
bauer <- data.frame(parent = c("Root", rep("D", K)),
                    child = c("D", paste0("s", 1:K)),
                    s = c(prod_cnt, rep(sdp, K)),
                    sh = c(0, rep(sp, K)),
                    typeDep = "MN")
fbauer <- allFitnessEffects(bauer)
(b1 <- evalAllGenotypes(fbauer, order = FALSE, addwt = TRUE))

## How does the fitness landscape look like?
plot(b1, use_ggrepel = TRUE) ## avoid overlapping labels


## ----hiddenbau, echo=FALSE-----------------------------------------------
set.seed(2)
RNGkind("L'Ecuyer-CMRG")


## ----exusagebau2---------------------------------------------------------
## For reproducibility
set.seed(2)
RNGkind("L'Ecuyer-CMRG")

totalpops <- 5
initSize <- 100
sb1 <- oncoSimulPop(totalpops, fbauer, model = "Exp",
                    initSize = initSize,
                    onlyCancer = FALSE, 
					mc.cores = 2, ## adapt to your hardware
                    seed = NULL) ## for reproducibility
					
## What proportion of the simulations reach 4x initSize?
sum(summary(sb1)[, "TotalPopSize"] > (4 * initSize))/totalpops


## ----hidden-rng-exbau, echo = FALSE--------------------------------------
set.seed(NULL)


## ----hiddenbau22, echo=FALSE---------------------------------------------
set.seed(2)
RNGkind("L'Ecuyer-CMRG")


## ----hhhhbbbb22----------------------------------------------------------

totalpops <- 5
initSize <- 100
sb2 <- oncoSimulPop(totalpops, fbauer, model = "Exp",
                    initSize = initSize,
                    onlyCancer = TRUE,
					detectionSize = 10 * initSize,
					mc.cores = 2, ## adapt to your hardware
                    seed = NULL) ## for reproducibility
				
## How long did it take to reach cancer?
unlist(lapply(sb2, function(x) x$FinalTime))


## ----hidden-rng-exbau22, echo = FALSE------------------------------------
set.seed(NULL)


## ----oex1intro-----------------------------------------------------------
## Order effects involving three genes.

## Genotype "D, M" has different fitness effects
## depending on whether M or D mutated first.
## Ditto for genotype "F, D, M".

## Meaning of specification: X > Y means
## that X is mutated before Y.


o3 <- allFitnessEffects(orderEffects = c(
                            "F > D > M" = -0.3,
                            "D > F > M" = 0.4,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.5))

## With the above specification, let's double check
## the fitness of the possible genotypes

(oeag <- evalAllGenotypes(o3, addwt = TRUE, order = TRUE))



## ----hiddoef, echo=FALSE-------------------------------------------------
set.seed(2)
RNGkind("L'Ecuyer-CMRG")


## ----exusageoe2----------------------------------------------------------
## For reproducibility
set.seed(2)
RNGkind("L'Ecuyer-CMRG")

totalpops <- 5
soe1 <- oncoSimulPop(totalpops, o3, model = "Exp",
                    initSize = 500,
                    onlyCancer = FALSE,
					mc.cores = 2, ## adapt to your hardware
                    seed = NULL) ## for reproducibility
					
## What proportion of the simulations do not end up extinct?
sum(summary(soe1)[, "TotalPopSize"] > 0)/totalpops



## ----hidden-rng-exoef, echo = FALSE--------------------------------------
set.seed(NULL)


## ---- results="hide",message=FALSE, echo=TRUE, include=TRUE--------------
library(OncoSimulR)
library(graph)
library(igraph)
igraph_options(vertex.color = "SkyBlue2")


## ---- echo=FALSE, results='hide'-----------------------------------------
options(width = 68)


## ------------------------------------------------------------------------
packageVersion("OncoSimulR")


## ---- fig.width=6.5, fig.height=10---------------------------------------
## 1. Fitness effects: here we specify an 
##    epistatic model with modules.
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
evalAllGenotypes(sv2, addwt = TRUE)

## 2. Simulate the data. Here we use the "McFL" model and set
##    explicitly parameters for mutation rate, initial size, size
##    of the population that will end the simulations, etc

RNGkind("Mersenne-Twister")
set.seed(983)
ep1 <- oncoSimulIndiv(sv2, model = "McFL",
                      mu = 5e-6,
                      sampleEvery = 0.025,
                      keepEvery = 0.5,
                      initSize = 2000,
                      finalTime = 3000,
                      onlyCancer = FALSE)


## ----iep1x1,fig.width=6.5, fig.height=4.5, fig.cap="Plot of drivers of an epistasis simulation."----
## 3. We will not analyze those data any further. We will only plot
## them.  For the sake of a small plot, we thin the data.
plot(ep1, show = "drivers", xlim = c(0, 1500),
     thinData = TRUE, thinData.keep = 0.5)


## ----fepancr1, fig.width=5-----------------------------------------------
## 1. Fitness effects: 
pancr <- allFitnessEffects(
    data.frame(parent = c("Root", rep("KRAS", 4), 
                   "SMAD4", "CDNK2A", 
                   "TP53", "TP53", "MLL3"),
               child = c("KRAS","SMAD4", "CDNK2A", 
                   "TP53", "MLL3",
                   rep("PXDN", 3), rep("TGFBR2", 2)),
               s = 0.1,
               sh = -0.9,
               typeDep = "MN"),
    drvNames = c("KRAS", "SMAD4", "CDNK2A", "TP53", 
	             "MLL3", "TGFBR2", "PXDN"))


## ----figfpancr1, fig.width=5, fig.cap="Plot of DAG corresponding to fitnessEffects object."----
## Plot the DAG of the fitnessEffects object
plot(pancr)


## ------------------------------------------------------------------------
## 2. Simulate from it. We change several possible options. 

set.seed(4) ## Fix the seed, so we can repeat it
ep2 <- oncoSimulIndiv(pancr, model = "McFL",
                     mu = 1e-6,
                     sampleEvery = 0.02,
                     keepEvery = 1,
                     initSize = 1000,
                     finalTime = 10000,
                     onlyCancer = FALSE)


## ----iep2x2, fig.width=6.5, fig.height=5, fig.cap= "Plot of genotypes of a simulation from a DAG."----
## 3. What genotypes and drivers we get? And play with limits
##    to show only parts of the data. We also aggressively thin
##    the data.
par(cex = 0.7)
plot(ep2, show = "genotypes", xlim = c(1000, 8000), 
     ylim = c(0, 2400),
     thinData = TRUE, thinData.keep = 0.03)


## ------------------------------------------------------------------------
citation("OncoSimulR")


## ----colnames_benchmarks, echo = FALSE, eval = TRUE----------------------

data(benchmark_1)
data(benchmark_1_0.05)
data(benchmark_2)
data(benchmark_3)

colnames(benchmark_1)[
    match(c(
	"time_per_simul",
    "size_mb_per_simul", "NumClones.Median", "NumIter.Median",
	"FinalTime.Median", "TotalPopSize.Median", "TotalPopSize.Mean",
	"TotalPopSize.Max.", "keepEvery",  "Attempts.Median",
	"Attempts.Mean", "Attempts.Max.",
	"PDBaseline", "n2", "onlyCancer"),
	 colnames(benchmark_1)
	)] <- c("Elapsed Time, average per simulation (s)",
	              "Object Size, average per simulation (MB)",
				  "Number of Clones, median",
				  "Number of Iterations, median",
				  "Final Time, median",
				  "Total Population Size, median",
				   "Total Population Size, mean",
				  "Total Population Size, max.",
				  "keepEvery",
				  "Attempts until Cancer, median",
				  "Attempts until Cancer, mean",
				  "Attempts until Cancer, max.",
				  "PDBaseline", "n2", "onlyCancer"
				  )
				  
	
colnames(benchmark_1_0.05)[
    match(c("time_per_simul",
    "size_mb_per_simul", "NumClones.Median", "NumIter.Median",
	"FinalTime.Median", "TotalPopSize.Median", "TotalPopSize.Mean", 
	"TotalPopSize.Max.",
	"keepEvery",
	"PDBaseline", "n2", "onlyCancer", "Attempts.Median"),
	colnames(benchmark_1_0.05))] <- c("Elapsed Time, average per simulation (s)",
	              "Object Size, average per simulation (MB)",
				  "Number of Clones, median",
				  "Number of Iterations, median",
				  "Final Time, median",
				  "Total Population Size, median",
				  "Total Population Size, mean",
				  "Total Population Size, max.",
				  "keepEvery",
				  "PDBaseline", "n2", "onlyCancer",
				  "Attempts until Cancer, median"
				  )


colnames(benchmark_2)[match(c("Model", "fitness", "time_per_simul",
    "size_mb_per_simul", "NumClones.Median", "NumIter.Median",
	"FinalTime.Median", "TotalPopSize.Median", "TotalPopSize.Mean", 
	"TotalPopSize.Max."), colnames(benchmark_2))] <-  c("Model",
				  "Fitness",
	"Elapsed Time, average per simulation (s)",
	              "Object Size, average per simulation (MB)",
				  "Number of Clones, median",
				  "Number of Iterations, median",
				  "Final Time, median",
				  "Total Population Size, median",
				  "Total Population Size, mean",
				  "Total Population Size, max."
				  )	
				  
colnames(benchmark_3)[match(c("Model", "fitness", "time_per_simul",
    "size_mb_per_simul", "NumClones.Median", "NumIter.Median",
	"FinalTime.Median", "TotalPopSize.Median", "TotalPopSize.Mean", 
	"TotalPopSize.Max."), colnames(benchmark_3))] <-  c("Model",
				  "Fitness",
	"Elapsed Time, average per simulation (s)",
	              "Object Size, average per simulation (MB)",
				  "Number of Clones, median",
				  "Number of Iterations, median",
				  "Final Time, median",
				  "Total Population Size, median",
				  "Total Population Size, mean",
				  "Total Population Size, max."
				  )					  


## ----timing1, eval=FALSE-------------------------------------------------
## ## Specify fitness
## pancr <- allFitnessEffects(
##     data.frame(parent = c("Root", rep("KRAS", 4),
##                    "SMAD4", "CDNK2A",
##                    "TP53", "TP53", "MLL3"),
##                child = c("KRAS","SMAD4", "CDNK2A",
##                    "TP53", "MLL3",
##                    rep("PXDN", 3), rep("TGFBR2", 2)),
##                s = 0.1,
##                sh = -0.9,
##                typeDep = "MN"),
##     drvNames = c("KRAS", "SMAD4", "CDNK2A", "TP53",
## 	             "MLL3", "TGFBR2", "PXDN"))
## 
## Nindiv <- 100 ## Number of simulations run.
##               ## Increase this number to decrease sampling variation
## 
## ## keepEvery = 1
## t_exp1 <- system.time(
##     exp1 <- oncoSimulPop(Nindiv, pancr,
##                             detectionProb = "default",
##                             detectionSize = NA,
##                             detectionDrivers = NA,
##                             finalTime = NA,
##                             keepEvery = 1,
##                             model = "Exp",
##                             mc.cores = 1))["elapsed"]/Nindiv
## 
## 
## t_mc1 <- system.time(
##     mc1 <- oncoSimulPop(Nindiv, pancr,
##                            detectionProb = "default",
##                            detectionSize = NA,
##                            detectionDrivers = NA,
##                            finalTime = NA,
##                            keepEvery = 1,
##                            model = "McFL",
##                            mc.cores = 1))["elapsed"]/Nindiv
## 
## ## keepEvery = NA
## t_exp2 <- system.time(
##     exp2 <- oncoSimulPop(Nindiv, pancr,
##                             detectionProb = "default",
##                             detectionSize = NA,
##                             detectionDrivers = NA,
##                             finalTime = NA,
##                             keepEvery = NA,
##                             model = "Exp",
##                             mc.cores = 1))["elapsed"]/Nindiv
## 
## 
## t_mc2 <- system.time(
##     mc2 <- oncoSimulPop(Nindiv, pancr,
##                            detectionProb = "default",
##                            detectionSize = NA,
##                            detectionDrivers = NA,
##                            finalTime = NA,
##                            keepEvery = NA,
##                            model = "McFL",
##                            mc.cores = 1))["elapsed"]/Nindiv
## 
## 


## ---- eval=FALSE---------------------------------------------------------
## cat("\n\n\n t_exp1 = ", t_exp1, "\n")
## object.size(exp1)/(Nindiv * 1024^2)
## cat("\n\n")
## summary(unlist(lapply(exp1, "[[", "NumClones")))
## summary(unlist(lapply(exp1, "[[", "NumIter")))
## summary(unlist(lapply(exp1, "[[", "FinalTime")))
## summary(unlist(lapply(exp1, "[[", "TotalPopSize")))


## ----bench1, eval=TRUE, echo = FALSE-------------------------------------

panderOptions('table.split.table', 99999999)
panderOptions('table.split.cells', 900)  ## For HTML
## panderOptions('table.split.cells', 8) ## For PDF

set.alignment('right')
panderOptions('round', 2)
panderOptions('big.mark', ',')
panderOptions('digits', 2)
				          
pander(benchmark_1[1:4, c("Elapsed Time, average per simulation (s)", 
 	              "Object Size, average per simulation (MB)",
 				  "Number of Clones, median",
 				  "Number of Iterations, median",
 				  "Final Time, median",
 				  "Total Population Size, median",
 				  "Total Population Size, max.",
 				  "keepEvery")],
				  justify = c('left', rep('right', 8)), ##  o.w. hlines not right
				  ## caption = "\\label{tab:bench1}Benchmarks of Exp and McFL  models using the default `detectionProb` with two settings of `keepEvery`."
				  )


## ----timing2, eval = FALSE-----------------------------------------------
## t_exp3 <- system.time(
##     exp3 <- oncoSimulPop(Nindiv, pancr,
##                             detectionProb = c(PDBaseline = 5e4,
##                                               p2 = 0.1, n2 = 5e5,
##                                               checkSizePEvery = 20),
##                             detectionSize = NA,
##                             detectionDrivers = NA,
##                             finalTime = NA,
##                             keepEvery = 1,
##                             model = "Exp",
##                             mc.cores = 1))["elapsed"]/Nindiv
## 
## t_exp4 <- system.time(
##     exp4 <- oncoSimulPop(Nindiv, pancr,
##                             detectionProb = c(PDBaseline = 5e4,
##                                               p2 = 0.1, n2 = 5e5,
##                                               checkSizePEvery = 20),
##                             detectionSize = NA,
##                             detectionDrivers = NA,
##                             finalTime = NA,
##                             keepEvery = NA,
##                             model = "Exp",
##                             mc.cores = 1))["elapsed"]/Nindiv
## 
## 
## 
## t_exp5 <- system.time(
##     exp5 <- oncoSimulPop(Nindiv, pancr,
##                             detectionProb = c(PDBaseline = 5e5,
##                                               p2 = 0.1, n2 = 5e7),
##                             detectionSize = NA,
##                             detectionDrivers = NA,
##                             finalTime = NA,
##                             keepEvery = 1,
##                             model = "Exp",
##                             mc.cores = 1))["elapsed"]/Nindiv
## 
## t_exp6 <- system.time(
##     exp6 <- oncoSimulPop(Nindiv, pancr,
##                             detectionProb = c(PDBaseline = 5e5,
##                                               p2 = 0.1, n2 = 5e7),
##                             detectionSize = NA,
##                             detectionDrivers = NA,
##                             finalTime = NA,
##                             keepEvery = NA,
##                             model = "Exp",
##                             mc.cores = 1))["elapsed"]/Nindiv
## 


## ----bench1b, eval=TRUE, echo = FALSE------------------------------------
panderOptions('table.split.table', 99999999)
panderOptions('table.split.cells', 900)  ## For HTML
## panderOptions('table.split.cells', 8) ## For PDF

set.alignment('right')
panderOptions('round', 2)
panderOptions('big.mark', ',')
panderOptions('digits', 2)

pander(benchmark_1[5:8, c("Elapsed Time, average per simulation (s)",
 	              "Object Size, average per simulation (MB)",
 				  "Number of Clones, median",
 				  "Number of Iterations, median",
 				  "Final Time, median",
 				  "Total Population Size, median",
 				  "Total Population Size, max.",
 				  "keepEvery",
				  "PDBaseline",
				  "n2")], 
				  justify = c('left', rep('right', 10)), ##  o.w. hlines not right
## 				  round = c(rep(2, 3), rep(0, 7)),
## 				  digits = c(rep(2, 3), rep(1, 7)),
	  ## caption = "\\label{tab:bench1b}Benchmarks of Exp and McFL models modifying the default `detectionProb` with two settings of `keepEvery`."
    )



## ----bench1c, eval=TRUE, echo = FALSE------------------------------------
panderOptions('table.split.table', 99999999)
panderOptions('table.split.cells', 900)  ## For HTML
## panderOptions('table.split.cells', 12) ## For PDF
set.alignment('right')
panderOptions('round', 2)
panderOptions('big.mark', ',')
panderOptions('digits', 2)

pander(benchmark_1[1:8, c(
"Attempts until Cancer, median", 
"Attempts until Cancer, mean", 
"Attempts until Cancer, max.", 
				  "PDBaseline",
				  "n2")], 
				  justify = c('left', rep('right', 5)), ##  o.w. hlines not right
## 				  round = c(rep(2, 3), rep(0, 7)),
## 				  digits = c(rep(2, 3), rep(1, 7)),
	  ## caption = "\\label{tab:bench1c}Number of attempts until cancer."
    )
## ## data(benchmark_1)
## knitr::kable(benchmark_1[1:8, c("Attempts.Median",
##                                 "PDBaseline", "n2"), drop = FALSE], 
##     booktabs = TRUE,
## 	row.names = TRUE,
## 	col.names = c("Attempts until cancer", "PDBaseline", "n2"),
##     caption = "Median number of attempts until cancer.", 
## 	align = "r")
	


## ----bench1d, eval=TRUE, echo = FALSE------------------------------------
panderOptions('table.split.table', 99999999)
panderOptions('table.split.cells', 900)  ## For HTML
## panderOptions('table.split.cells', 8) ## For PDF
panderOptions('table.split.cells', 15) ## does not fit otherwise
set.alignment('right')
panderOptions('round', 3)

pander(benchmark_1[9:16, 
    c("Elapsed Time, average per simulation (s)",
 	              "Object Size, average per simulation (MB)",
 				  "Number of Clones, median",
 				  "Number of Iterations, median",
 				  "Final Time, median",
 				  "Total Population Size, median",
				  "Total Population Size, mean",
 				  "Total Population Size, max.",
 				  "keepEvery",
				  "PDBaseline",
				  "n2")],
				  justify = c('left', rep('right', 11)), ##  o.w. hlines not right
## caption = "\\label{tab:timing3} Benchmarks of models in Table \\@ref(tab:bench1) and \\@ref(tab:bench1b) when run with `onlyCancer = FALSE`."
				  )	
	


## ----bench1dx0, eval=TRUE, echo = FALSE----------------------------------
panderOptions('table.split.table', 99999999)
## panderOptions('table.split.cells', 900)  ## For HTML
panderOptions('table.split.cells', 19)

set.alignment('right') 
panderOptions('round', 3)
	
pander(benchmark_1[ , c("Elapsed Time, average per simulation (s)",
 	              "Object Size, average per simulation (MB)", 
				  "Number of Clones, median", 
				  "Number of Iterations, median", 
				  "Final Time, median", "Total Population Size, median", 
				  "Total Population Size, mean", "Total Population Size, max.",
 	              "keepEvery", "PDBaseline", "n2", "onlyCancer")], 
				  justify = c('left', rep('right', 12)), ##  o.w. hlines not right
				  ## caption = "\\label{tab:allr1bck}Benchmarks of all models in Tables \\@ref(tab:bench1), \\@ref(tab:bench1b),  and \\@ref(tab:timing3)."  
				  )  


## ----bench1dx, eval=TRUE, echo = FALSE-----------------------------------
## data(benchmark_1_0.05)
## knitr::kable(benchmark_1_0.05[, c("time_per_simul",
##     "size_mb_per_simul", "NumClones.Median", "NumIter.Median",
## 	"FinalTime.Median", "TotalPopSize.Median", "TotalPopSize.Mean", 
## 	"TotalPopSize.Max.",
## 	"keepEvery",
## 	"PDBaseline", "n2", "onlyCancer")], 
##     booktabs = TRUE,
## 	col.names = c("Elapsed Time, average per simulation (s)",
## 	              "Object Size, average per simulation (MB)",
## 				  "Number of Clones, median",
## 				  "Number of Iterations, median",
## 				  "Final Time, median",
## 				  "Total Population Size, median",
## 				  "Total Population Size, mean",
## 				  "Total Population Size, max.",				  
## 				  "keepEvery",
## 				  "PDBaseline", "n2", "onlyCancer"
## 				  ),
## ##    caption = "Benchmarks of models in Table \@ref(tab:bench1) and
## ##   \@ref(tab:bench1b) when run with `onlyCancer = FALSE`", 
## 	align = "c")
	
panderOptions('table.split.table', 99999999)
## panderOptions('table.split.cells', 900)  ## For HTML
panderOptions('table.split.cells', 19)

set.alignment('right') 
panderOptions('round', 3)
	
pander(benchmark_1_0.05[ , c("Elapsed Time, average per simulation (s)",
 	              "Object Size, average per simulation (MB)", 
				  "Number of Clones, median", 
				  "Number of Iterations, median", 
				  "Final Time, median", 
				  "Total Population Size, median", 
				  "Total Population Size, mean", "Total Population Size, max.",
 	              "keepEvery", "PDBaseline", "n2", "onlyCancer")], 
				  justify = c('left', rep('right', 12)), ##  o.w. hlines not right
 	              ## caption = "\\label{tab:timing3xf}Benchmarks of all models in Table \\@ref(tab:allr1bck) using $s=0.05$ (instead of $s=0.1$)."  
)  
				  


## ----fitusualb, echo = TRUE, eval = FALSE--------------------------------
## pancr <- allFitnessEffects(
##     data.frame(parent = c("Root", rep("KRAS", 4),
##                    "SMAD4", "CDNK2A",
##                    "TP53", "TP53", "MLL3"),
##                child = c("KRAS","SMAD4", "CDNK2A",
##                    "TP53", "MLL3",
##                    rep("PXDN", 3), rep("TGFBR2", 2)),
##                s = 0.1,
##                sh = -0.9,
##                typeDep = "MN"),
##     drvNames = c("KRAS", "SMAD4", "CDNK2A", "TP53",
## 	             "MLL3", "TGFBR2", "PXDN"))
## 
## 
## ## Random fitness landscape with 6 genes
## ## At least 50 accessible genotypes
## rfl6 <- rfitness(6, min_accessible_genotypes = 50)
## attributes(rfl6)$accessible_genotypes ## How many accessible
## rf6 <- allFitnessEffects(genotFitness = rfl6)
## 
## 
## ## Random fitness landscape with 12 genes
## ## At least 200 accessible genotypes
## rfl12 <- rfitness(12, min_accessible_genotypes = 200)
## attributes(rfl12)$accessible_genotypes ## How many accessible
## rf12 <- allFitnessEffects(genotFitness = rfl12)
## 
## 
## 
## 
## ## Independent genes; positive fitness from exponential distribution
## ## with mean around 0.1, and negative from exponential with mean
## ## around -0.02. Half of genes positive fitness effects, half
## ## negative.
## 
## ng <- 200 re_200 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10),
## -rexp(ng/2, 50)))
## 
## ng <- 500
## re_500 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10),
##                                            -rexp(ng/2, 50)))
## 
## ng <- 2000
## re_2000 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10),
##                                             -rexp(ng/2, 50)))
## 
## ng <- 4000
## re_4000 <- allFitnessEffects(noIntGenes = c(rexp(ng/2, 10),
##                                             -rexp(ng/2, 50)))
## 


## ----exp-usual-r, eval = FALSE, echo = TRUE------------------------------
## 
## oncoSimulPop(Nindiv,
##             fitness,
##             detectionProb = NA,
##             detectionSize = 1e6,
##             initSize = 500,
##             detectionDrivers = NA,
##             keepPhylog = TRUE,
##             model = "Exp",
##             errorHitWallTime = FALSE,
##             errorHitMaxTries = FALSE,
##             finalTime = 5000,
##             onlyCancer = FALSE,
##             mc.cores = 1,
##             sampleEvery = 0.5,
## 			keepEvery = 1)


## ----mc-usual-r, eval = FALSE, echo = TRUE-------------------------------
## initSize <- 1000
## oncoSimulPop(Nindiv,
##               fitness,
##                detectionProb = c(
##                    PDBaseline = 1.4 * initSize,
##                    n2 = 2 * initSize,
##                    p2 = 0.1,
##                    checkSizePEvery = 4),
##                initSize = initSize,
##                detectionSize = NA,
##                detectionDrivers = NA,
##                keepPhylog = TRUE,
##                model = "McFL",
##                errorHitWallTime = FALSE,
##                errorHitMaxTries = FALSE,
##                finalTime = 5000,
##                max.wall.time = 10,
##                onlyCancer = FALSE,
##                mc.cores = 1,
## 			   keepEvery = 1)
## 


## ----benchustable, eval=TRUE, echo = FALSE-------------------------------
## data(benchmark_2)

## knitr::kable(benchmark_2[, c("Model", "fitness", "time_per_simul",
##     "size_mb_per_simul", "NumClones.Median", "NumIter.Median",
## 	"FinalTime.Median", "TotalPopSize.Median", "TotalPopSize.Mean", 
## 	"TotalPopSize.Max.")], 
##     booktabs = TRUE,
## 	col.names = c("Model",
## 				  "Fitness",
## 	"Elapsed Time, average per simulation (s)",
## 	              "Object Size, average per simulation (MB)",
## 				  "Number of Clones, median",
## 				  "Number of Iterations, median",
## 				  "Final Time, median",
## 				  "Total Population Size, median",
## 				  "Total Population Size, mean",
## 				  "Total Population Size, max."
## 				  ),
## 	align = "c")

panderOptions('table.split.table', 99999999)
panderOptions('table.split.cells', 900)  ## For HTML
## panderOptions('table.split.cells', 8) ## For PDF

## set.alignment('right', row.names = 'center')
panderOptions('table.alignment.default', 'right')

panderOptions('round', 3)

pander(benchmark_2[ , c(
    "Model", "Fitness",
    "Elapsed Time, average per simulation (s)",
 	              "Object Size, average per simulation (MB)",
 				  "Number of Clones, median",
 				  "Number of Iterations, median",
 				  "Final Time, median",
 				  "Total Population Size, median",
 				  "Total Population Size, mean",				  
 				  "Total Population Size, max.")], 
				  justify = c('left', 'left', rep('right', 8)),
				  ## caption = "\\label{tab:timingusual}Benchmarks under some common use cases, set 1." 
				  )	
	


## ----benchustable2, eval=TRUE, echo = FALSE------------------------------
## data(benchmark_3)

## knitr::kable(benchmark_3[, c("Model", "fitness", "time_per_simul",
##     "size_mb_per_simul", "NumClones.Median", "NumIter.Median",
## 	"FinalTime.Median", "TotalPopSize.Median", "TotalPopSize.Mean", 
## 	"TotalPopSize.Max.")], 
##     booktabs = TRUE,
## 	col.names = c("Model",
## 				  "Fitness", "Elapsed Time, average per simulation (s)",
## 	              "Object Size, average per simulation (MB)",
## 				  "Number of Clones, median",
## 				  "Number of Iterations, median",
## 				  "Final Time, median",
## 				  "Total Population Size, median",
## 				  "Total Population Size, mean",
## 				  "Total Population Size, max."
## 				  ),
## 	align = "c")

panderOptions('table.split.table', 99999999)
panderOptions('table.split.cells', 900)  ## For HTML
## panderOptions('table.split.cells', 8) ## For PDF


panderOptions('round', 3)
panderOptions('table.alignment.default', 'right')

pander(benchmark_3[ , c(
    "Model", "Fitness",
    "Elapsed Time, average per simulation (s)",
 	              "Object Size, average per simulation (MB)",
 				  "Number of Clones, median",
 				  "Number of Iterations, median",
 				  "Final Time, median",
 				  "Total Population Size, median",
 				  "Total Population Size, mean",				  
 				  "Total Population Size, max.")],
				  justify = c('left', 'left', rep('right', 8)),
				  ## caption = "\\label{tab:timingusual2}Benchmarks under some common use cases, set 2."
				  )	


## ----exp10000, echo = TRUE, eval = FALSE---------------------------------
## ng <- 10000
## u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2),
##                                       rep(-0.1, ng/2)))
## 
## t_e_10000 <- system.time(
##     e_10000 <- oncoSimulPop(5, u, model = "Exp", mu = 1e-7,
##                             detectionSize = 1e6,
##                             detectionDrivers = NA,
##                             detectionProb = NA,
##                             keepPhylog = TRUE,
##                             onlyCancer = FALSE,
##                             mutationPropGrowth = TRUE,
##                             mc.cores = 1))


## ----exp10000-out, echo = TRUE, eval = FALSE-----------------------------
## t_e_10000
## ##    user  system elapsed
## ##   4.368   0.196   4.566
## 
## summary(e_10000)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      5017      1180528       415116       143    7547
## ## 2      3726      1052061       603612       131    5746
## ## 3      4532      1100721       259510       132    6674
## ## 4      4150      1283115       829728        99    6646
## ## 5      4430      1139185       545958       146    6748
## 
## print(object.size(e_10000), units = "MB")
## ## 863.9 Mb
## 


## ----exp10000b, eval = FALSE, echo = TRUE--------------------------------
## t_e_10000b <- system.time(
##     e_10000b <- oncoSimulPop(5,
##                              u,
##                              model = "Exp",
##                              mu = 1e-7,
##                              detectionSize = 1e6,
##                              detectionDrivers = NA,
##                              detectionProb = NA,
##                              keepPhylog = TRUE,
##                              onlyCancer = FALSE,
##                              keepEvery = NA,
##                              mutationPropGrowth = TRUE,
##                              mc.cores = 1
##                              ))
## 


## ----exp10000b-out, echo = TRUE, eval = FALSE----------------------------
## t_e_10000b
## ##    user  system elapsed
## ##   5.484   0.100   5.585
## 
## summary(e_10000b)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      2465      1305094       727989        91    6447
## ## 2      2362      1070225       400329       204    8345
## ## 3      2530      1121164       436721       135    8697
## ## 4      2593      1206293       664494       125    8149
## ## 5      2655      1186994       327835       191    8572
## 
## print(object.size(e_10000b), units = "MB")
## ## 488.3 Mb
## 


## ----exp50000, echo = TRUE, eval = FALSE---------------------------------
## ng <- 50000
## u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2),
##                                       rep(-0.1, ng/2)))
## t_e_50000 <- system.time(
##     e_50000 <- oncoSimulPop(5,
##                             u,
##                             model = "Exp",
##                             mu = 1e-7,
##                             detectionSize = 1e6,
##                             detectionDrivers = NA,
##                             detectionProb = NA,
##                             keepPhylog = TRUE,
##                             onlyCancer = FALSE,
##                             keepEvery = NA,
##                             mutationPropGrowth = FALSE,
##                             mc.cores = 1
##                             ))
## 
## 
## t_e_50000
## ##    user  system elapsed
## ##  44.192   1.684  45.891
## 
## summary(e_50000)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      7367      1009949       335455     75.00   18214
## ## 2      8123      1302324       488469     63.65   17379
## ## 3      8408      1127261       270690     72.57   21144
## ## 4      8274      1138513       318152     80.59   20994
## ## 5      7520      1073131       690814     70.00   18569
## 
## print(object.size(e_50000), units = "MB")
## ## 7598.6 Mb


## ----exp50000np, echo = TRUE, eval = FALSE-------------------------------
## ng <- 50000
## u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2),
##                                       rep(-0.1, ng/2)))
## t_e_50000np <- system.time(
##     e_50000np <- oncoSimulPop(5,
##                               u,
##                               model = "Exp",
##                               mu = 1e-7,
##                               detectionSize = 1e6,
##                               detectionDrivers = NA,
##                               detectionProb = NA,
##                               keepPhylog = TRUE,
##                               onlyCancer = FALSE,
##                               keepEvery = 1,
##                               mutationPropGrowth = FALSE,
##                               mc.cores = 1
##                               ))
## 
## t_e_50000np
## ##   user  system elapsed
## ## 42.316   2.764  45.079
## 
## summary(e_50000np)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1     13406      1027949       410074     71.97   19469
## ## 2     12469      1071325       291852     66.00   17834
## ## 3     11821      1089834       245720     90.00   16711
## ## 4     14008      1165168       505607     77.61   19675
## ## 5     14759      1074621       205954     87.68   20597
## 
## print(object.size(e_50000np), units = "MB")
## ## 12748.4 Mb
## 


## ----exp50000mpg, echo = TRUE, eval = FALSE------------------------------
## 
## ng <- 50000
## u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2),
##                                       rep(-0.1, ng/2)))
## 
## t_e_50000c <- system.time(
##     e_50000c <- oncoSimulPop(5,
##                              u,
##                              model = "Exp",
##                              mu = 1e-7,
##                              detectionSize = 1e6,
##                              detectionDrivers = NA,
##                              detectionProb = NA,
##                              keepPhylog = TRUE,
##                              onlyCancer = FALSE,
##                              keepEvery = NA,
##                              mutationPropGrowth = TRUE,
##                              mc.cores = 1
##                              ))
## 
## t_e_50000c
## ##    user  system elapsed
## ## 84.228   2.416  86.665
## 
## summary(e_50000c)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1     11178      1241970       344479     84.74   27137
## ## 2     12820      1307086       203544     91.94   33448
## ## 3     10592      1126091       161057     83.81   26064
## ## 4     11883      1351114       148986     65.68   25396
## ## 5     10518      1101392       253523     99.79   26082
## 
## print(object.size(e_50000c), units = "MB")
## ## 10904.9 Mb
## 


## ----sizedetail, eval = FALSE, echo = TRUE-------------------------------
## r1 <- oncoSimulIndiv(u,
##                      model = "Exp",
##                      mu = 1e-7,
##                      detectionSize = 1e6,
##                      detectionDrivers = NA,
##                      detectionProb = NA,
##                      keepPhylog = TRUE,
##                      onlyCancer = FALSE,
##                      mutationPropGrowth = TRUE
##                      )
## summary(r1)[c(1, 8)]
## ##   NumClones  FinalTime
## ## 1      3887        345
## 
## print(object.size(r1), units = "MB")
## ## 160 Mb
## 
## ## Size of the two largest objects inside:
## sizes <- lapply(r1, function(x) object.size(x)/(1024^2))
## sort(unlist(sizes), decreasing = TRUE)[1:2]
## ## Genotypes pops.by.time
## ##       148.28        10.26
## 
## dim(r1$Genotypes)
## ## [1] 10000  3887


## ----mc50000_1, echo = TRUE, eval = FALSE--------------------------------
## ng <- 50000
## u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2),
##                                       rep(-0.1, ng/2)))
## 
## t_mc_50000_nmpg <- system.time(
##     mc_50000_nmpg <- oncoSimulPop(5,
##                                   u,
##                                   model = "McFL",
##                                   mu = 1e-7,
##                                   detectionSize = 1e6,
##                                   detectionDrivers = NA,
##                                   detectionProb = NA,
##                                   keepPhylog = TRUE,
##                                   onlyCancer = FALSE,
##                                   keepEvery = NA,
##                                   mutationPropGrowth = FALSE,
##                                   mc.cores = 1
##                                   ))
## t_mc_50000_nmpg
## ##   user  system elapsed
## ##  30.46    0.54   31.01
## 
## 
## summary(mc_50000_nmpg)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      1902      1002528       582752     284.2   31137
## ## 2      2159      1002679       404858     274.8   36905
## ## 3      2247      1002722       185678     334.5   42429
## ## 4      2038      1009606       493574     218.4   32519
## ## 5      2222      1004661       162628     291.0   38470
## 
## print(object.size(mc_50000_nmpg), units = "MB")
## ## 2057.6 Mb
## 


## ----mc50000_kp, echo = TRUE, eval = FALSE-------------------------------
## t_mc_50000_nmpg_k <- system.time(
##     mc_50000_nmpg_k <- oncoSimulPop(5,
##                                     u,
##                                     model = "McFL",
##                                     mu = 1e-7,
##                                     detectionSize = 1e6,
##                                     detectionDrivers = NA,
##                                     detectionProb = NA,
##                                     keepPhylog = TRUE,
##                                     onlyCancer = FALSE,
##                                     keepEvery = 1,
##                                     mutationPropGrowth = FALSE,
##                                     mc.cores = 1
##                                     ))
## 
## t_mc_50000_nmpg_k
## ##    user  system elapsed
## ##  30.000   1.712  31.714
## 
## summary(mc_50000_nmpg_k)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      8779      1000223       136453     306.7   38102
## ## 2      7442      1006563       428150     345.3   35139
## ## 3      8710      1003509       224543     252.3   35659
## ## 4      8554      1002537       103889     273.7   36783
## ## 5      8233      1003171       263005     301.8   35236
## 
## print(object.size(mc_50000_nmpg_k), units = "MB")
## ## 8101.4 Mb


## ----mc50000_popx, echo = TRUE, eval = FALSE-----------------------------
## ng <- 50000
## u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng/2),
##                                       rep(-0.1, ng/2)))
## 
## t_mc_50000_nmpg_3e6 <- system.time(
##     mc_50000_nmpg_3e6 <- oncoSimulPop(5,
##                                       u,
##                                       model = "McFL",
##                                       mu = 1e-7,
##                                       detectionSize = 3e6,
##                                       detectionDrivers = NA,
##                                       detectionProb = NA,
##                                       keepPhylog = TRUE,
##                                       onlyCancer = FALSE,
##                                       keepEvery = NA,
##                                       mutationPropGrowth = FALSE,
##                                       mc.cores = 1
##                                       ))
## t_mc_50000_nmpg_3e6
## ##    user  system elapsed
## ##  77.240   1.064  78.308
## 
## summary(mc_50000_nmpg_3e6)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      5487      3019083       836793     304.5   65121
## ## 2      4812      3011816       789146     286.3   53087
## ## 3      4463      3016896      1970957     236.6   45918
## ## 4      5045      3028142       956026     360.3   63464
## ## 5      4791      3029720       916692     358.1   55012
## 
## print(object.size(mc_50000_nmpg_3e6), units = "MB")
## ## 4759.3 Mb
## 												


## ----mc50000_mux, echo = TRUE, eval = FALSE------------------------------
## 
## t_mc_50000_nmpg_5mu <- system.time(
##     mc_50000_nmpg_5mu <- oncoSimulPop(5,
##                                       u,
##                                       model = "McFL",
##                                       mu = 5e-7,
##                                       detectionSize = 1e6,
##                                       detectionDrivers = NA,
##                                       detectionProb = NA,
##                                       keepPhylog = TRUE,
##                                       onlyCancer = FALSE,
##                                       keepEvery = NA,
##                                       mutationPropGrowth = FALSE,
##                                       mc.cores = 1
##                                       ))
## 
## t_mc_50000_nmpg_5mu
## ##    user  system elapsed
## ## 167.332   1.796 169.167
## 
## summary(mc_50000_nmpg_5mu)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      7963      1004415       408352     99.03   57548
## ## 2      8905      1010751       120155    130.30   74738
## ## 3      8194      1005465       274661     96.98   58546
## ## 4      9053      1014049       119943    112.23   75379
## ## 5      8982      1011817        95047     99.95   76757
## 
## print(object.size(mc_50000_nmpg_5mu), units = "MB")
## ## 8314.4 Mb


## ----mcf5muk, echo = TRUE, eval = FALSE----------------------------------
## t_mc_50000_nmpg_5mu_k <- system.time(
##     mc_50000_nmpg_5mu_k <- oncoSimulPop(5,
##                                         u,
##                                         model = "McFL",
##                                         mu = 5e-7,
##                                         detectionSize = 1e6,
##                                         detectionDrivers = NA,
##                                         detectionProb = NA,
##                                         keepPhylog = TRUE,
##                                         onlyCancer = FALSE,
##                                         keepEvery = 1,
##                                         mutationPropGrowth = FALSE,
##                                         mc.cores = 1
##                                         ))
## 												
## t_mc_50000_nmpg_5mu_k
## ##    user  system elapsed
## ## 174.404   5.068 179.481
## 
## summary(mc_50000_nmpg_5mu_k)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1     25294      1001597       102766     123.4   74524
## ## 2     23766      1006679       223010     124.3   71808
## ## 3     21755      1001379       203638     114.8   62609
## ## 4     24889      1012103       161003     119.3   75031
## ## 5     21844      1002927       255388     108.8   64556
## 
## print(object.size(mc_50000_nmpg_5mu_k), units = "MB")
## ## 22645.8 Mb
## 


## ----mc50000_2, echo = TRUE, eval = FALSE--------------------------------
## 
## t_mc_50000 <- system.time(
##     mc_50000 <- oncoSimulPop(5,
##                              u,
##                              model = "McFL",
##                              mu = 1e-7,
##                              detectionSize = 1e6,
##                              detectionDrivers = NA,
##                              detectionProb = NA,
##                              keepPhylog = TRUE,
##                              onlyCancer = FALSE,
##                              keepEvery = NA,
##                              mutationPropGrowth = TRUE,
##                              mc.cores = 1
##                              ))
## 
## t_mc_50000
## ##    user  system elapsed
## ## 303.352   2.808 306.223
## 
## summary(mc_50000)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1     13928      1010815       219814     210.9   91255
## ## 2     12243      1003267       214189     178.1   67673
## ## 3     13880      1014131       124354     161.4   88322
## ## 4     14104      1012941        75521     205.7   98583
## ## 5     12428      1005594       232603     167.4   70359
## 
## print(object.size(mc_50000), units = "MB")
## ## 12816.6 Mb
## 


## ----mcf5muk005, echo = TRUE, eval = FALSE-------------------------------
## t_mc_50000_nmpg_5mu_k <- system.time(
##     mc_50000_nmpg_5mu_k <- oncoSimulPop(2,
##                                         u,
##                                         model = "McFL",
##                                         mu = 5e-7,
##                                         detectionSize = 1e6,
##                                         detectionDrivers = NA,
##                                         detectionProb = NA,
##                                         keepPhylog = TRUE,
##                                         onlyCancer = FALSE,
##                                         keepEvery = 1,
##                                         mutationPropGrowth = FALSE,
##                                         mc.cores = 1
##                                         ))
## t_mc_50000_nmpg_5mu_k
## ##    user  system elapsed
## ## 305.512   5.164 310.711
## 
## summary(mc_50000_nmpg_5mu_k)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1     61737      1003273       104460  295.8731  204214
## ## 2     65072      1000540       133068  296.6243  210231
## 
## print(object.size(mc_50000_nmpg_5mu_k), units = "MB")
## ## 24663.6 Mb
## 


## ----mc50000_1_005, echo = TRUE, eval = FALSE----------------------------
## t_mc_50000_nmpg <- system.time(
##     mc_50000_nmpg <- oncoSimulPop(5,
##                                   u,
##                                   model = "McFL",
##                                   mu = 1e-7,
##                                   detectionSize = 1e6,
##                                   detectionDrivers = NA,
##                                   detectionProb = NA,
##                                   keepPhylog = TRUE,
##                                   onlyCancer = FALSE,
##                                   keepEvery = NA,
##                                   mutationPropGrowth = FALSE,
##                                   mc.cores = 1
##                                   ))
## t_mc_50000_nmpg
## ##    user  system elapsed
## ## 111.236   0.596 111.834
## 
## summary(mc_50000_nmpg)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      2646      1000700       217188   734.475  108566
## ## 2      2581      1001626       209873   806.500  107296
## ## 3      2903      1001409       125148   841.700  120859
## ## 4      2310      1000146       473948   906.300   91519
## ## 5      2704      1001290       448409   838.800  103556
## 
## print(object.size(mc_50000_nmpg), units = "MB")
## ## 2638.3 Mb
## 


## ----filog_exp50000_1, echo = TRUE, eval = FALSE-------------------------
## head(e_50000[[1]]$other$PhylogDF)
## ##   parent child   time
## ## 1         3679 0.8402
## ## 2         4754 1.1815
## ## 3        20617 1.4543
## ## 4        15482 2.3064
## ## 5         4431 3.7130
## ## 6        41915 4.0628
## 
## tail(e_50000[[1]]$other$PhylogDF)
## ##                          parent                            child time
## ## 20672               3679, 20282               3679, 20282, 22359 75.0
## ## 20673        3679, 17922, 22346        3679, 17922, 22346, 35811 75.0
## ## 20674                2142, 3679                2142, 3679, 25838 75.0
## ## 20675        3679, 17922, 19561        3679, 17922, 19561, 43777 75.0
## ## 20676 3679, 15928, 19190, 20282 3679, 15928, 19190, 20282, 49686 75.0
## ## 20677         2142, 3679, 16275         2142, 3679, 16275, 24201 75.0


## ----noplotlconephylog, echo = TRUE, eval = FALSE------------------------
## plotClonePhylog(e_50000[[1]]) ## plot not shown
## 


## ----ex-large-pop-size, eval = FALSE, echo = TRUE------------------------
## ng <- 50
## u <- allFitnessEffects(noIntGenes = c(rep(0.1, ng)))


## ----ex-large-mf, eval = FALSE, echo = TRUE------------------------------
## t_mc_k_50_1e11 <- system.time(
##     mc_k_50_1e11 <- oncoSimulPop(5,
##                                  u,
##                                  model = "McFL",
##                                  mu = 1e-7,
##                                  detectionSize = 1e11,
##                                  initSize = 1e5,
##                                  detectionDrivers = NA,
##                                  detectionProb = NA,
##                                  keepPhylog = TRUE,
##                                  onlyCancer = FALSE,
##                                  mutationPropGrowth = FALSE,
##                                  keepEvery = 1,
##                                  finalTime = 5000,
##                                  mc.cores = 1,
##                                  max.wall.time = 600
##                                  ))
## 
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## 
## t_mc_k_50_1e11
## ## user  system elapsed
## ## 613.612   0.040 613.664
## 
## summary(mc_k_50_1e11)[, c(1:3, 8, 9)]
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      5491 100328847809  44397848771  1019.950  942764
## ## 2      3194 100048090441  34834178374   789.675  888819
## ## 3      5745 100054219162  24412502660   927.950  929231
## ## 4      4017 101641197799  60932177160   750.725  480938
## ## 5      5393 100168156804  41659212367   846.250  898245
## 
## ## print(object.size(mc_k_50_1e11), units = "MB")
## ## 177.8 Mb
## 


## ----ex-large-exp, eval = FALSE, echo = TRUE-----------------------------
## t_exp_k_50_1e11 <- system.time(
##     exp_k_50_1e11 <- oncoSimulPop(5,
##                                   u,
##                                   model = "Exp",
##                                   mu = 1e-7,
##                                   detectionSize = 1e11,
##                                   initSize = 1e5,
##                                   detectionDrivers = NA,
##                                   detectionProb = NA,
##                                   keepPhylog = TRUE,
##                                   onlyCancer = FALSE,
##                                   mutationPropGrowth = FALSE,
##                                   keepEvery = 1,
##                                   finalTime = 5000,
##                                   mc.cores = 1,
##                                   max.wall.time = 600,
##                                   errorHitWallTime = FALSE,
##                                   errorHitMaxTries = FALSE
##                                   ))
## 
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Hitted wall time. Exiting.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Hitted wall time. Exiting.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Recoverable exception ti set to DBL_MIN. Rerunning.
## ## Hitted wall time. Exiting.
## ## Hitted wall time. Exiting.
## 
## t_exp_k_50_1e11
## ##     user   system  elapsed
## ## 2959.068    0.128 2959.556
## try(summary(exp_k_50_1e11)[, c(1:3, 8, 9)])
## ##   NumClones TotalPopSize LargestClone FinalTime NumIter
## ## 1      6078  65172752616  16529682757  235.7590 1883438
## ## 2      5370 106476643712  24662446729  232.0000 2516675
## ## 3      2711  21911284363  17945303353  224.8608  543698
## ## 4      2838  13241462284   2944300245  216.8091  372298
## ## 5      7289  76166784312  10941729810  240.0217 1999489
## 
## print(object.size(exp_k_50_1e11), units = "MB")
## ## 53.5 Mb
## 


## ------------------------------------------------------------------------
m4 <- data.frame(G = c("WT", "A", "B", "A, B"), F = c(1, 2, 3, 4))


## ------------------------------------------------------------------------
fem4 <- allFitnessEffects(genotFitness = m4)


## ------------------------------------------------------------------------
try(plot(fem4))


## ---- fig.width=6.5, fig.height = 6.5------------------------------------
plotFitnessLandscape(evalAllGenotypes(fem4))


## ------------------------------------------------------------------------
evalAllGenotypes(fem4, addwt = TRUE)


## ------------------------------------------------------------------------
plotFitnessLandscape(evalAllGenotypes(fem4))


## ------------------------------------------------------------------------
m6 <- cbind(c(1, 1), c(1, 0), c(2, 3))
fem6 <- allFitnessEffects(genotFitness = m6)
evalAllGenotypes(fem6, addwt = TRUE)
## plot(fem6)


## ---- fig.width=6.5, fig.height = 6.5------------------------------------
plotFitnessLandscape(evalAllGenotypes(fem6))


## ------------------------------------------------------------------------
r <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                 Fitness = c("10*f_", 
                             "10*f_1", 
                             "50*f_2", 
                             "200*(f_1 + f_2) + 50*f_1_2"))

afe <- allFitnessEffects(genotFitness = r, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(2500, 2000, 5500, 700))

plotFitnessLandscape(evalAllGenotypes(afe))


## ----mcflparam-----------------------------------------------------------
sp <- 1e-3 
spp <- -sp/(1 + sp)


## ------------------------------------------------------------------------

ai1 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c(0.05, -.2, .1)), order = FALSE)


## ------------------------------------------------------------------------
ai1


## ------------------------------------------------------------------------
all(ai1[, "Fitness"]  == c( (1 + .05), (1 - .2), (1 + .1),
       (1 + .05) * (1 - .2),
       (1 + .05) * (1 + .1),
       (1 - .2) * (1 + .1),
       (1 + .05) * (1 - .2) * (1 + .1)))



## ------------------------------------------------------------------------
(ai2 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c(0.05, -.2, .1)), order = TRUE,
    addwt = TRUE))



## ---- fig.height=4-------------------------------------------------------
data(examplesFitnessEffects)
plot(examplesFitnessEffects[["cbn1"]])


## ------------------------------------------------------------------------

cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = -0.9,
                 typeDep = "MN")

cbn1 <- allFitnessEffects(cs)



## ---- fig.height=3-------------------------------------------------------
plot(cbn1)


## ---- fig.height=5-------------------------------------------------------
plot(cbn1, "igraph")


## ---- fig.height=5-------------------------------------------------------
library(igraph) ## to make the reingold.tilford layout available
plot(cbn1, "igraph", layout = layout.reingold.tilford)


## ------------------------------------------------------------------------
gfs <- evalAllGenotypes(cbn1, order = FALSE, addwt = TRUE)

gfs[1:15, ]


## ------------------------------------------------------------------------

c1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.1, -.2), c(-.05, -.06, -.07)),
                 typeDep = "MN")

try(fc1 <- allFitnessEffects(c1))



## ------------------------------------------------------------------------
c1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "MN")

cbn2 <- allFitnessEffects(c1)



## ------------------------------------------------------------------------
gcbn2 <- evalAllGenotypes(cbn2, order = FALSE)
gcbn2[1:10, ]


## ------------------------------------------------------------------------
gcbn2o <- evalAllGenotypes(cbn2, order = TRUE, max = 1956)
gcbn2o[1:10, ]


## ------------------------------------------------------------------------
all.equal(
        gcbn2[c(1:21, 22, 28, 41, 44, 56, 63 ) , "Fitness"],
        c(1.01, 1.02, 0.1, 1.03, 1.04, 0.05,
          1.01 * c(1.02, 0.1, 1.03, 1.04, 0.05),
          1.02 * c(0.10, 1.03, 1.04, 0.05),
          0.1 * c(1.03, 1.04, 0.05),
          1.03 * c(1.04, 0.05),
          1.04 * 0.05,
          1.01 * 1.02 * 1.1,
          1.01 * 0.1 * 0.05,
          1.03 * 1.04 * 0.05,
          1.01 * 1.02 * 1.1 * 0.05,
          1.03 * 1.04 * 1.2 * 0.1, ## notice this
          1.01 * 1.02 * 1.03 * 1.04 * 1.1 * 1.2
          ))


## ------------------------------------------------------------------------
gcbn2[56, ] 
all.equal(gcbn2[56, "Fitness"], 1.03 * 1.04 * 1.2 * 0.10)


## ------------------------------------------------------------------------

s1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "SM")

smn1 <- allFitnessEffects(s1)



## ---- fig.height=3-------------------------------------------------------
plot(smn1)


## ------------------------------------------------------------------------
gsmn1 <- evalAllGenotypes(smn1, order = FALSE)



## ------------------------------------------------------------------------
gcbn2[c(8, 12, 22), ]
gsmn1[c(8, 12, 22), ]

gcbn2[c(20:21, 28), ]
gsmn1[c(20:21, 28), ]


## ------------------------------------------------------------------------

x1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "XMPN")

xor1 <- allFitnessEffects(x1)



## ---- fig.height=3-------------------------------------------------------
plot(xor1)


## ------------------------------------------------------------------------

gxor1 <- evalAllGenotypes(xor1, order = FALSE)



## ------------------------------------------------------------------------
gxor1[c(22, 41), ] 
c(1.01 * 1.02 * 0.1, 1.03 * 1.04 * 0.05)


## ------------------------------------------------------------------------
gxor1[28, ] 
1.01 * 1.1 * 1.2


## ------------------------------------------------------------------------
gxor1[44, ] 
1.01 * 1.02 * 0.1 * 1.2


## ------------------------------------------------------------------------

p3 <- data.frame(
    parent = c(rep("Root", 4), "a", "b", "d", "e", "c", "f"),
    child = c("a", "b", "d", "e", "c", "c", "f", "f", "g", "g"),
    s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
    sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
    typeDep = c(rep("--", 4), 
                "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
fp3 <- allFitnessEffects(p3)


## ---- fig.height=3-------------------------------------------------------
plot(fp3)


## ---- fig.height=6-------------------------------------------------------
plot(fp3, "igraph", layout.reingold.tilford)


## ------------------------------------------------------------------------

gfp3 <- evalAllGenotypes(fp3, order = FALSE)



## ------------------------------------------------------------------------
gfp3[c(9, 24, 29, 59, 60, 66, 119, 120, 126, 127), ]

c(1.01 * 1.1, 1.03 * .05, 1.01 * 1.02 * 0.1, 0.1 * 0.05 * 1.3,
  1.03 * 1.04 * 1.2, 1.01 * 1.02 * 0.1 * 0.05,
  0.1 * 1.03 * 1.04 * 1.2 * 1.3,
  1.01 * 1.02 * 0.1 * 1.03 * 1.04 * 1.2,
  1.02 * 1.1 * 1.03 * 1.04 * 1.2 * 1.3,
  1.01 * 1.02 * 1.03 * 1.04 * 0.1 * 1.2 * 1.3)



## ------------------------------------------------------------------------
s <- 0.2
sboth <- (1/(1 + s)) - 1
m0 <- allFitnessEffects(data.frame(
    parent = c("Root", "Root", "a1", "a2"),
    child = c("a1", "a2", "b", "b"),
    s = s,
    sh = -1,
    typeDep = "OR"),
                        epistasis = c("a1:a2" = sboth))
evalAllGenotypes(m0, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
s <- 0.2
m1 <- allFitnessEffects(data.frame(
    parent = c("Root", "A"),
    child = c("A", "B"),
    s = s,
    sh = -1,
    typeDep = "OR"),
                        geneToModule = c("Root" = "Root",
                                         "A" = "a1, a2",
                                         "B" = "b1"))
evalAllGenotypes(m1, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
fnme <- allFitnessEffects(epistasis = c("A" = 0.1,
                                        "B" = 0.2),
                          geneToModule = c("A" = "a1, a2",
                                           "B" = "b1"))
evalAllGenotypes(fnme, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
fnme2 <- allFitnessEffects(epistasis = c("A" = 0.1,
                                        "B" = 0.2),
                          geneToModule = c(
                              "Root" = "Root",
                              "A" = "a1, a2",
                              "B" = "b1"))
evalAllGenotypes(fnme, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
p4 <- data.frame(
    parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
    child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
    s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
    sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
    typeDep = c(rep("--", 4), 
                "XMPN", "XMPN", "MN", "MN", "SM", "SM"))

fp4m <- allFitnessEffects(
    p4,
    geneToModule = c("Root" = "Root", "A" = "a1",
                     "B" = "b1, b2", "C" = "c1",
                     "D" = "d1, d2", "E" = "e1",
                     "F" = "f1, f2", "G" = "g1"))


## ---- fig.height=3-------------------------------------------------------
plot(fp4m)


## ---- fig.height=3-------------------------------------------------------
plot(fp4m, expandModules = TRUE)


## ---- fig.height=6-------------------------------------------------------
plot(fp4m, "igraph", layout = layout.reingold.tilford, 
     expandModules = TRUE)



## ------------------------------------------------------------------------
gfp4 <- evalAllGenotypes(fp4m, order = FALSE, max = 1024)


## ------------------------------------------------------------------------
gfp4[c(12, 20, 21, 40, 41, 46, 50, 55, 64, 92,
       155, 157, 163, 372, 632, 828), ]

c(1.01 * 1.02, 1.02, 1.02 * 1.1, 0.1 * 1.3, 1.03, 
  1.03 * 1.04, 1.04 * 0.05, 0.05 * 1.3,  
  1.01 * 1.02 * 0.1, 1.02 * 1.1, 0.1 * 0.05 * 1.3,
  1.03 * 0.05, 1.03 * 0.05, 1.03 * 1.04 * 1.2, 1.03 * 1.04 * 1.2, 
  1.02 * 1.1 * 1.03 * 1.04 * 1.2 * 1.3)



## ------------------------------------------------------------------------
o3 <- allFitnessEffects(orderEffects = c(
                            "F > D > M" = -0.3,
                            "D > F > M" = 0.4,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.5),
                        geneToModule =
                            c("M" = "m",
                              "F" = "f",
                              "D" = "d") )


(ag <- evalAllGenotypes(o3, addwt = TRUE, order = TRUE))


## ------------------------------------------------------------------------

ofe1 <- allFitnessEffects(
    orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
    geneToModule =
        c("F" = "f1, f2",
          "D" = "d1, d2") )

ag <- evalAllGenotypes(ofe1, order = TRUE)



## ------------------------------------------------------------------------
ag[5:16,]


## ------------------------------------------------------------------------
ag[c(17, 39, 19, 29), ]


## ------------------------------------------------------------------------
ag[c(17, 39, 19, 29), "Fitness"] == c(1.4, 0.7, 1.4, 0.7)


## ------------------------------------------------------------------------
ag[c(43, 44),]
ag[c(43, 44), "Fitness"] == c(1.4, 1.4)


## ------------------------------------------------------------------------
all(ag[41:52, "Fitness"] == 1.4)


## ------------------------------------------------------------------------
all(ag[53:64, "Fitness"] == 0.7)


## ------------------------------------------------------------------------

ofe2 <- allFitnessEffects(
    orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
    geneToModule =
        c("F" = "f1, f2, f3",
          "D" = "d1, d2") )
ag2 <- evalAllGenotypes(ofe2, max = 325, order = TRUE)



## ------------------------------------------------------------------------
all(ag2[grep("^d.*f.*", ag2[, 1]), "Fitness"] == 1.4)
all(ag2[grep("^f.*d.*", ag2[, 1]), "Fitness"] == 0.7)
oe <- c(grep("^f.*d.*", ag2[, 1]), grep("^d.*f.*", ag2[, 1]))
all(ag2[-oe, "Fitness"] == 1)


## ------------------------------------------------------------------------

foi1 <- allFitnessEffects(
    orderEffects = c("D>B" = -0.2, "B > D" = 0.3),
    noIntGenes = c("A" = 0.05, "C" = -.2, "E" = .1))



## ------------------------------------------------------------------------
foi1[c("geneModule", "long.geneNoInt")]


## ------------------------------------------------------------------------
agoi1 <- evalAllGenotypes(foi1,  max = 325, order = TRUE)
head(agoi1)


## ------------------------------------------------------------------------
rn <- 1:nrow(agoi1)
names(rn) <- agoi1[, 1]

agoi1[rn[LETTERS[1:5]], "Fitness"] == c(1.05, 1, 0.8, 1, 1.1)



## ------------------------------------------------------------------------
agoi1[grep("^A > [BD]$", names(rn)), "Fitness"] == 1.05
agoi1[grep("^C > [BD]$", names(rn)), "Fitness"] == 0.8
agoi1[grep("^E > [BD]$", names(rn)), "Fitness"] == 1.1
agoi1[grep("^[BD] > A$", names(rn)), "Fitness"] == 1.05
agoi1[grep("^[BD] > C$", names(rn)), "Fitness"] == 0.8
agoi1[grep("^[BD] > E$", names(rn)), "Fitness"] == 1.1


## ------------------------------------------------------------------------
all.equal(agoi1[230:253, "Fitness"] ,
          rep((1 - 0.2) * 1.05 * 0.8 * 1.1, 24))


## ------------------------------------------------------------------------
sa <- 0.2
sb <- 0.3
sab <- 0.7

e2 <- allFitnessEffects(epistasis =
                            c("A: -B" = sa,
                              "-A:B" = sb,
                              "A : B" = sab))
evalAllGenotypes(e2, order = FALSE, addwt = TRUE)



## ------------------------------------------------------------------------
s2 <- ((1 + sab)/((1 + sa) * (1 + sb))) - 1

e3 <- allFitnessEffects(epistasis =
                            c("A" = sa,
                              "B" = sb,
                              "A : B" = s2))
evalAllGenotypes(e3, order = FALSE, addwt = TRUE)



## ------------------------------------------------------------------------
sa <- 0.1
sb <- 0.15
sc <- 0.2
sab <- 0.3
sbc <- -0.25
sabc <- 0.4

sac <- (1 + sa) * (1 + sc) - 1

E3A <- allFitnessEffects(epistasis =
                            c("A:-B:-C" = sa,
                              "-A:B:-C" = sb,
                              "-A:-B:C" = sc,
                              "A:B:-C" = sab,
                              "-A:B:C" = sbc,
                              "A:-B:C" = sac,
                              "A : B : C" = sabc)
                                                )

evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)



## ------------------------------------------------------------------------

sa <- 0.1
sb <- 0.15
sc <- 0.2
sab <- 0.3
Sab <- ( (1 + sab)/((1 + sa) * (1 + sb))) - 1
Sbc <- ( (1 + sbc)/((1 + sb) * (1 + sc))) - 1
Sabc <- ( (1 + sabc)/
          ( (1 + sa) * (1 + sb) * (1 + sc) *
            (1 + Sab) * (1 + Sbc) ) ) - 1

E3B <- allFitnessEffects(epistasis =
                             c("A" = sa,
                               "B" = sb,
                               "C" = sc,
                               "A:B" = Sab,
                               "B:C" = Sbc,
                               ## "A:C" = sac, ## not needed now
                               "A : B : C" = Sabc)
                                                )
evalAllGenotypes(E3B, order = FALSE, addwt = FALSE)



## ------------------------------------------------------------------------
all(evalAllGenotypes(E3A, order = FALSE, addwt = FALSE) == 
    evalAllGenotypes(E3B, order = FALSE, addwt = FALSE))


## ------------------------------------------------------------------------

sa <- 0.2
sb <- 0.3
sab <- 0.7

em <- allFitnessEffects(epistasis =
                            c("A: -B" = sa,
                              "-A:B" = sb,
                              "A : B" = sab),
                        geneToModule = c("A" = "a1, a2",
                                         "B" = "b1, b2"))
evalAllGenotypes(em, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
s2 <- ((1 + sab)/((1 + sa) * (1 + sb))) - 1

em2 <- allFitnessEffects(epistasis =
                            c("A" = sa,
                              "B" = sb,
                              "A : B" = s2),
                         geneToModule = c("A" = "a1, a2",
                                         "B" = "b1, b2")
                         )
evalAllGenotypes(em2, order = FALSE, addwt = TRUE)



## ------------------------------------------------------------------------

fnme <- allFitnessEffects(epistasis = c("A" = 0.1,
                                        "B" = 0.2),
                          geneToModule = c("A" = "a1, a2",
                                           "B" = "b1, b2, b3"))

evalAllGenotypes(fnme, order = FALSE, addwt = TRUE)



## ------------------------------------------------------------------------
fnme <- allFitnessEffects(epistasis = c("A" = 0.1,
                                        "B" = 0.2,
                                        "A : B" = 0.0),
                          geneToModule = c("A" = "a1, a2",
                                           "B" = "b1, b2, b3"))

evalAllGenotypes(fnme, order = FALSE, addwt = TRUE)



## ------------------------------------------------------------------------
s <- 0.2
sv <- allFitnessEffects(epistasis = c("-A : B" = -1,
                                      "A : -B" = -1,
                                      "A:B" = s))


## ------------------------------------------------------------------------
(asv <- evalAllGenotypes(sv, order = FALSE, addwt = TRUE))


## ------------------------------------------------------------------------
evalAllGenotypes(sv, order = TRUE, addwt = TRUE)


## ------------------------------------------------------------------------
sa <- -0.1
sb <- -0.2
sab <- 0.25
sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                             "A : -B" = sa,
                             "A:B" = sab),
                         geneToModule = c(
                             "A" = "a1, a2",
                             "B" = "b"))
evalAllGenotypes(sv2, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
evalAllGenotypes(sv2, order = TRUE, addwt = TRUE)


## ------------------------------------------------------------------------
sa <- 0.1
sb <- 0.2
sab <- -0.8
sm1 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                             "A : -B" = sa,
                             "A:B" = sab))
evalAllGenotypes(sm1, order = FALSE, addwt = TRUE)



## ------------------------------------------------------------------------
evalAllGenotypes(sm1, order = TRUE, addwt = TRUE)


## ------------------------------------------------------------------------
evalAllGenotypes(sv, order = FALSE, addwt = TRUE, model = "Bozic")


## ------------------------------------------------------------------------
s <- 0.2
svB <- allFitnessEffects(epistasis = c("-A : B" = -Inf,
                                      "A : -B" = -Inf,
                                      "A:B" = s))
evalAllGenotypes(svB, order = FALSE, addwt = TRUE, model = "Bozic")


## ------------------------------------------------------------------------

s <- 1
svB1 <- allFitnessEffects(epistasis = c("-A : B" = -Inf,
                                       "A : -B" = -Inf,
                                       "A:B" = s))

evalAllGenotypes(svB1, order = FALSE, addwt = TRUE, model = "Bozic")


s <- 3
svB3 <- allFitnessEffects(epistasis = c("-A : B" = -Inf,
                                       "A : -B" = -Inf,
                                       "A:B" = s))

evalAllGenotypes(svB3, order = FALSE, addwt = TRUE, model = "Bozic")




## ------------------------------------------------------------------------
i1 <- allFitnessEffects(noIntGenes = c(1, 0.5))
evalAllGenotypes(i1, order = FALSE, addwt = TRUE, 
                 model = "Bozic")
				 
i1_b <- oncoSimulIndiv(i1, model = "Bozic")



## ------------------------------------------------------------------------
evalAllGenotypes(i1, order = FALSE, addwt = TRUE, 
                 model = "Exp")
i1_e <- oncoSimulIndiv(i1, model = "Exp")
summary(i1_e)


## ------------------------------------------------------------------------
p4 <- data.frame(
    parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
    child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
    s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
    sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
    typeDep = c(rep("--", 4), 
                "XMPN", "XMPN", "MN", "MN", "SM", "SM"))

oe <- c("C > F" = -0.1, "H > I" = 0.12)
sm <- c("I:J"  = -1)
sv <- c("-K:M" = -.5, "K:-M" = -.5)
epist <- c(sm, sv)

modules <- c("Root" = "Root", "A" = "a1",
             "B" = "b1, b2", "C" = "c1",
             "D" = "d1, d2", "E" = "e1",
             "F" = "f1, f2", "G" = "g1",
             "H" = "h1, h2", "I" = "i1",
             "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")

set.seed(1) ## for reproducibility
noint <- rexp(5, 10)
names(noint) <- paste0("n", 1:5)

fea <- allFitnessEffects(rT = p4, epistasis = epist,
                         orderEffects = oe,
                         noIntGenes = noint,
                         geneToModule = modules)


## ---- fig.height=5.5-----------------------------------------------------
plot(fea)


## ---- fig.height=5.5-----------------------------------------------------
plot(fea, "igraph")


## ---- fig.height=5.5-----------------------------------------------------
plot(fea, expandModules = TRUE)


## ---- fig.height=6.5-----------------------------------------------------
plot(fea, "igraph", expandModules = TRUE)


## ------------------------------------------------------------------------

evalGenotype("k1 > i1 > h2", fea) ## 0.5
evalGenotype("k1 > h1 > i1", fea) ## 0.5 * 1.12

evalGenotype("k2 > m1 > h1 > i1", fea) ## 1.12

evalGenotype("k2 > m1 > h1 > i1 > c1 > n3 > f2", fea) 
## 1.12 * 0.1 * (1 + noint[3]) * 0.05 * 0.9



## ------------------------------------------------------------------------

randomGenotype <- function(fe, ns = NULL) {
    gn <- setdiff(c(fe$geneModule$Gene,
                    fe$long.geneNoInt$Gene), "Root")
    if(is.null(ns)) ns <- sample(length(gn), 1)
    return(paste(sample(gn, ns), collapse = " > "))
}

set.seed(2) ## for reproducibility

evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  k2 > i1 > c1 > n1 > m1
##  Individual s terms are : 0.0755182 -0.9
##  Fitness:  0.107552 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  n2 > h1 > h2
##  Individual s terms are : 0.118164
##  Fitness:  1.11816 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  d2 > k2 > c1 > f2 > n4 > m1 > n3 > f1 > b1 > g1 > n5 > h1 > j2
##  Individual s terms are : 0.0145707 0.0139795 0.0436069 0.02 0.1 0.03 -0.95 0.3 -0.1
##  Fitness:  0.0725829 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  h2 > c1 > f1 > n2 > b2 > a1 > n1 > i1
##  Individual s terms are : 0.0755182 0.118164 0.01 0.02 -0.9 -0.95 -0.1 0.12
##  Fitness:  0.00624418 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  h2 > j1 > m1 > d2 > i1 > b2 > k2 > d1 > b1 > n3 > n1 > g1 > h1 > c1 > k1 > e1 > a1 > f1 > n5 > f2
##  Individual s terms are : 0.0755182 0.0145707 0.0436069 0.01 0.02 -0.9 0.03 0.04 0.2 0.3 -1 -0.1 0.12
##  Fitness:  0 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  n1 > m1 > n3 > i1 > j1 > n5 > k1
##  Individual s terms are : 0.0755182 0.0145707 0.0436069 -1
##  Fitness:  0 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  d2 > n1 > g1 > f1 > f2 > c1 > b1 > d1 > k1 > a1 > b2 > i1 > n4 > h2 > n2
##  Individual s terms are : 0.0755182 0.118164 0.0139795 0.01 0.02 -0.9 0.03 -0.95 0.3 -0.5
##  Fitness:  0.00420528 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  j1 > f1 > j2 > a1 > n4 > c1 > n3 > k1 > d1 > h1
##  Individual s terms are : 0.0145707 0.0139795 0.01 0.1 0.03 -0.95 -0.5
##  Fitness:  0.0294308 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  n5 > f2 > f1 > h2 > n4 > c1 > n3 > b1
##  Individual s terms are : 0.0145707 0.0139795 0.0436069 0.02 0.1 -0.95
##  Fitness:  0.0602298 
evalGenotype(randomGenotype(fea), fea, echo = TRUE, verbose = TRUE)
## Genotype:  h1 > d1 > f2
##  Individual s terms are : 0.03 -0.95
##  Fitness:  0.0515 




## ------------------------------------------------------------------------

muvar2 <- c("U" = 1e-6, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3,
            "D" = 1e-4)
ni1 <- rep(0, 5)
names(ni1) <- names(muvar2) ## We use the same names, of course
fe1 <- allFitnessEffects(noIntGenes = ni1)
bb <- oncoSimulIndiv(fe1, 
                     mu = muvar2, onlyCancer = FALSE,
                     initSize = 1e5,
                     finalTime = 25,
                     seed =NULL)



## ------------------------------------------------------------------------
fe2 <- allFitnessEffects(noIntGenes =
                         c(a1 = 0.1, a2 = 0.2,
                           b1 = 0.01, b2 = 0.3, b3 = 0.2,
                           c1 = 0.3, c2 = -0.2))

fm2 <- allMutatorEffects(epistasis = c("A" = 5,
                                       "B" = 10,
                                       "C" = 3),
                         geneToModule = c("A" = "a1, a2",
                                          "B" = "b1, b2, b3",
                                          "C" = "c1, c2"))

## Show the fitness effect of a specific genotype
evalGenotype("a1, c2", fe2, verbose = TRUE)

## Show the mutator effect of a specific genotype
evalGenotypeMut("a1, c2", fm2, verbose = TRUE)

## Fitness and mutator of a specific genotype
evalGenotypeFitAndMut("a1, c2", fe2, fm2, verbose = TRUE)


## ---- eval=FALSE---------------------------------------------------------
## ## Show only all the fitness effects
## evalAllGenotypes(fe2, order = FALSE)
## 
## ## Show only all mutator effects
## evalAllGenotypesMut(fm2)
## 
## ## Show all fitness and mutator
## evalAllGenotypesFitAndMut(fe2, fm2, order = FALSE)


## ------------------------------------------------------------------------

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


## ------------------------------------------------------------------------
## These only affect mutation, not fitness
evalGenotypeFitAndMut("a1, a2", fe3, fm3, verbose = TRUE)
evalGenotypeFitAndMut("a1, b3", fe3, fm3, verbose = TRUE)

## These only affect fitness: the mutator multiplier is 1
evalGenotypeFitAndMut("g1", fe3, fm3, verbose = TRUE)                      
evalGenotypeFitAndMut("g3, g9", fe3, fm3, verbose = TRUE)

## These affect both
evalGenotypeFitAndMut("g3, g9, a2, b3", fe3, fm3, verbose = TRUE)


## ------------------------------------------------------------------------
set.seed(1) ## so that it is easy to reproduce
mue1 <- oncoSimulIndiv(fe3, muEF = fm3, 
                       mu = 1e-6,
                       initSize = 1e5,
                       model = "McFL",
                       detectionSize = 5e6,
                       finalTime = 500,
                       onlyCancer = FALSE)


## ---- eval=FALSE---------------------------------------------------------
## ## We do not show this in the vignette to avoid cluttering it
## ## with output
## mue1


## ---- eval=FALSE---------------------------------------------------------
## 
## set.seed(1) ## for reproducibility
## ## 17 genes, 7 with no direct fitness effects
## ni <- c(rep(0, 7), runif(10, min = -0.01, max = 0.1))
## names(ni) <- c("a1", "a2", "b1", "b2", "b3", "c1", "c2",
##                paste0("g", 1:10))
## 
## ## Next is for nicer figure labeling.
## ## Consider as drivers genes with s >0
## gp <- which(ni > 0)
## 
## fe3 <- allFitnessEffects(noIntGenes = ni,
##                          drvNames = names(ni)[gp])
## 
## 
## set.seed(12)
## mue1 <- oncoSimulIndiv(fe3, muEF = fm3,
##                        mu = 1e-6,
##                        initSize = 1e5,
##                        model = "McFL",
##                        detectionSize = 5e6,
##                        finalTime = 270,
##                        keepPhylog = TRUE,
##                        onlyCancer = FALSE)
## mue1
## ## If you decrease N even further it gets even more cluttered
## op <- par(ask = TRUE)
## plotClonePhylog(mue1, N = 10, timeEvents = TRUE)
## plot(mue1, plotDrivers = TRUE, addtot = TRUE,
##      plotDiversity = TRUE)
## 	
## ## The stacked plot is slow; be patient
## ## Most clones have tiny population sizes, and their lines
## ## are piled on top of each other
## plot(mue1, addtot = TRUE,
##      plotDiversity = TRUE, type = "stacked")
## par(op)


## ------------------------------------------------------------------------

d1 <- -0.05 ## single mutant fitness 0.95
d2 <- -0.08 ## double mutant fitness 0.92
d3 <- 0.2   ## triple mutant fitness 1.2
s2 <- ((1 + d2)/(1 + d1)^2) - 1
s3 <- ( (1 + d3)/((1 + d1)^3 * (1 + s2)^3) ) - 1

wb <- allFitnessEffects(
    epistasis = c(
        "A" = d1,
        "B" = d1,
        "C" = d1,
        "A:B" = s2,
        "A:C" = s2,
        "B:C" = s2,
        "A:B:C" = s3))


## ---- fig.width=6.5, fig.height=5----------------------------------------
plotFitnessLandscape(wb, use_ggrepel = TRUE) 


## ---- fig.width=6.5, fig.height=5----------------------------------------
(ewb <- evalAllGenotypes(wb, order = FALSE))
plot(ewb, use_ggrepel = TRUE) 



## ----wasthis111, fig.width=9.5, fig.height=9.5---------------------------
par(cex = 0.7)
pancr <- allFitnessEffects(
    data.frame(parent = c("Root", rep("KRAS", 4), 
                   "SMAD4", "CDNK2A", 
                   "TP53", "TP53", "MLL3"),
               child = c("KRAS","SMAD4", "CDNK2A", 
                   "TP53", "MLL3",
                   rep("PXDN", 3), rep("TGFBR2", 2)),
               s = 0.1,
               sh = -0.9,
               typeDep = "MN"))
plot(evalAllGenotypes(pancr, order = FALSE), use_ggrepel = TRUE)



## ------------------------------------------------------------------------
K <- 4
sp <- 1e-5
sdp <- 0.015
sdplus <- 0.05
sdminus <- 0.1
cnt <- (1 + sdplus)/(1 + sdminus)
prod_cnt <- cnt - 1
bauer <- data.frame(parent = c("Root", rep("D", K)),
                    child = c("D", paste0("s", 1:K)),
                    s = c(prod_cnt, rep(sdp, K)),
                    sh = c(0, rep(sp, K)),
                    typeDep = "MN")
fbauer <- allFitnessEffects(bauer)
(b1 <- evalAllGenotypes(fbauer, order = FALSE, addwt = TRUE))



## ---- fig.height=3-------------------------------------------------------
plot(fbauer)


## ---- fig.width=6, fig.height=6------------------------------------------
plot(b1, use_ggrepel = TRUE)


## ------------------------------------------------------------------------
m1 <- expand.grid(D = c(1, 0), s1 = c(1, 0), s2 = c(1, 0),
                  s3 = c(1, 0), s4 = c(1, 0))

fitness_bauer <- function(D, s1, s2, s3, s4, 
                          sp = 1e-5, sdp = 0.015, sdplus = 0.05,
                          sdminus = 0.1) {
    if(!D) {
        b <- 0.5 * ( (1 + sp)^(sum(c(s1, s2, s3, s4))))
    } else {
        b <- 0.5 * 
            (((1 + sdplus)/(1 + sdminus)  *
              (1 + sdp)^(sum(c(s1, s2, s3, s4)))))
    }
    fitness <- b - (1 - b)
    our_fitness <- 1 + fitness ## prevent negative fitness and
    ## make wt fitness = 1
    return(our_fitness)
}

m1$Fitness <- 
    apply(m1, 1, function(x) do.call(fitness_bauer, as.list(x)))

bauer2 <- allFitnessEffects(genotFitness = m1)


## ------------------------------------------------------------------------
evalAllGenotypes(bauer2, order = FALSE, addwt = TRUE)


## ---- echo=FALSE, fig.height=4, fig.width=4------------------------------

df1 <- data.frame(x = c(1, 1.2, 1.4), f = c(1, 1.2, 1.2),
                 names = c("wt", "A", "B"))
plot(df1[, 2] ~ df1[, 1], axes = TRUE, xlab= "", 
     ylab = "Fitness", xaxt = "n", yaxt = "n", ylim = c(1, 1.21))
segments(1, 1, 1.2, 1.2)
segments(1, 1, 1.4, 1.2)
text(1, 1, "wt", pos = 4)
text(1.2, 1.2, "A", pos = 2)
text(1.4, 1.2, "B", pos = 2)
## axis(1,  tick = FALSE, labels = FALSE)
## axis(2,  tick = FALSE, labels = FALSE)


## ------------------------------------------------------------------------
s <- 0.1 ## or whatever number
m1a1 <- allFitnessEffects(data.frame(parent = c("Root", "Root"),
                                     child = c("A", "B"),
                                     s = s,
                                     sh = 0,
                                     typeDep = "MN"))
evalAllGenotypes(m1a1, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
s <- 0.1
sab <- 0.3
m1a2 <- allFitnessEffects(epistasis = c("A:-B" = s,
                                        "-A:B" = s,
                                        "A:B" = sab))
evalAllGenotypes(m1a2, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
sab3 <- ((1 + sab)/((1 + s)^2)) - 1
m1a3 <- allFitnessEffects(data.frame(parent = c("Root", "Root"),
                                     child = c("A", "B"),
                                     s = s,
                                     sh = 0,
                                     typeDep = "MN"),
                          epistasis = c("A:B" = sab3))
evalAllGenotypes(m1a3, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
all.equal(evalAllGenotypes(m1a2, order = FALSE, addwt = TRUE),
          evalAllGenotypes(m1a3, order = FALSE, addwt = TRUE))


## ---- echo=FALSE, fig.width=4, fig.height=4------------------------------

df1 <- data.frame(x = c(1, 1.2, 1.2, 1.4), f = c(1, 0.4, 0.3, 1.3),
                 names = c("wt", "A", "B", "AB"))
plot(df1[, 2] ~ df1[, 1], axes = TRUE, xlab= "", ylab = "Fitness",
     xaxt = "n", yaxt = "n", ylim = c(0.29, 1.32))
segments(1, 1, 1.2, 0.4)
segments(1, 1, 1.2, 0.3)
segments(1.2, 0.4, 1.4, 1.3)
segments(1.2, 0.3, 1.4, 1.3)
text(x = df1[, 1], y = df1[, 2], labels = df1[, "names"], pos = c(4, 2, 2, 2))
## text(1, 1, "wt", pos = 4)
## text(1.2, 1.2, "A", pos = 2)
## text(1.4, 1.2, "B", pos = 2)


## ------------------------------------------------------------------------
sa <- -0.6
sb <- -0.7
sab <- 0.3
m1b1 <- allFitnessEffects(epistasis = c("A:-B" = sa,
                                        "-A:B" = sb,
                                        "A:B" = sab))
evalAllGenotypes(m1b1, order = FALSE, addwt = TRUE)


## ---- echo=FALSE, fig.width=4, fig.height=4------------------------------

df1 <- data.frame(x = c(1, 1.2, 1.2, 1.4), f = c(1, 1.2, 0.7, 1.5),
                 names = c("wt", "A", "B", "AB"))
plot(df1[, 2] ~ df1[, 1], axes = TRUE, xlab = "", ylab = "Fitness",
     xaxt = "n", yaxt = "n", ylim = c(0.69, 1.53))
segments(1, 1, 1.2, 1.2)
segments(1, 1, 1.2, 0.7)
segments(1.2, 1.2, 1.4, 1.5)
segments(1.2, 0.7, 1.4, 1.5)
text(x = df1[, 1], y = df1[, 2], labels = df1[, "names"], pos = c(3, 3, 3, 2))
## text(1, 1, "wt", pos = 4)
## text(1.2, 1.2, "A", pos = 2)
## text(1.4, 1.2, "B", pos = 2)



## ------------------------------------------------------------------------
sa <- 0.2
sb <- -0.3
sab <- 0.5
m1c1 <- allFitnessEffects(epistasis = c("A:-B" = sa,
                                        "-A:B" = sb,
                                        "A:B" = sab))
evalAllGenotypes(m1c1, order = FALSE, addwt = TRUE)


## ---- echo=FALSE, fig.width=4.5, fig.height=3.5--------------------------

df1 <- data.frame(x = c(1, 2, 3, 4), f = c(1.1, 1, 0.95, 1.2),
                 names = c("u", "wt", "i", "v"))
plot(df1[, 2] ~ df1[, 1], axes = FALSE, xlab = "", ylab = "")
par(las = 1)
axis(2)
axis(1, at = c(1, 2, 3, 4), labels = df1[, "names"], ylab = "")
box()
arrows(c(2, 2, 3), c(1, 1, 0.95),
       c(1, 3, 4), c(1.1, 0.95, 1.2))
## text(1, 1, "wt", pos = 4)
## text(1.2, 1.2, "A", pos = 2)
## text(1.4, 1.2, "B", pos = 2)


## ------------------------------------------------------------------------
su <- 0.1
si <- -0.05
fvi <- 1.2 ## the fitness of the vi mutant
sv <- (fvi/(1 + si)) - 1
sui <- suv <- -1
od <- allFitnessEffects(
    data.frame(parent = c("Root", "Root", "i"),
               child = c("u", "i", "v"),
               s = c(su, si, sv),
               sh = -1,
               typeDep = "MN"),
    epistasis = c(
        "u:i" = sui,
        "u:v" = suv))


## ---- fig.width=3, fig.height=3------------------------------------------
plot(od)


## ------------------------------------------------------------------------
evalAllGenotypes(od, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
u <- 0.1; i <- -0.05; vi <- (1.2/0.95) - 1; ui <- uv <- -Inf
od2 <- allFitnessEffects(
    epistasis = c("u" = u,  "u:i" = ui,
                  "u:v" = uv, "i" = i,
                  "v:-i" = -Inf, "v:i" = vi))
evalAllGenotypes(od2, addwt = TRUE)



## ---- echo=FALSE, fig.width=4, fig.height=3------------------------------

df1 <- data.frame(x = c(1, 2, 3), f = c(1, 0.95, 1.2),
                 names = c("wt", "1", "2"))
plot(df1[, 2] ~ df1[, 1], axes = FALSE, xlab = "", ylab = "")
par(las = 1)
axis(2)
axis(1, at = c(1, 2, 3), labels = df1[, "names"], ylab = "")
box()
segments(c(1, 2), c(1, 0.95),
       c(2, 3), c(0.95, 1.2))
## text(1, 1, "wt", pos = 4)
## text(1.2, 1.2, "A", pos = 2)
## text(1.4, 1.2, "B", pos = 2)




## ---- echo=FALSE, fig.width=4, fig.height=3------------------------------

df1 <- data.frame(x = c(1, 2, 3, 4), f = c(1, 0.95, 0.92, 1.2),
                 names = c("wt", "1", "2", "3"))
plot(df1[, 2] ~ df1[, 1], axes = FALSE, xlab = "", ylab = "")
par(las = 1)
axis(2)
axis(1, at = c(1, 2, 3, 4), labels = df1[, "names"], ylab = "")
box()
segments(c(1, 2, 3), c(1, 0.95, 0.92),
       c(2, 3, 4), c(0.95, 0.92, 1.2))
## text(1, 1, "wt", pos = 4)
## text(1.2, 1.2, "A", pos = 2)
## text(1.4, 1.2, "B", pos = 2)


## ------------------------------------------------------------------------

d1 <- -0.05 ## single mutant fitness 0.95
d2 <- -0.08 ## double mutant fitness 0.92
d3 <- 0.2   ## triple mutant fitness 1.2

s2 <- ((1 + d2)/(1 + d1)^2) - 1
s3 <- ( (1 + d3)/((1 + d1)^3 * (1 + s2)^3) ) - 1

w <- allFitnessEffects(
    data.frame(parent = c("Root", "Root", "Root"),
               child = c("A", "B", "C"),
               s = d1,
               sh = -1,
               typeDep = "MN"),
    epistasis = c(
        "A:B" = s2,
        "A:C" = s2,
        "B:C" = s2,
        "A:B:C" = s3))


## ---- fig.width=4, fig.height=4------------------------------------------
plot(w)


## ------------------------------------------------------------------------
evalAllGenotypes(w, order = FALSE, addwt = TRUE)


## ------------------------------------------------------------------------
wb <- allFitnessEffects(
    epistasis = c(
        "A" = d1,
        "B" = d1,
        "C" = d1,
        "A:B" = s2,
        "A:C" = s2,
        "B:C" = s2,
        "A:B:C" = s3))

evalAllGenotypes(wb, order = FALSE, addwt = TRUE)


## ---- , fig.width=3, fig.height=3----------------------------------------
plot(wb)


## ------------------------------------------------------------------------
wc <- allFitnessEffects(
    epistasis = c(
        "A:-B:-C" = d1,
        "B:-C:-A" = d1,
        "C:-A:-B" = d1,
        "A:B:-C" = d2,
        "A:C:-B" = d2,
        "B:C:-A" = d2,
        "A:B:C" = d3))
evalAllGenotypes(wc, order = FALSE, addwt = TRUE)


## ---- fig.width=4--------------------------------------------------------

pancr <- allFitnessEffects(
    data.frame(parent = c("Root", rep("KRAS", 4), 
                   "SMAD4", "CDNK2A", 
                   "TP53", "TP53", "MLL3"),
               child = c("KRAS","SMAD4", "CDNK2A", 
                   "TP53", "MLL3",
                   rep("PXDN", 3), rep("TGFBR2", 2)),
               s = 0.1,
               sh = -0.9,
               typeDep = "MN"))

plot(pancr)


## ---- fig.height = 4-----------------------------------------------------
rv1 <- allFitnessEffects(data.frame(parent = c("Root", "A", "KRAS"),
                                    child = c("A", "KRAS", "FBXW7"),
                                    s = 0.1,
                                    sh = -0.01,
                                    typeDep = "MN"),
                         geneToModule = c("Root" = "Root",
                             "A" = "EVC2, PIK3CA, TP53",
                             "KRAS" = "KRAS",
                             "FBXW7" = "FBXW7"))

plot(rv1, expandModules = TRUE, autofit = TRUE)



## ---- fig.height=6-------------------------------------------------------
rv2 <- allFitnessEffects(
    data.frame(parent = c("Root", "1", "2", "3", "4"),
               child = c("1", "2", "3", "4", "ELF3"),
               s = 0.1,
               sh = -0.01,
               typeDep = "MN"),
    geneToModule = c("Root" = "Root",
                     "1" = "APC, FBXW7",
                     "2" = "ATM, FAM123B, PIK3CA, TP53",
                     "3" = "BRAF, KRAS, NRAS",
                     "4" = "SMAD2, SMAD4, SOX9",
                     "ELF3" = "ELF3"))

plot(rv2, expandModules = TRUE,   autofit = TRUE)


## ----prbau003, fig.height=6----------------------------------------------

o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("M" = "m",
                              "F" = "f",
                              "D" = "d") )

oneI <- oncoSimulIndiv(o3init, model = "McFL",
                       mu = 5e-5, finalTime = 500,
                       detectionDrivers = 3,
                       onlyCancer = FALSE,
                       initSize = 1000,
                       keepPhylog = TRUE,
                       initMutant = c("m > u > d")
                       )
plotClonePhylog(oneI, N = 0)


## ----prbau003bb, fig.height=6--------------------------------------------
## Note we also disable the stopping stochastically as a function of size
## to allow the population to grow large and generate may different
## clones.
ospI <- oncoSimulPop(2, 
                    o3init, model = "Exp",
                    mu = 5e-5, finalTime = 500,
                    detectionDrivers = 3,
                    onlyCancer = TRUE,
                    initSize = 10,
                    keepPhylog = TRUE,
                    initMutant = c("d > m > z"),
                    mc.cores = 2
                    )
op <- par(mar = rep(0, 4), mfrow = c(1, 2))
plotClonePhylog(ospI[[1]])
plotClonePhylog(ospI[[2]])
par(op)

ossI <- oncoSimulSample(2, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 500,
                        detectionDrivers = 2,
                        onlyCancer = TRUE,
                        initSize = 10,
                        initMutant = c("z > d"),
                        ## check presence of initMutant:
                        thresholdWhole = 1 
                    )

## No phylogeny is kept with oncoSimulSample, but look at the 
## OcurringDrivers and the sample
ossI$popSample
ossI$popSummary[, "OccurringDrivers", drop = FALSE]


## ----prbaux002, eval=FALSE-----------------------------------------------
## gi2 <- rep(0, 5)
## names(gi2) <- letters[1:5]
## oi2 <- allFitnessEffects(noIntGenes = gi2)
## s5 <- oncoSimulPop(200,
##                    oi2,
##                    model = "McFL",
##                    initSize = 1000,
##                    detectionProb = c(p2 = 0.1,
##                                      n2 = 2000,
##                                      PDBaseline = 1000,
##                                      checkSizePEvery = 2),
##                    detectionSize = NA,
##                    finalTime = NA,
##                    keepEvery = NA,
##                    detectionDrivers = NA)
## s5
## hist(unlist(lapply(s5, function(x) x$FinalTime)))


## ------------------------------------------------------------------------
u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf
od2 <- allFitnessEffects(
    epistasis = c("u" = u,  "u:i" = ui,
                  "u:v" = uv, "i" = i,
                  "v:-i" = -Inf, "v:i" = vi))


## ----simul-ochs----------------------------------------------------------
initS <- 20
## We use only a small number of repetitions for the sake
## of speed.
od100 <- oncoSimulPop(10, od2,
                      fixation = c("u", "v, i"),
                      model = "McFL",
                      mu = 1e-4,
                      detectionDrivers = NA,
                      finalTime = NA,
                      detectionSize = NA,
                      detectionProb = NA,
                      onlyCancer = TRUE,
                      initSize = initS, 
                      mc.cores = 2)


## ----ochs-freq-genots----------------------------------------------------
sampledGenotypes(samplePop(od100))                      


## ----sum-simul-ochs------------------------------------------------------
head(summary(od100)[, c(1:3, 8:9)])


## ----fixationG1----------------------------------------------------------
## Create a simple fitness landscape
rl1 <- matrix(0, ncol = 6, nrow = 9)
colnames(rl1) <- c(LETTERS[1:5], "Fitness")
rl1[1, 6] <- 1
rl1[cbind((2:4), c(1:3))] <- 1
rl1[2, 6] <- 1.4
rl1[3, 6] <- 1.32
rl1[4, 6] <- 1.32
rl1[5, ] <- c(0, 1, 0, 0, 1, 1.5)
rl1[6, ] <- c(0, 0, 1, 1, 0, 1.54)
rl1[7, ] <- c(1, 0, 1, 1, 0, 1.65)
rl1[8, ] <- c(1, 1, 1, 1, 0, 1.75)
rl1[9, ] <- c(1, 1, 1, 1, 1, 1.85)
class(rl1) <- c("matrix", "genotype_fitness_matrix")
## plot(rl1) ## to see the fitness landscape

## Gene combinations
local_max_g <- c("A", "B, E", "A, B, C, D, E")
## Specify the genotypes
local_max <- paste0("_,", local_max_g)

fr1 <- allFitnessEffects(genotFitness = rl1, drvNames = LETTERS[1:5])
initS <- 2000


######## Stop on gene combinations   #####
r1 <- oncoSimulPop(10,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = local_max_g, 
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)
sp1 <- samplePop(r1, "last", "singleCell")
sgsp1 <- sampledGenotypes(sp1)
## often you will stop on gene combinations that
## are not local maxima in the fitness landscape
sgsp1
sgsp1$Genotype %in% local_max_g

####### Stop on genotypes   ####

r2 <- oncoSimulPop(10,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = local_max, 
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)
## All final genotypes should be local maxima                       
sp2 <- samplePop(r2, "last", "singleCell")
sgsp2 <- sampledGenotypes(sp2)
sgsp2$Genotype %in% local_max_g





## ----fixationG2----------------------------------------------------------
## Create a simple fitness landscape
rl1 <- matrix(0, ncol = 6, nrow = 9)
colnames(rl1) <- c(LETTERS[1:5], "Fitness")
rl1[1, 6] <- 1
rl1[cbind((2:4), c(1:3))] <- 1
rl1[2, 6] <- 1.4
rl1[3, 6] <- 1.32
rl1[4, 6] <- 1.32
rl1[5, ] <- c(0, 1, 0, 0, 1, 1.5)
rl1[6, ] <- c(0, 0, 1, 1, 0, 1.54)
rl1[7, ] <- c(1, 0, 1, 1, 0, 1.65)
rl1[8, ] <- c(1, 1, 1, 1, 0, 1.75)
rl1[9, ] <- c(1, 1, 1, 1, 1, 1.85)
class(rl1) <- c("matrix", "genotype_fitness_matrix")
## plot(rl1) ## to see the fitness landscape

## The local fitness maxima are
## c("A", "B, E", "A, B, C, D, E")

fr1 <- allFitnessEffects(genotFitness = rl1, drvNames = LETTERS[1:5])
initS <- 2000

## Stop on genotypes

r3 <- oncoSimulPop(10,
                  fp = fr1,
                  model = "McFL",
                  initSize = initS,
                  mu = 1e-4,
                  detectionSize = NA,
                  sampleEvery = .03,
                  keepEvery = 1, 
                  finalTime = 50000,
                  fixation = c(paste0("_,",
                                   c("A", "B, E", "A, B, C, D, E")),
                               fixation_tolerance = 0.1,
                               min_successive_fixation = 200,
                               fixation_min_size = 3000),
                  detectionDrivers = NA,
                  detectionProb = NA,
                  onlyCancer = TRUE,
                  max.num.tries = 500,
                  max.wall.time = 20, 
                  errorHitMaxTries = TRUE,
                  keepPhylog = FALSE,
                  mc.cores = 2)



## ----prbaux001-----------------------------------------------------------
K <- 5
sd <- 0.1
sdp <- 0.15
sp <- 0.05
bauer <- data.frame(parent = c("Root", rep("p", K)),
                    child = c("p", paste0("s", 1:K)),
                    s = c(sd, rep(sdp, K)),
                    sh = c(0, rep(sp, K)),
                    typeDep = "MN")
fbauer <- allFitnessEffects(bauer, drvNames = "p")
set.seed(1)
## Use fairly large mutation rate
b1 <- oncoSimulIndiv(fbauer, mu = 5e-5, initSize = 1000,
                     finalTime = NA,
                     onlyCancer = TRUE,
                     detectionProb = "default")



## ----baux1,fig.width=6.5, fig.height=10, fig.cap="Three drivers' plots of a simulation of Bauer's model"----
par(mfrow = c(3, 1))
## First, drivers
plot(b1, type = "line", addtot = TRUE)
plot(b1, type = "stacked")
plot(b1, type = "stream")


## ----baux2,fig.width=6.5, fig.height=10, fig.cap="Three genotypes' plots of a simulation of Bauer's model"----
par(mfrow = c(3, 1))
## Next, genotypes
plot(b1, show = "genotypes", type = "line")
plot(b1, show = "genotypes", type = "stacked")
plot(b1, show = "genotypes", type = "stream")


## ---- fig.width=6--------------------------------------------------------

set.seed(678)
nd <- 70  
np <- 5000 
s <- 0.1  
sp <- 1e-3 
spp <- -sp/(1 + sp)
mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                          drvNames = seq.int(nd))
mcf1s <-  oncoSimulIndiv(mcf1,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionProb = "default",
                         detectionSize = NA, 
                         detectionDrivers = NA,
                         sampleEvery = 0.025,
                         keepEvery = 8,
                         initSize = 2000,
                         finalTime = 4000,
                         onlyCancer = FALSE)
summary(mcf1s)



## ----mcf1sx1,fig.width=6.5, fig.height=10--------------------------------
par(mfrow  = c(2, 1))
## I use thinData to make figures smaller and faster
plot(mcf1s, addtot = TRUE, lwdClone = 0.9, log = "", 
     thinData = TRUE, thinData.keep = 0.5)
plot(mcf1s, show = "drivers", type = "stacked",
     thinData = TRUE, thinData.keep = 0.3,
     legend.ncols = 2)


## ---- eval=FALSE---------------------------------------------------------
## 
## set.seed(123)
## nd <- 70
## np <- 50000
## s <- 0.1
## sp <- 1e-4 ## as we have many more passengers
## spp <- -sp/(1 + sp)
## mcfL <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
##                           drvNames = seq.int(nd))
## mcfLs <-  oncoSimulIndiv(mcfL,
##                          model = "McFL",
##                          mu = 1e-7,
##                          detectionSize = 1e8,
##                          detectionDrivers = 100,
##                          sampleEvery = 0.02,
##                          keepEvery = 2,
##                          initSize = 1000,
##                          finalTime = 2000,
##                          onlyCancer = FALSE)


## ---- mcflsx2,fig.width=6------------------------------------------------
data(mcfLs)
plot(mcfLs, addtot = TRUE, lwdClone = 0.9, log = "",
     thinData = TRUE, thinData.keep = 0.3,
     plotDiversity = TRUE)


## ------------------------------------------------------------------------
summary(mcfLs)
## number of passengers per clone
summary(colSums(mcfLs$Genotypes[-(1:70), ]))


## ----mcflsx3-------------------------------------------------------------
plot(mcfLs, type = "stacked", thinData = TRUE, 
     thinData.keep = 0.2,
     plotDiversity = TRUE,
     xlim = c(0, 1000))


## ------------------------------------------------------------------------
data(examplesFitnessEffects)
names(examplesFitnessEffects)


## ----smmcfls-------------------------------------------------------------
data(examplesFitnessEffects)
evalAllGenotypes(examplesFitnessEffects$cbn1, order = FALSE)[1:10, ]
sm <-  oncoSimulIndiv(examplesFitnessEffects$cbn1,
                      model = "McFL", 
                      mu = 5e-7,
                      detectionSize = 1e8, 
                      detectionDrivers = 2,
                      detectionProb = "default",
                      sampleEvery = 0.025,
                      keepEvery = 5,
                      initSize = 2000,
                      onlyCancer = TRUE)
summary(sm)


## ----smbosb1-------------------------------------------------------------
set.seed(1234)
evalAllGenotypes(examplesFitnessEffects$cbn1, order = FALSE, 
                 model = "Bozic")[1:10, ]
sb <-  oncoSimulIndiv(examplesFitnessEffects$cbn1,
                       model = "Bozic", 
                       mu = 5e-6,
                      detectionProb = "default",
                       detectionSize = 1e8, 
                       detectionDrivers = 4,
                       sampleEvery = 2,
                       initSize = 2000,
                       onlyCancer = TRUE)
summary(sb)


## ----sbx1,fig.width=6.5, fig.height=3.3----------------------------------
## Show drivers, line plot
par(cex = 0.75, las = 1)
plot(sb,show = "drivers", type = "line", addtot = TRUE,
     plotDiversity = TRUE)

## ----sbx2,fig.width=6.5, fig.height=3.3----------------------------------
## Drivers, stacked
par(cex = 0.75, las = 1)
plot(sb,show = "drivers", type = "stacked", plotDiversity = TRUE)

## ----sbx3,fig.width=6.5, fig.height=3.3----------------------------------
## Drivers, stream
par(cex = 0.75, las = 1)
plot(sb,show = "drivers", type = "stream", plotDiversity = TRUE)


## ----sbx4,fig.width=6.5, fig.height=3.3----------------------------------
## Genotypes, line plot
par(cex = 0.75, las = 1)
plot(sb,show = "genotypes", type = "line", plotDiversity = TRUE)

## ----sbx5,fig.width=6.5, fig.height=3.3----------------------------------
## Genotypes, stacked
par(cex = 0.75, las = 1)
plot(sb,show = "genotypes", type = "stacked", plotDiversity = TRUE)

## ----sbx6,fig.width=6.5, fig.height=3.3----------------------------------
## Genotypes, stream
par(cex = 0.75, las = 1)
plot(sb,show = "genotypes", type = "stream", plotDiversity = TRUE)


## ---- fig.width=6--------------------------------------------------------

set.seed(4321)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 6000,
                       onlyCancer = FALSE) 


## ----tmpmx1,fig.width=6.5, fig.height=4.1--------------------------------
par(las = 1, cex = 0.85)
plot(tmp, addtot = TRUE, log = "", plotDiversity = TRUE,
     thinData = TRUE, thinData.keep = 0.2)

## ----tmpmx2,fig.width=6.5, fig.height=4.1--------------------------------
par(las = 1, cex = 0.85)
plot(tmp, type = "stacked", plotDiversity = TRUE, 
     ylim = c(0, 5500), legend.ncols = 4,
     thinData = TRUE, thinData.keep = 0.2)


## ------------------------------------------------------------------------
evalAllGenotypes(examplesFitnessEffects[["o3"]], addwt = TRUE,
                 order = TRUE)


## ----tmpmx3,fig.width=6.5, fig.height=5.5--------------------------------
plot(tmp, show = "genotypes", ylim = c(0, 5500), legend.ncols = 3,
     thinData = TRUE, thinData.keep = 0.5)


## ------------------------------------------------------------------------
set.seed(15)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 20000,
                       onlyCancer = FALSE,
                       extraTime = 1500)
tmp


## ----tmpmdx5,fig.width=6.5, fig.height=4---------------------------------
par(las = 1, cex = 0.85)
plot(tmp, addtot = TRUE, log = "", plotDiversity = TRUE,
     thinData = TRUE, thinData.keep = 0.5)

## ----tmpmdx6,fig.width=6.5, fig.height=4---------------------------------
par(las = 1, cex = 0.85)
plot(tmp, type = "stacked", plotDiversity = TRUE,
     legend.ncols = 4, ylim = c(0, 5200), xlim = c(3400, 5000),
     thinData = TRUE, thinData.keep = 0.5)


## ----tmpmdx7,fig.width=6.5, fig.height=5.3-------------------------------
## Improve telling apart the most abundant 
## genotypes by sorting colors
## differently via breakSortColors
## Modify ncols of legend, so it is legible by not overlapping
## with plot
par(las = 1, cex = 0.85)
plot(tmp, show = "genotypes", breakSortColors = "distave",
     plotDiversity = TRUE, legend.ncols = 4,
     ylim = c(0, 5300), xlim = c(3400, 5000),
     thinData = TRUE, thinData.keep = 0.5)


## ---- eval=FALSE---------------------------------------------------------
## ## Convert the data
## lb1 <- OncoSimulWide2Long(b1)
## 
## ## Install the streamgraph package from GitHub and load
## library(devtools)
## devtools::install_github("hrbrmstr/streamgraph")
## library(streamgraph)
## 
## ## Stream plot for Genotypes
## sg_legend(streamgraph(lb1, Genotype, Y, Time, scale = "continuous"),
##               show=TRUE, label="Genotype: ")
## 
## ## Staked area plot and we use the pipe
## streamgraph(lb1, Genotype, Y, Time, scale = "continuous",
##             offset = "zero") %>%
##     sg_legend(show=TRUE, label="Genotype: ")


## ----pancrpopcreate------------------------------------------------------

pancrPop <- oncoSimulPop(4, pancr,
                         detectionSize = 1e7,
                         keepEvery = 10,
                         mc.cores = 2)

summary(pancrPop)
samplePop(pancrPop)



## ------------------------------------------------------------------------
library(parallel)

if(.Platform$OS.type == "windows") {
    mc.cores <- 1
} else {
    mc.cores <- 2
}

p2 <- mclapply(1:4, function(x) oncoSimulIndiv(pancr,
                                               detectionSize = 1e7,
                                               keepEvery = 10),
                                               mc.cores = mc.cores)
class(p2) <- "oncosimulpop"
samplePop(p2)


## ------------------------------------------------------------------------
tail(pancrPop[[1]]$pops.by.time)


## ------------------------------------------------------------------------
pancrPopNH <- oncoSimulPop(4, pancr,
                           detectionSize = 1e7,
                           keepEvery = NA,
                           mc.cores = 2)

summary(pancrPopNH)
samplePop(pancrPopNH)


## ------------------------------------------------------------------------
pancrPopNH[[1]]$pops.by.time


## ------------------------------------------------------------------------
pancrSamp <- oncoSimulSample(4, pancr)
pancrSamp$popSamp



## ------------------------------------------------------------------------

set.seed(15)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 20000,
                       onlyCancer = FALSE,
                       extraTime = 1500,
                       keepPhylog = TRUE)
tmp


## ------------------------------------------------------------------------
plotClonePhylog(tmp, N = 0)


## ------------------------------------------------------------------------
plotClonePhylog(tmp, N = 1)


## ----pcpkeepx1-----------------------------------------------------------
plotClonePhylog(tmp, N = 1, keepEvents = TRUE)


## ------------------------------------------------------------------------
plotClonePhylog(tmp, N = 1, timeEvents = TRUE)


## ---- fig.keep="none"----------------------------------------------------
get.adjacency(plotClonePhylog(tmp, N = 1, returnGraph = TRUE))



## ------------------------------------------------------------------------

set.seed(456)
mcf1s <-  oncoSimulIndiv(mcf1,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionSize = 1e8, 
                         detectionDrivers = 100,
                         sampleEvery = 0.025,
                         keepEvery = 2,
                         initSize = 2000,
                         finalTime = 1000,
                         onlyCancer = FALSE,
                         keepPhylog = TRUE)



## ------------------------------------------------------------------------
plotClonePhylog(mcf1s, N = 1)


## ------------------------------------------------------------------------
par(cex = 0.7)
plotClonePhylog(mcf1s, N = 1, timeEvents = TRUE)


## ------------------------------------------------------------------------
par(cex = 0.7)
plotClonePhylog(mcf1s, N = 1, t = c(800, 1000))


## ------------------------------------------------------------------------
par(cex = 0.7)
plotClonePhylog(mcf1s, N = 1, t = c(900, 1000), timeEvents = TRUE)


## ----fig.keep="none"-----------------------------------------------------
g1 <- plotClonePhylog(mcf1s, N = 1, t = c(900, 1000),
                      returnGraph = TRUE)


## ------------------------------------------------------------------------
plot(g1)


## ---- eval=FALSE---------------------------------------------------------
## op <- par(ask = TRUE)
## while(TRUE) {
##     tmp <- oncoSimulIndiv(smn1, model = "McFL",
##                           mu = 5e-5, finalTime = 500,
##                           detectionDrivers = 3,
##                           onlyCancer = FALSE,
##                           initSize = 1000, keepPhylog = TRUE)
##     plotClonePhylog(tmp, N = 0)
## }
## par(op)


## ------------------------------------------------------------------------

oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("F" = "f1, f2, f3",
                                "D" = "d1, d2") )
oiI1 <- oncoSimulIndiv(oi, model = "Exp")
oiP1 <- oncoSimulPop(4, oi,
                     keepEvery = 10,
                     mc.cores = 2,
                     keepPhylog = TRUE)



## ---- fig.height=9-------------------------------------------------------

op <- par(mar = rep(0, 4), mfrow = c(2, 1))
plotClonePhylog(oiP1[[1]])
plotClonePhylog(oiP1[[2]])
par(op)



## ------------------------------------------------------------------------
## A small example
rfitness(3)

## A 5-gene example, where the reference genotype is the
## one with all positions mutated, similar to Greene and Crona,
## 2014.  We will plot the landscape and use it for simulations
## We downplay the random component with a sd = 0.5

r1 <- rfitness(5, reference = rep(1, 5), sd = 0.6)
plot(r1)
oncoSimulIndiv(allFitnessEffects(genotFitness = r1))


## ----nkex1---------------------------------------------------------------
rnk <- rfitness(5, K = 3, model = "NK")
plot(rnk)
oncoSimulIndiv(allFitnessEffects(genotFitness = rnk))


## ----addex1--------------------------------------------------------------
radd <- rfitness(4, model = "Additive", mu = 0, sd = 0.5)
plot(radd)


## ----eggex1--------------------------------------------------------------
regg <- rfitness(4, model = "Eggbox", e = 1, E = 0.5)
plot(regg)


## ----isingex1------------------------------------------------------------
ris <- rfitness(g = 4, model = "Ising", i = 0.002, I = 0.45)
plot(ris)


## ----fullex1-------------------------------------------------------------
rnk <- rfitness(4, model = "Full", i = 0.5, I = 0.1,
                e = 2, s = 0.3, S = 0.08)
plot(rnk)


## ----magstats1-----------------------------------------------------------
rnk1 <- rfitness(6, K = 1, model = "NK")
Magellan_stats(rnk1)

rnk2 <- rfitness(6, K = 4, model = "NK")
Magellan_stats(rnk2)


## ----fdf1----------------------------------------------------------------
## Define fitness of the different genotypes
gffd <- data.frame(
    Genotype = c("WT", "A", "B", "C", "A, B"), 
    Fitness = c("1 + 1.5 * f_1_2",
                "1.3 + 1.5 * f_1_2",
                "1.4",
                "1.1 + 0.7*((f_1 + f_2) > 0.3) + f_1_2",
                "1.2 + sqrt(f_1 + f_3 + f_2) - 0.3 * (f_1_2 > 0.5)"),
    stringsAsFactors = FALSE)


## ----fdf1b---------------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 20, 20, 30, 0)))

evalAllGenotypes(allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 30, 40, 0, 10)))

evalAllGenotypes(allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 30, 40, 0, 100)))



## ----fdf1c---------------------------------------------------------------
afd <- allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

set.seed(1) ## for reproducibility
sfd <- oncoSimulIndiv(afd, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = FALSE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(sfd, show = "genotypes")


## ----hurlbutpay, eval=TRUE,echo=FALSE, fig.cap="Payoff matrix from Table 2 of Hurlbut et al., 2018, 'Game Theoretical Model of Cancer Dynamics withFour Cell Phenotypes', *Games*, 9, 61; doi:10.3390/g9030061."----
knitr::include_graphics("hurlbut.png")


## ----hurlbutfit1---------------------------------------------------------
options(stringsAsFactors = FALSE) ## Get rid of the messages

create_fe <- function(a, b, c, d, e, f, g,
                        gt = c("WT", "A", "P", "C")) {
  data.frame(Genotype = gt,
             Fitness = c(
                 paste0("1 + ",
                       d, " * f_1 ",
                       "- ", c, " * f_3"),
                 paste0("1 - ", a, 
				     " + ", d, " + ",
                       f, " * f_1 ",
                      "- ", c, " * f_3"),
                 paste0("1 + ", g, " + ",
                       d, " * f_1 ",
                       "- ", c, " * (1 + ",
                       g, ") * f_3"),
                 paste0("1 - ", b, " + ",
                       e, " * f_ + ",
                       "(", d, " + ", 
					   e, ") * f_1 + ",
                       e , " * f_2")),
             stringsAsFactors = FALSE)
}


## ----hbf1check-----------------------------------------------------------
create_fe("a", "b", "c", "d", "e", "f", "g")


## ----hurlbutfit2---------------------------------------------------------
## Different assumption about origins from mutation:
## WT -> P; P -> A,P; P -> C,P

create_fe2 <- function(a, b, c, d, e, f, g,
                        gt = c("WT", "A", "P", "C", "A, P", "A, C",
                                "C, P")) {
  data.frame(Genotype = gt,
             Fitness = c(
                 paste0("1 + ",
                       d, " * f_1_2 ",
                       "- ", c, " * f_2_3"),
                 "0",
                 paste0("1 + ", g, " + ",
                       d, " * f_1_2 ",
                       "- ", c, " * (1 + ",
                       g, ") * f_2_3"),
                 "0",
                 paste0("1 - ", a, " + ", 
				 d, " + ",
                       f, " * f_1_2 ",
                       "- ", c, " * f_2_3"),
                 "0",
                 paste0("1 - ", b, " + ",
                       e, " * f_ + ",
                       "(", d, " + ", 
					   e, ") * f_1_2 + ",
                       e , " * f_2")),
             stringsAsFactors = FALSE)
}

## And check:
create_fe2("a", "b", "c", "d", "e", "f", "g")


## ----hurl3a, message=FALSE-----------------------------------------------
## Figure 3a
afe_3_a <- allFitnessEffects(
		        genotFitness =
                       create_fe(0.02, 0.04, 0.08, 0.06,
                                 0.15, 0.1, 0.06),
                frequencyDependentFitness = TRUE,
                frequencyType = "rel")
set.seed(2)
s_3_a <- oncoSimulIndiv(afe_3_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
## plot(s_3_a, show = "genotypes",
##      xlim = c(40, 200),
##      col = c("black", "green", "red", "blue"))
plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))


## ----hurl3b, message=FALSE-----------------------------------------------
## Figure 3b
afe_3_b <- allFitnessEffects(
                genotFitness =
                       create_fe(0.02, 0.04, 0.08, 0.1,
                                 0.15, 0.1, 0.05),
                frequencyDependentFitness = TRUE,
                frequencyType = "rel")
set.seed(2)
s_3_b <- oncoSimulIndiv(afe_3_b,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
## plot(s_3_b, show = "genotypes", 
##      col = c("black", "green", "red", "blue"))
plot(s_3_b, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))


## ----hurl3b2-------------------------------------------------------------
## Figure 3b. Now with WT -> P; P -> A,P; P -> C,P
afe_3_b_2 <- allFitnessEffects(
                  genotFitness =
                      create_fe2(0.02, 0.04, 0.08, 0.1,
                                 0.15, 0.1, 0.05),
                  frequencyDependentFitness = TRUE,
                  frequencyType = "rel")
set.seed(2)
s_3_b_2 <- oncoSimulIndiv(afe_3_b_2,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 300,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE)
plot(s_3_b_2, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))



## ----fdf2, message=FALSE-------------------------------------------------
options(stringsAsFactors = FALSE) ## get rid of the messages
gffd3 <- data.frame(Genotype = c("WT", "A", "B"), 
                   Fitness = c("1",
                               "1 + 0.2 * (n_2 > 10)",
                               ".9 + 0.4 * (n_1 > 10)"
                               ))
afd3 <- allFitnessEffects(genotFitness = gffd3, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")


## ----fdf2b, message=FALSE------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness = gffd3, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs",
                         spPopSizes = c(100, 1, 11)))

evalAllGenotypes(allFitnessEffects(genotFitness = gffd3, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs",
                         spPopSizes = c(100, 11, 1)))


## ----fdf2c---------------------------------------------------------------
set.seed(1)
sfd3 <- oncoSimulIndiv(afd3,
                       model = "McFLD", 
                       onlyCancer = FALSE, 
                       finalTime = 200,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)


## ----echo=FALSE----------------------------------------------------------
op <- par(mfrow = c(1, 2))


## ------------------------------------------------------------------------
plot(sfd3, show = "genotypes", type = "line")
plot(sfd3, show = "genotypes")
sfd3


## ----echo=FALSE----------------------------------------------------------
par(op)



## ----fdf2d---------------------------------------------------------------
set.seed(1)
sfd4 <- oncoSimulIndiv(afd3,
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 200,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)


## ---- echo=FALSE---------------------------------------------------------
op <- par(mfrow = c(1, 2))


## ----polotfdfv6----------------------------------------------------------
plot(sfd4, show = "genotypes", type = "line")
plot(sfd4, show = "genotypes")
sfd4


## ---- echo=FALSE---------------------------------------------------------
par(op)


## ----fdfpopfinal---------------------------------------------------------
## Check final pop size corresponds to birth = death 
K <- 5000/(exp(1) - 1) 
K
log1p(4290/K)


## ----lotka1--------------------------------------------------------------
G_fe_LV <- function(r1, r2, K1, K2, a_12, a_21, awt = 1e-4,
                                 gt = c("WT", "S1", "S2")) {
    data.frame(Genotype = gt,
               Fitness = c(
                  paste0("max(0.1, 1 - ", awt, " * (n_2 + n_1))"),
                  paste0("1 + ", r1,
                         " * ( 1 - (n_1 + ", a_12, " * n_2)/", K1,
                         ")"),
                  paste0("1 + ", r2,
                         " * ( 1 - (n_2 + ", a_21, " * n_1)/", K2,
                         ")")
                  ))
}

## Show expressions for birth rates
G_fe_LV("r1", "r2", "K1", "K2", "a_12", "a_21", "awt")



## ---- echo = FALSE-------------------------------------------------------
set.seed(1)

## ----complv1, message=FALSE----------------------------------------------
fe_competition <-
    allFitnessEffects(
        genotFitness =
            G_fe_LV(1.5, 1.4, 10000, 4000, 0.6, 0.2,
                  gt = c("WT","S1", "S2")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")

competition <- oncoSimulIndiv(fe_competition,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 40000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)


## ----pp1wt---------------------------------------------------------------
plot(competition, show = "genotypes")


## ---- echo=FALSE---------------------------------------------------------
op <- par(mfrow = c(1, 2))


## ----pp1nowt-------------------------------------------------------------
plot(competition, show = "genotypes",
     xlim = c(80, 100))
plot(competition, show = "genotypes", type = "line",
     xlim = c(80, 100), ylim = c(1500, 12000))


## ---- echo=FALSE---------------------------------------------------------
par(op)	 


## ----echo = FALSE--------------------------------------------------------
set.seed(1)


## ----lotka2--------------------------------------------------------------

fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            G_fe_LV(1.5, 1.4, 10000, 4000, 0.6, -0.5, awt = 1,
                  gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")

pred_prey <- oncoSimulIndiv(fe_pred_prey,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 40000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)


## ----echo=FALSE----------------------------------------------------------
op <- par(mfrow = c(1, 2))


## ----prepreylv2----------------------------------------------------------
plot(pred_prey, show = "genotypes")
plot(pred_prey, show = "genotypes",
     xlim = c(50, 100))


## ----echo=FALSE----------------------------------------------------------
op <- par(mfrow = c(1, 2))


## ----echo=FALSE----------------------------------------------------------
set.seed(2)


## ----ppsmallk, message=FALSE---------------------------------------------
fe_pred_prey <-
    allFitnessEffects(
        genotFitness =
            G_fe_LV(1.5, 1.4, 100, 40, 0.6, -0.5, awt = 1,
                  gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")


pred_prey <- oncoSimulIndiv(fe_pred_prey,
                            model = "Exp",
                            onlyCancer = FALSE, 
                            finalTime = 200,
                            mu = 1e-3,
                            initSize = 1000, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)


## ----echo=FALSE----------------------------------------------------------
op <- par(mfrow = c(1, 2))


## ----prepreylv7----------------------------------------------------------
plot(pred_prey, show = "genotypes")
plot(pred_prey, show = "genotypes",
     xlim = c(50, 200))



## ----echo=FALSE----------------------------------------------------------
par(op)	 


## ----checklvpredprey-----------------------------------------------------
evalAllGenotypes(allFitnessEffects(
    genotFitness =
        G_fe_LV(1.5, 1.4, 100, 40,
              0.6, -0.5, awt = 0.1,
              gt = c("WT","prey", "Predator")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs",
    spPopSizes = c(0, 0, 20)))


evalAllGenotypes(allFitnessEffects(
    genotFitness =
        G_fe_LV(1.5, 1.4, 100, 40,
              0.6, -0.5, awt = 0.1,
              gt = c("WT","prey", "Predator")),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs",
    spPopSizes = c(0, 0, 40)))


## ----predprey2a----------------------------------------------------------
## Use e for epsilon and d for delta
C_fe_pred_prey <- function(r, a, c, e, d, awt = 0.1,
                           gt = c("WT", "prey", "Predator")) {
    data.frame(Genotype = gt,
               Fitness = c(
                   paste0("max(0.1, 1 - ", awt,
                          " * (n_2 + n_1))"),
                   paste0("1 + ", r, " - ", a,
                          " * ", c, " * n_2"),
                   paste0("1 + ", e, " * ", a,
                          " * ", c, " * n_1 - ", d)
               ))
}

C_fe_pred_prey("r", "a", "c", "e", "d")

## ----echo=FALSE----------------------------------------------------------
set.seed(2)

## ----predprey2b----------------------------------------------------------
fe_pred_prey2 <-
    allFitnessEffects(
        genotFitness =
            C_fe_pred_prey(r = .7, a = 1, c = 0.005,
                           e = 0.02, d = 0.4, awt = 0.001,
                           gt = c("WT","prey", "Predator")),
        frequencyDependentFitness = TRUE,
        frequencyType = "abs")
pred_prey2 <- oncoSimulIndiv(fe_pred_prey2,
                            model = "Exp",
                            sampleEvery = 0.01,
                            mu = 1e-3,
                            onlyCancer = FALSE, 
                            finalTime = 80,
                            initSize = 1e4, 
                            keepPhylog = TRUE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE)		
plot(pred_prey2, show = "genotypes")


## ----commens, eval=FALSE-------------------------------------------------
## fe_commens <-
##     allFitnessEffects(
##         genotFitness =
##             G_fe_LV(1.2, 1.3, 5000, 20000,
##                                  0, -0.2,
##                                  gt = c("WT","A", "Commensal")),
##         frequencyDependentFitness = TRUE,
##         frequencyType = "abs")
## 
## commens <- oncoSimulIndiv(fe_commens,
##                             model = "Exp",
##                             onlyCancer = FALSE,
##                             finalTime = 100,
##                             mu = 1e-4,
##                             initSize = 40000,
##                             keepPhylog = TRUE,
##                             seed = NULL,
##                             errorHitMaxTries = FALSE,
##                             errorHitWallTime = FALSE)
## 
## plot(commens, show = "genotypes")
## 
## plot(commens, show = "genotypes",
##      xlim = c(80, 100))
## 
## plot(commens, show = "genotypes", type = "line",
##      xlim = c(80, 100), ylim = c(2000, 22000))
## 


## ----fdfar, message=FALSE------------------------------------------------
rar <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                 Fitness = c("1",
                             "1.1 + .3*f_2",
                             "1.2 + .4*f_1",
                             "1.0 + .5 * (f_1 + f_2)"))
afear <- allFitnessEffects(genotFitness = rar, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel",
                         spPopSizes = c(100, 200, 300, 400))
evalAllGenotypes(afear)


rar2 <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                 Fitness = c("1",
                             "1.1 + .3*(n_2/N)",
                             "1.2 + .4*(n_1/N)",
                             "1.0 + .5 * ((n_1 + n_2)/N)"))
afear2 <- allFitnessEffects(genotFitness = rar2, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs",
                         spPopSizes = c(100, 200, 300, 400))
evalAllGenotypes(afear2)


## ----relarres, message=FALSE---------------------------------------------
set.seed(1)
tmp1 <- oncoSimulIndiv(afear, 
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 30,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)

set.seed(1)
tmp2 <- oncoSimulIndiv(afear2, 
                       model = "McFL", 
                       onlyCancer = FALSE, 
                       finalTime = 30,
                       mu = 1e-4,
                       initSize = 5000, 
                       keepPhylog = FALSE,
                       seed = NULL, 
                       errorHitMaxTries = FALSE, 
                       errorHitWallTime = FALSE)
stopifnot(identical(print(tmp1), print(tmp2)))


## ----relar3a, message=FALSE----------------------------------------------
rar3 <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                 Fitness = c("1",
                             "1.1 + .3*(n_2/N)",
                             "1.2 + .4*(n_1/N)",
                             "1.0 + .5 * ( n_1 > 20)"))
afear3 <- allFitnessEffects(genotFitness = rar3, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs",
                         spPopSizes = c(100, 200, 300, 400))
evalAllGenotypes(afear3)


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


## ----fdfmutex------------------------------------------------------------
## Relative
r1fd <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                 Fitness = c("1",
                             "1.4 + 1*(f_2)",
                             "1.4 + 1*(f_1)",
                             "1.6 + f_1 + f_2"))
afe4 <- allFitnessEffects(genotFitness = r1fd, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")


set.seed(1)
s1fd <- oncoSimulIndiv(afe4, 
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 50,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s1fd, show = "genotypes")


mtfd <- allMutatorEffects(epistasis = c("A" = 0.1,
                                      "B" = 10))
set.seed(1)
s2fd <- oncoSimulIndiv(afe4,
                     muEF = mtfd,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 50,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s2fd, show = "genotypes")


## ----echo=FALSE----------------------------------------------------------
op <- par(mfrow = c(1, 2))


## ----figmutfdf-----------------------------------------------------------
plotClonePhylog(s1fd, keepEvents = TRUE)
plotClonePhylog(s2fd, keepEvents = TRUE)


## ----echo=FALSE----------------------------------------------------------
par(op)


## ----exmutfdf2, eval=FALSE-----------------------------------------------
## 
## ## Absolute
## r5 <- data.frame(Genotype = c("WT", "A", "B", "A, B"),
##                  Fitness = c("1",
##                              "1.25 - .0025*(n_2)",
##                              "1.25 - .0025*(n_1)",
##                              "1.4"))
## afe5 <- allFitnessEffects(genotFitness = r5,
##                          frequencyDependentFitness = TRUE,
##                          frequencyType = "abs")
## set.seed(8)
## s5 <- oncoSimulIndiv(afe5,
##                      model = "McFL",
##                      onlyCancer = FALSE,
##                      finalTime = 100,
##                      mu = 1e-4,
##                      initSize = 5000,
##                      keepPhylog = TRUE,
##                      seed = NULL,
##                      errorHitMaxTries = FALSE,
##                      errorHitWallTime = FALSE)
## plot(s5, show = "genotypes")
## plot(s5, show = "genotypes", log = "y", type = "line")
## 
## mt <- allMutatorEffects(epistasis = c("A" = 0.1,
##                                       "B" = 10))
## set.seed(8)
## s6 <- oncoSimulIndiv(afe5,
##                      muEF = mt,
##                      model = "McFL",
##                      onlyCancer = FALSE,
##                      finalTime = 100,
##                      mu = 1e-4,
##                      initSize = 5000,
##                      keepPhylog = TRUE,
##                      seed = NULL,
##                      errorHitMaxTries = FALSE,
##                      errorHitWallTime = FALSE)
## plot(s6, show = "genotypes")
## plot(s6, show = "genotypes", log = "y", type = "line")
## 
## plotClonePhylog(s5, keepEvents = TRUE)
## plotClonePhylog(s6, keepEvents = TRUE)
## 


## ----noworkeval----------------------------------------------------------

evalAllGenotypes(allFitnessEffects(genotFitness = r1fd, 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(10, 20, 30, 40)))

## Fitness is wrong
evalAllGenotypesFitAndMut(allFitnessEffects(genotFitness = r1fd, 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(10, 20, 30, 40)),
                    mtfd)


## ----message=F-----------------------------------------------------------
crs <- function (a, b, c){
  data.frame(Genotype = c("WT", "C", "R"),
             Fitness = c(paste0("1 + ", a, " * f_2 - ", b, " * f_1"),
                         paste0("1 + ", b, " * f_ - ", c, " * f_2"),
                         paste0("1 + ", c, " * f_1 - ", a, " * f_")
             ))
}


## ----message=F-----------------------------------------------------------
crs("a", "b", "c")


## ----message=F-----------------------------------------------------------
afcrs1 <- allFitnessEffects(genotFitness = crs(1, 1, 1), 
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

resultscrs1 <- oncoSimulIndiv(afcrs1,
                             model = "McFL", 
                             onlyCancer = FALSE,
                             finalTime = 100, 
                             mu = 1e-2,
                             initSize = 4000, 
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
op <- par(mfrow = c(1, 2))
plot(resultscrs1, show = "genotypes", type = "line", cex.lab=1.1,
     las = 1)
plot(resultscrs1, show = "genotypes", type = "stacked")
par(op)	 


## ----message=F-----------------------------------------------------------
afcrs2 <- allFitnessEffects(genotFitness = crs(10, 1, 1), 
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")




## ----eval=FALSE----------------------------------------------------------
## resultscrs2 <- oncoSimulPop(10,
##                            afcrs2,
##                              model = "McFL",
##                              onlyCancer = FALSE,
##                              finalTime = 100,
##                              mu = 1e-2,
##                              initSize = 4000,
##                              keepPhylog = TRUE,
##                              seed = NULL,
##                              errorHitMaxTries = FALSE,
##                              errorHitWallTime = FALSE)


## ----message=F-----------------------------------------------------------
set.seed(1)

resultscrs2a <- oncoSimulIndiv(afcrs2,
                             model = "McFL", 
                             onlyCancer = FALSE,
                             finalTime = 100, 
                             mu = 1e-2,
                             initSize = 4000, 
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)

plot(resultscrs2a, show = "genotypes", type = "line")


## ----message=F-----------------------------------------------------------
set.seed(3)

resultscrs2b <- oncoSimulIndiv(afcrs2,
                             model = "McFL", 
                             onlyCancer = FALSE,
                             finalTime = 60, 
                             mu = 1e-2,
                             initSize = 4000, 
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)

plot(resultscrs2b, show = "genotypes", type = "line", cex.lab=1.1,
     las = 1)


## ----message=F-----------------------------------------------------------
afcrs3 <- allFitnessEffects(genotFitness = crs(1, 5, 5), 
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

resultscrs3 <- oncoSimulIndiv(afcrs3,
                             model = "McFL", 
                             onlyCancer = FALSE,
                             finalTime = 60, 
                             mu = 1e-2,
                             initSize = 4000, 
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)

plot(resultscrs3, show = "genotypes", type = "line", cex.lab=1.1,
     las = 1)


## ---- message=F----------------------------------------------------------
## Stablish Genotype-Fitnees mapping. D = Dove, H = Hawk
H_D_fitness <- function(c, v,
                    gt = c("WT", "H", "D")) {
  data.frame(Genotype = gt,
             Fitness = c(
               paste0("1"),
               paste0("1 + f_1 *", (v-c)/2, "+ f_2 *", v),
               paste0("1 + f_2 *", v/2)))
}

## Fitness Effects specification
HD_competition <-allFitnessEffects(
  genotFitness = H_D_fitness(10, 2, 
                         gt = c("WT", "H", "D")),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel",
  spPopSizes = c(5000, 5000, 5000))

## Plot fitness landscape of genotype "H, D" evaluation
data.frame("Doves_fitness" = evalGenotype(genotype = "D",
                                          fitnessEffects = HD_competition), 
           "Hawks_fitness" = evalGenotype(genotype = "H", 
                                          fitnessEffects = HD_competition)) 


## ---- message=F----------------------------------------------------------
## Simulated trajectories
## run only a few for the sake of speed
simulation <- oncoSimulPop(2,
                           mc.cores = 2,
                           HD_competition,
                           model = "McFL", # There is no collapse
                           onlyCancer = FALSE,
                           finalTime = 50,
                           mu = 1e-2, # Quick emergence of D and H
                           initSize = 4000,
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)



## Plot first trajectory as an example
plot(simulation[[1]], show = "genotypes", type = "line", 
    xlim = c(40, 50),
     lwdClone = 2, ylab = "Number of individuals", 
     main = "Hawk and Dove trajectory",  
     col = c("#a37acc", "#f8776d", "#7daf00"),
     font.main=2, font.lab=2,
     cex.main=1.4, cex.lab=1.1,
     las = 1)



## ------------------------------------------------------------------------
## Recover the final result from first simulation
result <- tail(simulation[[1]][[1]], 1)

## Get the number of organisms from each species
n_WT <- result[2]
n_D <- result[3]
n_H <- result[4]
total <- n_WT + n_D + n_H

## Dove percentage
data.frame("Doves" = round(n_D/total, 2)*100, 
           "Hawks" = round(n_H/total, 2)*100  )


## ---- message=FALSE------------------------------------------------------

# Definition of the function for creating the corresponding dataframe.
avc <- function (a, v, c) {
  data.frame(Genotype = c("WT", "GLY", "VOP", "DEF"),
             Fitness = c("1",
                         paste0("1 + ",a," * (f_1 + 1)"),
                         paste0("1 + ",a," * f_1 + ",v," * (f_2 + 1) - ",c),
                         paste0("1 + ",a," * f_1 + ",v," * f_2")
                         ))
                          }

# Specification of the different effects on fitness.
afavc <- allFitnessEffects(genotFitness = avc(2.5, 2, 1),
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

## For real, you would probably want to run
## this multiple times with oncoSimulPop
simulation <- oncoSimulIndiv(afavc,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 15,
                           mu = 1e-3,
                           initSize = 4000,
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)



## ---- message=FALSE------------------------------------------------------

# Representation of the plot of one simulation as an example (the others are
# highly similar).
plot(simulation, show = "genotypes", type = "line",
     ylab = "Number of individuals", main = "Fully glycolytic tumours",
     font.main=2, font.lab=2, cex.main=1.4, cex.lab=1.1, las = 1)



## ---- message=FALSE------------------------------------------------------

# Definition of the function for creating the corresponding dataframe.
avc <- function (a, v, c) {
  data.frame(Genotype = c("WT", "GLY", "VOP", "DEF"),
             Fitness = c("1",
                         paste0("1 + ",a," * (f_1 + 1)"),
                         paste0("1 + ",a," * f_1 + ",v, " * (f_2 + 1) - ",c),
                         paste0("1 + ",a," * f_1 + ",v, " * f_2")
                         ))
                          }

# Specification of the different effects on fitness.
afavc <- allFitnessEffects(genotFitness = avc(2.5, 7, 1),
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

simulation <- oncoSimulIndiv(afavc,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 15,
                           mu = 1e-4,
                           initSize = 4000,
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)



## ---- message=FALSE------------------------------------------------------
## We get a huge number of VOP very quickly
## (too quickly?)
plot(simulation, show = "genotypes", type = "line",
     ylab = "Number of individuals", main = "Fully angiogenic tumours",
     font.main=2, font.lab=2, cex.main=1.4, cex.lab=1.1, las = 1)



## ---- message=FALSE------------------------------------------------------

# Definition of the function for creating the corresponding dataframe.
avc <- function (a, v, c) {
  data.frame(Genotype = c("WT", "GLY", "VOP", "DEF"),
             Fitness = c("1",
                         paste0("1 + ",a," * (f_1 + 1)"),
                         paste0("1 + ",a," * f_1 + ",v," * (f_2 + 1) - ",c),
                         paste0("1 + ",a," * f_1 + ",v," * f_2")
                         ))
                          }

# Specification of the different effects on fitness.
afavc <- allFitnessEffects(genotFitness = avc(7.5, 2, 1),
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

# Launching of the simulation (20 times).
simulation <- oncoSimulIndiv(afavc,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 25,
                           mu = 1e-4,
                           initSize = 4000,
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)



## ---- message=FALSE------------------------------------------------------

# Representation of the plot of one simulation as an example (the others are
# highly similar).
plot(simulation, show = "genotypes", type = "line",
     ylab = "Number of individuals", main = "Heterogeneous tumours",
     font.main=2, font.lab=2, cex.main=1.4, cex.lab=1.1, las = 1)



## ----example5_fitness, message=FALSE-------------------------------------
fitness_rel <- function(a, b, r, g, gt = c("WT", "S", "I", "D")) {
    data.frame(
      Genotype = gt,
      Fitness = c("1",
                  paste0("1 + ", a, " * f_3"),
                  paste0("1 + 1 - ", g),
                  paste0("1 + (1 - f_2 - f_3) * (1 - ", b, " + ", 
                         a, ") + f_2 * (1 - ", b, " + ", r,
                         ") + f_3 * (1 - 2 * ", b, ") + 1 - ", b, 
                         " + ", a, " + f_2 * (", r, " - ",
                         a, ") - f_3 * (", b, " + ", a, ")"))
                  )
}


## ----example5scen1, message=FALSE----------------------------------------
scen1 <- allFitnessEffects(genotFitness = fitness_rel(a = 0.5, b = 0.7,
                                                     r = 0.1, g = 0.8),
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

set.seed(1)
simulScen1 <- oncoSimulIndiv(scen1,
                              model = "McFL",
                              onlyCancer = FALSE,
                              finalTime = 70,
                              mu = 1e-4,
                              initSize = 5000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
op <- par(mfrow = c(1, 2))
plot(simulScen1, show = "genotypes", type = "line",
     main = "First scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
plot(simulScen1, show = "genotypes",
     main = "First scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
par(op)


## ----example5scen2, message=FALSE----------------------------------------
scen2 <- allFitnessEffects(genotFitness = fitness_rel(a = 0.3, b = 0.7,
                                                     r = 0.1, g = 0.7),
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

set.seed(1)

simulScen2 <- oncoSimulIndiv(scen2,
                              model = "McFL",
                              onlyCancer = FALSE,
                              finalTime = 70,
                              mu = 1e-4,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
op <- par(mfrow = c(1, 2))
plot(simulScen2, show = "genotypes", type = "line",
     main = "Second scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
plot(simulScen2, show = "genotypes",
     main = "Second scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
par(op)


## ----example5scen3, message=FALSE----------------------------------------
scen3 <- allFitnessEffects(genotFitness = fitness_rel(a = 0.2, b = 0.3,
                                                     r = 0.1, g = 0.3),
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

set.seed(1)
simulScen3 <- oncoSimulIndiv(scen3,
                              model = "McFL",
                              onlyCancer = FALSE,
                              finalTime = 50,
                              mu = 1e-4,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
op <- par(mfrow = c(1, 2))
plot(simulScen3, show = "genotypes", type = "line",
     main = "Third scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
plot(simulScen3, show = "genotypes",
     main = "Third scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
par(op)


## ----example5scen4, message=FALSE----------------------------------------
scen4 <- allFitnessEffects(genotFitness = fitness_rel(a = 0.2, b = 0.4,
                                                     r = 0.1, g = 0.3),
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

## Set a different seed to show the results better since
## with set.seed(1) the progression of I cells was not shown
set.seed(2)

simulScen4 <- oncoSimulIndiv(scen4,
                              model = "McFL",
                              onlyCancer = FALSE,
                              finalTime = 40,
                              mu = 1e-4,
                              initSize = 4000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)
op <- par(mfrow = c(1, 2))
plot(simulScen4, show = "genotypes", type = "line",
     main = "Fourth scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
plot(simulScen4, show = "genotypes",
     main = "Fourth scenario",
     cex.main = 1.4, cex.lab = 1.1,
     las = 1)
par(op)


## ------------------------------------------------------------------------
f_cells <- function(c1, c2, c3, r11, r12, r13, 
                    r21, r22, r23, r31, r32, r33, M, awt = 1e-4,
                                 gt = c("WT", "OC", "OB", "MM")) {
    data.frame(Genotype = gt,
               Fitness = c(
                  paste0("max(0.1, 1 - ", awt, " * (f_1+f_2+f_3)*N)"),
                  paste0("1", "+(((f_1*(", M, "-1)+1)*", c1, ")/", M, ")*",r11,
                          "+((((1-f_3)*(", M, "-1)-f_1*(", M, "-1)-1)*", c2, ")/", M, ")*", r12,
                          "+(((", M, "-(1-f_3)*(", M, "-1))*", c3, ")/", M, ")*", r13,
                          "-", c1
                          ),
                  paste0("1", "+(((f_2*(", M, "-1)+1)*", c2, ")/", M, ")*", r22,
                          "+((((1-f_1)*(", M, "-1)-f_2*(", M, "-1)-1)*", c3, ")/", M, ")*", r23,
                          "+(((", M, "-(1-f_1)*(", M, "-1))*", c1, ")/", M, ")*", r21,
                          "-", c2
                          ),
                  paste0("1", "+(((f_3*(", M, "-1)+1)*", c3, ")/", M, ")*", r33,
                          "+((((1-f_2)*(", M, "-1)-f_3*(", M, "-1)-1)*", c1, ")/", M, ")*", r31,
                          "+(((", M, "-(1-f_2)*(", M, "-1))*", c2, ")/", M, ")*", r32,
                          "-", c3
                          )
                  )
               ,stringsAsFactors = FALSE
               )
}



## ------------------------------------------------------------------------
N <- 40000
M <- 10
c1 <- 1
c2 <- 1.2
c3 <- 1.4
r11 <- 0
r12 <- 1
r13 <- 2.5
r21 <- 1
r22 <- 0
r23 <- -0.3
r31 <- 2.5
r32 <- 0
r33 <- 0

fe_cells <-
    allFitnessEffects(
        genotFitness =
            f_cells(c1, c2, c3, r11, r12, r13, 
                    r21, r22, r23, r31, r32, r33, M,
                    gt = c("WT", "OC", "OB", "MM")),
        frequencyDependentFitness = TRUE,
        frequencyType = "rel")

## Simulated trajectories
set.seed(2)
simulation <- oncoSimulIndiv(fe_cells,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 30,
                           mu = c("OC"=1e-1, "OB"=1e-1, "MM"=1e-4),
                           initSize = N,
                           keepPhylog = FALSE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)

                                        #Plot trajectorie
op <- par(mfrow = c(1, 2))
plot(simulation, show = "genotypes", type = "line")
plot(simulation, show = "genotypes")
par(op)



## ------------------------------------------------------------------------
N <- 40000
M <- 10
c1 <- 1
c2 <- 1
c3 <- 1
r11 <- 0
r12 <- 1
r13 <- 0.5
r21 <- 1
r22 <- 0
r23 <- -0.3
r31 <- 0.5
r32 <- 0
r33 <- 0

fe_cells <-
    allFitnessEffects(
        genotFitness =
            f_cells(c1, c2, c3, r11, r12, r13, 
                    r21, r22, r23, r31, r32, r33, M,
                    gt = c("WT", "OC", "OB", "MM")),
        frequencyDependentFitness = TRUE,
        frequencyType = "rel")

## Simulated trajectories
set.seed(1)
simulation <- oncoSimulIndiv(fe_cells,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 25,
                           mu = c("OC"=1e-1, "OB"=1e-1, "MM"=1e-4),
                           initSize = N,
                           keepPhylog = FALSE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)

#Plot trajectorie
plot(simulation, show = "genotypes")



## ----parka , message=FALSE-----------------------------------------------
park1<- data.frame(Genotype = c("WT", "A", "B", "A,B"), 
                 Fitness = c("1",
			     "1 + 3*(f_1 + f_2 + f_1_2)",
                 "1 + 2*(f_1 + f_2 + f_1_2)", ## We establish
                 ## the fitness of B smaller than the one of A because
                 ## it is an indirect cause of the disease and not a direct one. 
                 "1.5 + 4.5*(f_1 + f_2 + f_1_2)")) ## The baseline
                  ## of the fitness is higher in the
                  ## AB population (their growth is favored).

parkgen1<- allFitnessEffects(genotFitness = park1, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel") 


## ----parkb---------------------------------------------------------------
set.seed(1)
fpark1 <- oncoSimulIndiv(parkgen1, 
                     model = "McFL",
                     onlyCancer = FALSE,
                     finalTime = 100,
                     mu = 1e-4,
                     initSize = 5000,
                     keepPhylog = TRUE,
                     seed = NULL,
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)




## ----parkplot1, message=FALSE--------------------------------------------
plot(fpark1, show = "genotypes",  type = "line",
     col = c("black", "green", "red", "blue"))


## ----parkplot2, message=FALSE--------------------------------------------
plotClonePhylog(fpark1, N = 0, keepEvents=TRUE, timeEvents=TRUE)


## ----wuAMicrobesfit1-----------------------------------------------------
create_fe <- function(bG, cG, iPA, cI, cS, bPA, ab,
                      gt = c("WT", "CT", "PA")) {
  data.frame(Genotype = gt,
             Fitness = c(
               paste0("1 + ", bG, " * (f_ + f_1) - ", cS, 
                      " * (f_ +  f_1 + f_2) - ", cI, "(f_2 > 0.2) - ", cG, 
                      " - ", ab),
               paste0("1 + ", bG, " * (f_ + f_1) - ", cS, 
                      " * (f_ +  f_1 + f_2) - ", cG),
               paste0("1 +", bPA, " - ", cS, " * (f_ +  f_1 + f_2) - ", iPA,
                      " *(f_(f_2 > 0.2))")),
             stringsAsFactors = FALSE)
}



## ----wuAMicrobescheck----------------------------------------------------
create_fe("bG", "cG", "iPA", "cI", "cS", "bPA", "ab")


## ----wuAMicrobes2a-------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(7, 1, 9, 0.5, 2, 5, 0), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 1000, 1000)))



## ----wuAMicrobes2b-------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(7, 1, 9, 0.5, 2, 5, 2), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 1000, 1000)))



## ----wuAMicrobes2c-------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(7, 1, 9, 0.5, 2, 5, 2), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(0, 0, 1000)))



## ----wuAMicrobes2d-------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(7, 1, 9, 0.5, 2, 5, 2), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(100, 0, 1000)))



## ----wuAMicrobes2e-------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(7, 1, 9, 0.5, 2, 5, 5), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 0, 0)))



## ----wuAMicrobes2f-------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
				   create_fe(7, 1, 9, 0.5, 2, 5, 5), 
				   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 0, 1000)))



## ----woAntib1, eval = FALSE----------------------------------------------
## woAntib <- allFitnessEffects(
##   genotFitness = create_fe(7, 1, 9, 0.5, 2, 5, 0),
##   frequencyDependentFitness = TRUE,
##   frequencyType = "rel")
## 
## ## We do not run this for speed but load it below
## set.seed(2)
## woAntibS <- oncoSimulIndiv(woAntib,
##                         model = "McFL",
##                         onlyCancer = FALSE,
##                         finalTime = 2000,
##                         mu = 1e-4,
##                         initSize = 1000,
##                         keepPhylog = FALSE,
##                         seed = NULL,
##                         errorHitMaxTries = FALSE,
##                         errorHitWallTime = FALSE,
##                         keepEvery = 2 ## store a smaller object
##                         )
## 


## ----woAntib1data--------------------------------------------------------
## Load stored results
data(woAntibS)


## ----woAntib2------------------------------------------------------------

plot(woAntibS, show = "genotypes", type = "line",
     col = c("black", "green", "red"))



## ----wiAntib1------------------------------------------------------------
wiAntib <- allFitnessEffects(
  genotFitness = create_fe(7, 1, 9, 0.5, 2, 5, 2),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel")


set.seed(2)
wiAntibS <- oncoSimulIndiv(wiAntib,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 100,
                        mu = 1e-4,
                        initSize = 1000, 
                        keepPhylog = TRUE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)

plot(wiAntibS, show = "genotypes", type = "line",
     col = c("black", "green", "red"))



## ----breastCfit1---------------------------------------------------------
create_fe <- function(cG, bG, cS, cMMph, cMMTC, bR, cD,
                      gt = c("WT", "Mph", "BTC", "MTC")) {
  data.frame(Genotype = gt,
             Fitness = c(
               paste0("1 + ", bG, "*(f_ + f_1) - ", cG, " - ", cS, "*(f_ + f_2)"),
               paste0("1 + ", bG, "*(f_ + f_1) - ", cG, " - ", cMMph),
               paste0("1 + ", bR, " + ", bG, "*(f_ + f_1) - ", cS, "* (f_ + f_2) -",
                      cD , " * f_1"),
               paste0("1 + ", bR, " + ", bG, " *(f_ + f_1) -", cMMTC, " - ",
                      cD , " * f_1")
               ),
             stringsAsFactors = FALSE)
}



## ----breastCcheck--------------------------------------------------------
create_fe("cG", "bG","cS", "cMMph", "cMMTC", "bR", "cD")


## ----breastC2a-----------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(2, 5, 1, 0.8, 1, 1, 9), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 0, 0, 0)))



## ----breastC2b-----------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(2, 5, 1, 0.8, 1, 1, 9), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 1000, 0, 0)))



## ----breastC2c-----------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(2, 5, 1, 0.8, 1, 1, 9), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 1000, 100, 100)))



## ----cancercontrol1------------------------------------------------------
afe_3_a <- allFitnessEffects(
  genotFitness =
    create_fe(0.5, 4, 1, 0.2, 1, 0.5, 4),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel")

set.seed(2)

s_3_a <- oncoSimulIndiv(afe_3_a,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 50,
                        mu = 1e-4,
                        initSize = 10000, 
                        keepPhylog = TRUE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)

plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow"))



## ----cancerNM1-----------------------------------------------------------
afe_3_a <- allFitnessEffects(
  genotFitness =
    create_fe(1, 4, 0.5, 1, 1.5, 1, 4),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel")

set.seed(2)

s_3_a <- oncoSimulIndiv(afe_3_a,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 50,
                        mu = 1e-4,
                        initSize = 10000, 
                        keepPhylog = TRUE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)


plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow"))



## ----cancerM1------------------------------------------------------------

afe_3_a <- allFitnessEffects(
  genotFitness =
    create_fe(0.5, 4, 2, 0.5, 0.5, 1, 4),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel")

set.seed(2)

s_3_a <- oncoSimulIndiv(afe_3_a,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 50,
                        mu = 1e-4,
                        initSize = 10000, 
                        keepPhylog = TRUE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)


plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow"))




## ----breastCQfit1--------------------------------------------------------
create_fe <- function(cG, bG, cS, cMMph, cMMTC, bR, cD, Q,
                      gt = c("WT", "BTC", "R", "MTC", "Mph", "BTC,R", "MTC,R")) {
  data.frame(Genotype = gt,
             Fitness = c(
               
               paste0("1 + ", bG, "(f_ + f_2 + f_4) - ", cG, " - ", cS, "(f_ + f_1 + f_2 + f_1_2) -",
                      "0.01*", Q),
               paste0("1 + ", bR, " + ", bG, "(f_ + f_2 + f_4) - ", cS, " (f_ + f_1 + f_2 + f_1_2) -",
                      cD , " * f_4 -", Q),
               paste0("1 + ", bG, "(f_ + f_2 + f_4) - ", cG, " - ", cS, "(f_ + f_1 + f_2 + f_1_2)"),
               paste0("1 + ", bR, " + ", bG, " *(f_ + f_2 + f_4) -", cMMTC, " - ",
                      cD , " * f_4 -", Q),
               paste0("1 + ", bG, "(f_ + f_2 + f_4) - ", cG, " - ", cMMph, "- 0.01*",Q),
               paste0("1 + ", bR, " + ", bG, "(f_ + f_2 + f_4) - ", cS, " (f_ + f_1 + f_2 + f_1_2) -",
                      cD , " * f_4"),
               paste0("1 + ", bR, " + ", bG, " *(f_ + f_2 + f_4) -", cMMTC, " - ",
                      cD , " * f_4")
             ),
             stringsAsFactors = FALSE)
}



## ----breastCQcheck-------------------------------------------------------
create_fe("cG", "bG","cS", "cMMph", "cMMTC", "bR", "cD", "Q")


## ----breastCQ2a----------------------------------------------------------
evalAllGenotypes(allFitnessEffects(genotFitness =  
                                   create_fe(2,5,1,0.8,1,1,9,5), 
                                   frequencyDependentFitness = TRUE,
                                   frequencyType = "rel",
                                   spPopSizes = c(1000, 100, 0, 100, 1000, 0, 0)))



## ----cancerwoQ1----------------------------------------------------------

afe_3_a <- allFitnessEffects(
  genotFitness =
    create_fe(2, 5, 1, 0.8, 1, 1, 9, 0),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel")

#Set mutation rates
muvar2 <- c("Mph" = 1e-2, "BTC" = 1e-3, "MTC"=1e-3, "R" = 1e-7)

set.seed(2)
s_3_a <- oncoSimulIndiv(afe_3_a,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 20,
                        mu = muvar2,
                        initSize = 10000, 
                        keepPhylog = TRUE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)


plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow", "orange", "brown"))



## ----cancerwQlRm1--------------------------------------------------------

afe_3_a <- allFitnessEffects(
  genotFitness =
    create_fe(2, 5, 1, 0.8, 1, 1, 9, 2),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel")

muvar2 <- c("Mph" = 1e-2, "BTC" = 1e-3, "MTC"=1e-3, "R" = 1e-7)

set.seed(2)

s_3_a <- oncoSimulIndiv(afe_3_a,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 20,
                        mu = muvar2,
                        initSize = 10000, 
                        keepPhylog = TRUE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)


plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow", "orange", "brown"))



## ----cancerwQHRm1--------------------------------------------------------

afe_3_a <- allFitnessEffects(
  genotFitness =
    create_fe(2, 5, 1, 0.8, 1, 1, 9, 2),
  frequencyDependentFitness = TRUE,
  frequencyType = "rel")

muvar2 <- c("Mph" = 1e-2, "BTC" = 1e-3, "MTC"=1e-3, "R" = 1e-5)

set.seed(2)

s_3_a <- oncoSimulIndiv(afe_3_a,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 50,
                        mu = muvar2,
                        initSize = 10000, 
                        keepPhylog = TRUE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)


plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow", "orange", "brown"))



## ----lod_pom_ex----------------------------------------------------------
pancr <- allFitnessEffects(
    data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                          "TP53", "TP53", "MLL3"),
               child = c("KRAS","SMAD4", "CDNK2A", 
                         "TP53", "MLL3",
                         rep("PXDN", 3), rep("TGFBR2", 2)),
               s = 0.05, sh = -0.3, typeDep = "MN"))

pancr16 <- oncoSimulPop(16, pancr, model = "Exp",
                        mc.cores = 2)

## Look a the first POM 
str(POM(pancr16)[1:3])

LOD(pancr16)[1:2]

## The diversity of LOD (lod_single) and POM might or might not
## be identical
diversityPOM(POM(pancr16))
diversityLOD(LOD(pancr16))

## Show the genotypes and their diversity (which might, or might
## not, differ from the diversity of LOD and POM)
sampledGenotypes(samplePop(pancr16))



## ------------------------------------------------------------------------
## No seed fixed, so reruns will give different DAGs.
(a1 <- simOGraph(10))
library(graph) ## for simple plotting
plot(as(a1, "graphNEL"))


## ----simographindirect, eval=FALSE,echo=TRUE-----------------------------
## g2 <- simOGraph(4, out = "rT", removeDirectIndirect = FALSE)
## 
## fe_from_d <- allFitnessEffects(g2)
## fitness_d <- evalAllGenotypes(fe_from_d)
## 
## fe_from_t <- allFitnessEffects(genotFitness =
##                           OncoSimulR:::allGenotypes_to_matrix(fitness_d))
## 						
## ## Compare
## fitness_d
## (fitness_t <- evalAllGenotypes(fe_from_t))
## 
## identical(fitness_d, fitness_t)
## 						
## 
## ## ... but to be safe use fe_from_t as the fitnessEffects object for simulations
## 


## ------------------------------------------------------------------------
## This code will only be evaluated under Windows
if(.Platform$OS.type == "windows")
    try(pancrError <- oncoSimulPop(10, pancr,
                               initSize = 1e-5,
                               detectionSize = 1e7,
                               keepEvery = 10,
                               mc.cores = 2))


## ------------------------------------------------------------------------
## Do not run under Windows
if(.Platform$OS.type != "windows")
    pancrError <- oncoSimulPop(10, pancr,
                               initSize = 1e-5,
                               detectionSize = 1e7,
                               keepEvery = 10,
                               mc.cores = 2)


## ---- eval=FALSE---------------------------------------------------------
## pancrError[[1]]


## ---- eval=FALSE---------------------------------------------------------
## pancrError[[1]][1]


## ----ex-tomlin1exc-------------------------------------------------------
sd <- 0.1 ## fitness effect of drivers
sm <- 0 ## fitness effect of mutator
nd <- 20 ## number of drivers
nm <- 5  ## number of mutators
mut <- 50 ## mutator effect  THIS IS THE DIFFERENCE

fitnessGenesVector <- c(rep(sd, nd), rep(sm, nm))
names(fitnessGenesVector) <- 1:(nd + nm)
mutatorGenesVector <- rep(mut, nm)
names(mutatorGenesVector) <- (nd + 1):(nd + nm)

ft <- allFitnessEffects(noIntGenes = fitnessGenesVector,
                        drvNames = 1:nd)
mt <- allMutatorEffects(noIntGenes = mutatorGenesVector)



## ----ex-tomlinexc2-------------------------------------------------------
ddr <- 4
set.seed(2)
RNGkind("L'Ecuyer-CMRG")
st <- oncoSimulPop(4, ft, muEF = mt,
                   detectionDrivers = ddr,
                   finalTime = NA,
                   detectionSize = NA,
                   detectionProb = NA,
                   onlyCancer = TRUE,
                   keepEvery = NA, 
                   mc.cores = 2, ## adapt to your hardware
                   seed = NULL) ## for reproducibility
set.seed(NULL) ## return things to their "usual state"



## ---- fig.height=3-------------------------------------------------------
## Node 2 and 3 depend on 1, and 4 depends on no one
p1 <- cbind(c(1L, 1L, 0L), c(2L, 3L, 4L))
plotPoset(p1, addroot = TRUE)


## ---- fig.height=3-------------------------------------------------------
## A simple way to create a poset where no gene (in a set of 15)
## depends on any other.
p4 <- cbind(0L, 15L)
plotPoset(p4, addroot = TRUE)


## ---- fig.height=3-------------------------------------------------------
pancreaticCancerPoset <- cbind(c(1, 1, 1, 1, 2, 3, 4, 4, 5),
                               c(2, 3, 4, 5, 6, 6, 6, 7, 7))
storage.mode(pancreaticCancerPoset) <- "integer"
plotPoset(pancreaticCancerPoset,
          names = c("KRAS", "SMAD4", "CDNK2A", "TP53",
                    "MLL3","PXDN", "TGFBR2"))



## ---- echo=FALSE,results='hide',error=FALSE------------------------------
options(width=60)


## ------------------------------------------------------------------------
## use poset p1101
data(examplePosets)
p1101 <- examplePosets[["p1101"]]

## Bozic Model
b1 <- oncoSimulIndiv(p1101, keepEvery = 15)
summary(b1)


## ----pb2bothx1,fig.height=5.5, fig.width=5.5-----------------------------
b2 <- oncoSimulIndiv(p1101, keepEvery = 1)
summary(b2)
plot(b2)


## ----pbssttx1,eval=FALSE-------------------------------------------------
## plot(b2, type = "stacked")


## ---- echo=FALSE,eval=TRUE-----------------------------------------------
set.seed(1) ## for reproducibility. Once I saw it not reach cancer.

## ------------------------------------------------------------------------

m2 <- oncoSimulIndiv(examplePosets[["p1101"]], model = "McFL", 
                     numPassengers = 0, detectionDrivers = 8, 
                     mu = 5e-7, initSize = 4000, 
                     sampleEvery = 0.025,
                     finalTime = 25000, keepEvery = 5, 
                     detectionSize = 1e6) 


## ----m2x1,fig.width=6.5, fig.height=10-----------------------------------
par(mfrow = c(2, 1))
plot(m2, addtot = TRUE, log = "",
     thinData = TRUE, thinData.keep = 0.5)
plot(m2, type = "stacked",
     thinData = TRUE, thinData.keep = 0.5)


## ------------------------------------------------------------------------
b3 <- oncoSimulIndiv(p1101, onlyCancer = FALSE)
summary(b3)

b4 <- oncoSimulIndiv(p1101, onlyCancer = FALSE)
summary(b4)


## ----b3b4x1ch1, fig.width=8, fig.height=4--------------------------------
par(mfrow = c(1, 2))
par(cex = 0.8) ## smaller font
plot(b3)
plot(b4)


## ----ch2-----------------------------------------------------------------
p1 <- oncoSimulPop(4, p1101, mc.cores = 2)
par(mfrow = c(2, 2))
plot(p1, ask = FALSE)


## ----p1multx1,eval=FALSE-------------------------------------------------
## par(mfrow = c(2, 2))
## plot(p1, type = "stream", ask = FALSE)


## ----p1multstx1,eval=FALSE-----------------------------------------------
## par(mfrow = c(2, 2))
## plot(p1, type = "stacked", ask = FALSE)


## ------------------------------------------------------------------------

m1 <- oncoSimulPop(100, examplePosets[["p1101"]], 
                   numPassengers = 0, mc.cores = 2)



## ------------------------------------------------------------------------
genotypes <- samplePop(m1)


## ----fxot1,fig.width=4, fig.height=4-------------------------------------
colSums(genotypes)/nrow(genotypes)

require(Oncotree)
ot1 <- oncotree.fit(genotypes)
plot(ot1)


## ----fxot2,fig.width=4, fig.height=4-------------------------------------
genotypesSC <- samplePop(m1, typeSample = "single")
colSums(genotypesSC)/nrow(genotypesSC)

ot2 <- oncotree.fit(genotypesSC)
plot(ot2)


## ------------------------------------------------------------------------
sessionInfo()


## ---- echo=FALSE, eval=TRUE----------------------------------------------
## reinitialize the seed
set.seed(NULL)

