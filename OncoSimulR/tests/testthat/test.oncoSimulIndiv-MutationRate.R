inittime <- Sys.time()
cat(paste("\n Starting oncoSimulIndiv-MutationRate tests", date(), "\n"))

test_that("printing oncosimul object", {
  
  fl <- data.frame(Genotype = c("WT", "A", "B"),
                   Fitness = c("1 + f_1*f_2",
                               "1 + 1.5*f_*f_2",
                               "1 + 1.5*f_*f_1"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = fl, frequencyDependentFitness = TRUE)
  
  muexpression = "if(T>300) 10000; else 1;"
  out <- oncoSimulIndiv(fe,
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
  
  expect_output(print(out), "Individual OncoSimul trajectory with call")
  
})

#test_that("find number in muExpression for oncoSimulIndiv", {
#  
# 
#  Rcpp::sourceCpp(file = "/home/jmiguel/OncoSimul/miscell-files/test-findNumber.cpp")
#  muexpression = "if(T>300) 10000; else 1;"
#  expect_identical(readr::parse_number(muexpression), findNumber(muexpression))
#  
#})

#test_that("find number when used relative or absolute frequencies", {
#  
#  Rcpp::sourceCpp("/home/jmiguel/OncoSimul/miscell-files/test-findNumberFvars.cpp")
#  muexpressionRel = "if(f_>0.3) 100; else 1;"
#  muexpressionAbs = "if(n_1>10) 100; else 1;"
#  expect_identical("f_", findVars(muexpressionRel))
#  expect_identical("n_1", findVars(muexpressionAbs))
#})

test_that("print oncosimul object when fvars", {
  
  fl <- data.frame(
    Genotype = c("WT", "A", "B"),
    Fitness = c("1 + f_1*f_2",
                "1 + 1.5*f_*f_2",
                "1 + 1.5*f_*f_1"),
    stringsAsFactors = FALSE
  )
  
  fe <- allFitnessEffects(genotFitness = fl, frequencyDependentFitness = TRUE)
  
  muexpression = "if(f_ > 0.3) 100; else 1;"
  sim <- oncoSimulIndiv(fe,
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
  
  expect_output(print(sim), "Individual OncoSimul trajectory with call")
  
  rar3 <- data.frame(Genotype = c("WT", "A", "B", "C"), 
                     Fitness = c("1",
                                 "1.1 + .3*(n_2/N)",
                                 "1.2 + .4*(n_1/N)",
                                 "1.0 + .5 * ( n_1 > 20)"),
                     stringsAsFactors = FALSE)
  
  afear3 <- allFitnessEffects(genotFitness = rar3, frequencyDependentFitness = TRUE)
  
  muexpression = "if(n_1 > 10) 100; else 1;"
  tmp <- oncoSimulIndiv(afear3, 
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
  
  expect_output(print(tmp), "Individual OncoSimul trajectory with call")
  
})

test_that("print oncosimul object when fdf is not used: mutators", {
  
  set.seed(1) ## for reproducibility -> MUTATORS
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
  
  expect_output(print(mue1), "Individual OncoSimul trajectory with call")
  
})

test_that("print oncosimul object when fdf is not used: no-interactions genes", {
  
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
  
  expect_output(print(ep1), "Individual OncoSimul trajectory with call")

})

cat(paste("\n Ending test.oncoSimulIndiv-MutationRate at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)