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

test_that("find number in muExpression for oncoSimulIndiv", {
  
  Rcpp::sourceCpp(file = "OncoSimul/miscell-files/test-findNumber.cpp")
  muexpression = "if(T>300) 10000; else 1;"
  expect_identical(readr::parse_number(muexpression), findNumber(muexpression))
  
})

test_that("find number when used relative or absolute frequencies", {
  
  Rcpp::sourceCpp(file = "OncoSimul/miscell-files/test-findNumberFvars.cpp")
  muexpressionRel = "if(f_>0.3) 100; else 1;"
  muexpressionAbs = "if(n_1>10) 100; else 1;"
  expect_identical("f_", findVars(muexpressionRel))
  expect_identical("n_1", findVars(muexpressionAbs))
})

test_that("print oncosimul object when fvars", {
  
  fl <- data.frame(
    Genotype = c("WT", "A", "B"),
    Fitness = c("1 + f_1*f_2",
                "1 + 1.5*f_*f_2",
                "1 + 1.5*f_*f_1"),
    stringsAsFactors = FALSE
  )
  
  fe <- allFitnessEffects(genotFitness = fl, frequencyDependentFitness = TRUE)
  
  muexpression = "if(f_>0.3) 100; else 1;"
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

cat(paste("\n Ending test.oncoSimulIndiv-MutationRate at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)