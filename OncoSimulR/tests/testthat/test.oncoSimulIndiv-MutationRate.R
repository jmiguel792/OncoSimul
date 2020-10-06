inittime <- Sys.time()
cat(paste("\n Starting oncoSimulIndiv-MutationRate tests", date(), "\n"))


test_that("testing output classes when FDF", {
  
  fl <- data.frame(Genotype = c("WT", "A", "B"),
                   Fitness = c("1 + f_1*f_2",
                               "1 + 1.5*f_*f_2",
                               "1 + 1.5*f_*f_1"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = fl, frequencyDependentFitness = TRUE)
  
  expect_message(allFitnessEffects(genotFitness = fl, frequencyDependentFitness = TRUE),
                 "frequencyType set to 'auto'")
  
  set.seed(2)
  
  muExpr1 <- 100 
  muExpr2 <- "if(T>200) 10000; else 1;"
  
  expect_error(out1 <- oncoSimulIndiv(fe,
                                      model = "McFL", 
                                      onlyCancer = FALSE, 
                                      finalTime = 500,
                                      mu = 1e-6,
                                      muFactor = muExpr1,
                                      initSize = 5000, 
                                      keepPhylog = FALSE,
                                      seed = NULL, 
                                      errorHitMaxTries = FALSE, 
                                      errorHitWallTime = FALSE),
               "muFactor must be a string")
  
  out2 <- oncoSimulIndiv(fe,
                        model = "McFL", 
                        onlyCancer = FALSE, 
                        finalTime = 500,
                        mu = 1e-6,
                        muFactor = muExpr2,
                        initSize = 5000, 
                        keepPhylog = FALSE,
                        seed = NULL, 
                        errorHitMaxTries = FALSE, 
                        errorHitWallTime = FALSE)
  
  expect_identical(class(fl), "data.frame")
  expect_identical(class(fe), "fitnessEffects")
  expect_identical(class(out2), c("oncosimul", "oncosimul2"))
  expect_output(print(out2), "Individual OncoSimul trajectory", fixed = TRUE)
  
  expect_message(out2 <- oncoSimulIndiv(fe,
                                        model = "McFL", 
                                        onlyCancer = FALSE, 
                                        finalTime = 500,
                                        mu = 1e-6,
                                        muFactor = muExpr2,
                                        initSize = 5000, 
                                        keepPhylog = FALSE,
                                        seed = NULL, 
                                        errorHitMaxTries = FALSE, 
                                        errorHitWallTime = FALSE), 
                 "Exprtk expression for mutation rate defined.")
  
  expect_identical(out2$NumClones, 3)
  expect_true(out2$TotalPopSize > 5000)
  expect_true(out2$geneNames[2] %in% LETTERS)
  expect_false(out2$FinalTime > 500)
  
})

test_that("testing output when FDF", {
  
  fl <- data.frame(Genotype = c("WT", "A", "B"),
                   Fitness = c("1 + f_1*f_2",
                               "1 + 1.5*f_*f_2",
                               "1 + 1.5*f_*f_1"),
                   stringsAsFactors = FALSE)
  
  fe <- allFitnessEffects(genotFitness = fl, frequencyDependentFitness = TRUE)
  
  set.seed(2)
  sim1 <- oncoSimulIndiv(fe,
                         model = "McFL",
                         onlyCancer = FALSE, 
                         finalTime = 500,
                         mu = 1e-6,
                         initSize = 5000, 
                         keepPhylog = FALSE,
                         seed = NULL, 
                         errorHitMaxTries = FALSE, 
                         errorHitWallTime = FALSE)
  
  muExpr <- "if(T>200) 10000; else 1;"
  
  set.seed(2)
  sim3 <- oncoSimulIndiv(fe,
                         model = "McFL", 
                         onlyCancer = FALSE, 
                         finalTime = 500,
                         mu = 1e-6,
                         muFactor = muExpr,
                         initSize = 5000, 
                         keepPhylog = FALSE,
                         seed = NULL, 
                         errorHitMaxTries = FALSE, 
                         errorHitWallTime = FALSE)
  
  wt_sim1 <- sim1$pops.by.time[nrow(sim1$pops.by.time), ][2]
  wt_sim3 <- sim3$pops.by.time[nrow(sim3$pops.by.time), ][2]
  
  A_sim1 <- sim1$pops.by.time[nrow(sim1$pops.by.time), ][3]
  A_sim3 <- sim3$pops.by.time[nrow(sim3$pops.by.time), ][3]
    
  B_sim1 <- sim1$pops.by.time[nrow(sim1$pops.by.time), ][4]
  B_sim3 <- sim3$pops.by.time[nrow(sim3$pops.by.time), ][4]
  
  expect_true(wt_sim1 > wt_sim3)
  expect_true(A_sim1 < A_sim3)
  expect_true(B_sim1 < B_sim3)
  
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
  
  muexpression1 = "if(f_ > 0.3 or f_1 > 0.3 or f_2 > 0.3) 100; else 1;"
  
  sim <- oncoSimulIndiv(fe,
                         model = "McFL", 
                         onlyCancer = FALSE, 
                         finalTime = 500,
                         mu = 1e-6,
                         muFactor = muexpression1,
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
  
  muexpression2 = "if(n_ > 10 or n_1 > 10 or n_2 > 10) 100; else 1;"
  tmp <- oncoSimulIndiv(afear3, 
                         model = "McFL", 
                         onlyCancer = FALSE, 
                         finalTime = 100,
                         mu = 1e-4,
                         muFactor = muexpression2,
                         initSize = 5000, 
                         keepPhylog = FALSE,
                         seed = NULL, 
                         errorHitMaxTries = FALSE, 
                         errorHitWallTime = FALSE)
  
  expect_output(print(tmp), "Individual OncoSimul trajectory with call")
  
})

test_that("print oncosimul object when fdf is not used", {
  
  pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A",
                                                   "TP53", "TP53", "MLL3"),
                                        child = c("KRAS","SMAD4", "CDNK2A",
                                                  "TP53", "MLL3",
                                                  rep("PXDN", 3), rep("TGFBR2", 2)),
                                        s = 0.05,
                                        sh = -0.3,
                                        typeDep = "MN"))
  
  muExpr1 <- "if(T>4000) (1/sqrt(T) + 1/N)*100; else 1;"
  ep2 <- oncoSimulIndiv(pancr, model = "McFL",
                        mu = 1e-6,
                        muFactor = muExpr1,
                        sampleEvery = 0.02,
                        keepEvery = 1,
                        initSize = 1000,
                        finalTime = 10000,
                        onlyCancer = FALSE)
  
  expect_output(print(ep2), "Individual OncoSimul trajectory with call")
  
  muExpr2 <- "if(T>200) log(T); else 1 + 1/N;"
  pancr2 <- oncoSimulIndiv(pancr, model = "Exp", muFactor = muExpr2)
  
  expect_output(print(pancr2), "Individual OncoSimul trajectory with call")
  
})

cat(paste("\n Ending test.oncoSimulIndiv-MutationRate at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)