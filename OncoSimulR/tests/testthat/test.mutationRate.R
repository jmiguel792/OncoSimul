## Test for mutator rate functionality: updating if needed by exprtk
inittime <- Sys.time()
cat(paste("\n Starting test.mutatorRate.R test at", date()))
date()

test_that("eval mutator rate functionality with FDF situation and mu size of 1", {
  
  fl <- data.frame(
    Genotype = c("WT", "A", "B", "A, B"),
    Fitness = c("1 + f_1*f_2",
                "1 + 1.5*f_*f_2",
                "1 + 1.5*f_*f_1",
                "1 + 1.5*f_1*f_2"),
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
  
  muexpression <- "if(T>=0) 100; else 1;"
  set.seed(2)
  expect_message(sim2 <- oncoSimulIndiv(fe,
                         model = "McFL", 
                         onlyCancer = FALSE, 
                         finalTime = 500,
                         mu = 1e-6,
                         muFactor = muexpression,
                         initSize = 5000, 
                         keepPhylog = FALSE,
                         seed = NULL, 
                         errorHitMaxTries = FALSE, 
                         errorHitWallTime = FALSE),
                 "Exprtk expression for mutation rate defined.")
  
  expect_equal(
    sim1$pops.by.time[nrow(sim1$pops.by.time), ],
               sim2$pops.by.time[nrow(sim2$pops.by.time), ]
    )
})


test_that("expecting errors with mutation rate funcionality", {
  
  fl <- data.frame(
    Genotype = c("WT", "A", "B", "A, B"),
    Fitness = c("1 + f_1*f_2",
                "1 + 1.5*f_*f_2",
                "1 + 1.5*f_*f_1",
                "1 + 1.5*f_1*f_2"),
    stringsAsFactors = FALSE
  )
  
  fe <- allFitnessEffects(genotFitness = fl, 
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")
  
  muexpression1 <- 10
  set.seed(2)
  expect_error(sim2 <- oncoSimulIndiv(fe,
                                        model = "McFL", 
                                        onlyCancer = FALSE, 
                                        finalTime = 500,
                                        mu = 1e-6,
                                        muFactor = muexpression1,
                                        initSize = 5000, 
                                        keepPhylog = FALSE,
                                        seed = NULL, 
                                        errorHitMaxTries = FALSE, 
                                        errorHitWallTime = FALSE),
                 "muFactor must be a string")
})

test_that("eval mutator rate functionality with FDF situation and mu size > 1", {
  
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
  
  afe_3_a <- allFitnessEffects(
    genotFitness =
      create_fe(2, 5, 1, 0.8, 1, 1, 9, 0),
    frequencyDependentFitness = TRUE,
    frequencyType = "rel")
  
  muvar2 <- c("Mph" = 1e-4, "BTC" = 1e-5, "MTC"=1e-5, "R" = 1e-6)
  
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
  
  muvar3 <- c("Mph" = 1e-5, "BTC" = 1e-6, "MTC"=1e-6, "R" = 1e-7)
  muexpression <- "if(T>=0) 10; else 1;"
  
  set.seed(2)
  s_3_b <- oncoSimulIndiv(afe_3_a,
                          model = "McFL", 
                          onlyCancer = FALSE, 
                          finalTime = 20,
                          mu = muvar3,
                          muFactor = muexpression,
                          initSize = 10000, 
                          keepPhylog = TRUE,
                          seed = NULL, 
                          errorHitMaxTries = FALSE, 
                          errorHitWallTime = FALSE)
  
  expect_equal(
    s_3_a$pops.by.time[nrow(s_3_a$pops.by.time), ],
    s_3_b$pops.by.time[nrow(s_3_b$pops.by.time), ]
  )
  
})

test_that("eval mutator rate functionality with No-FDF situation and mu size == 1", {
  
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
  set.seed(5)
  ep1 <- oncoSimulIndiv(sv2, model = "McFL",
                        mu = 5e-4,
                        muFactor = "None",
                        sampleEvery = 0.025,
                        keepEvery = 0.5,
                        initSize = 2000,
                        finalTime = 3000,
                        onlyCancer = FALSE)
  
  RNGkind("Mersenne-Twister")
  set.seed(5)
  ep2 <- oncoSimulIndiv(sv2, model = "McFL",
                        mu = 5e-5,
                        muFactor = "if(T>=0) 10; else 1;",
                        sampleEvery = 0.025,
                        keepEvery = 0.5,
                        initSize = 2000,
                        finalTime = 3000,
                        onlyCancer = FALSE)
  
  expect_equal(
    ep1$pops.by.time[nrow(ep1$pops.by.time), ],
    ep2$pops.by.time[nrow(ep2$pops.by.time), ]
  )
  
})

test_that("test funcionality at T>X which only happens once", {
  
  fl <- data.frame(
    Genotype = c("WT", "A", "B", "A, B"),
    Fitness = c("1 + f_1*f_2",
                "1 + 1.5*f_*f_2",
                "1 + 1.5*f_*f_1",
                "1 + 1.5*f_1*f_2"),
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
  
  muexpression = "if(T>200) 100; else 1;"
  set.seed(2)
  expect_message(sim2 <- oncoSimulIndiv(fe,
                         model = "McFL", 
                         onlyCancer = FALSE, 
                         finalTime = 1000,
                         mu = 1e-6,
                         muFactor = muexpression,
                         initSize = 5000, 
                         keepPhylog = FALSE,
                         seed = NULL, 
                         errorHitMaxTries = FALSE, 
                         errorHitWallTime = FALSE),
                 "Exprtk expression for mutation rate defined.")
  
  expect_true(
    sim1$pops.by.time[nrow(sim1$pops.by.time), 1] == sim2$pops.by.time[nrow(sim2$pops.by.time), 1]
    )
  
  expect_false(
    sim1$pops.by.time[nrow(sim1$pops.by.time), 2] == sim2$pops.by.time[nrow(sim2$pops.by.time), 2]
  )
  
  expect_false(
    sim1$pops.by.time[nrow(sim1$pops.by.time), 3] == sim2$pops.by.time[nrow(sim2$pops.by.time), 3]
  )
  
  expect_false(
    sim1$pops.by.time[nrow(sim1$pops.by.time), 4] == sim2$pops.by.time[nrow(sim2$pops.by.time), 4]
  )
  
  expect_false(
    sim1$pops.by.time[nrow(sim1$pops.by.time), 5] == sim2$pops.by.time[nrow(sim2$pops.by.time), 5]
  )
  
})

test_that("eval mutator rate functionality with FDF situation and mu size > 1 at T>X", {
  
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
  
  afe_3_a <- allFitnessEffects(
    genotFitness =
      create_fe(2, 5, 1, 0.8, 1, 1, 9, 0),
    frequencyDependentFitness = TRUE,
    frequencyType = "rel")
  
  muvar2 <- c("Mph" = 1e-4, "BTC" = 1e-5, "MTC"=1e-5, "R" = 1e-6)
  
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
  
  muvar3 <- c("Mph" = 1e-5, "BTC" = 1e-6, "MTC"=1e-6, "R" = 1e-7)
  muexpression <- "if(T>=5) 10; else 1;"
  
  set.seed(2)
  s_3_b <- oncoSimulIndiv(afe_3_a,
                          model = "McFL", 
                          onlyCancer = FALSE, 
                          finalTime = 20,
                          mu = muvar3,
                          muFactor = muexpression,
                          initSize = 10000, 
                          keepPhylog = TRUE,
                          seed = NULL, 
                          errorHitMaxTries = FALSE, 
                          errorHitWallTime = FALSE)
  
  expect_false(s_3_a$TotalPopSize == s_3_b$TotalPopSize)
  
})