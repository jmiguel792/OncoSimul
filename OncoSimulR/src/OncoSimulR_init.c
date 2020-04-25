// This file was automatically generated by 'Kmisc::registerFunctions()'

#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP OncoSimulR_nr_BNB_Algo5(SEXP rFESEXP, SEXP mu_SEXP, SEXP deathSEXP, SEXP initSizeSEXP, SEXP sampleEverySEXP, SEXP detectionSizeSEXP, SEXP finalTimeSEXP, SEXP initSpSEXP, SEXP initItSEXP, SEXP seedSEXP, SEXP verbositySEXP, SEXP speciesFSSEXP, SEXP ratioForceSEXP, SEXP typeFitness_SEXP, SEXP maxramSEXP, SEXP mutationPropGrowthSEXP, SEXP initMutant_SEXP, SEXP maxWallTimeSEXP, SEXP keepEverySEXP,  SEXP KSEXP, SEXP detectionDriversSEXP, SEXP onlyCancerSEXP, SEXP errorHitWallTimeSEXP, SEXP maxNumTriesSEXP, SEXP errorHitMaxTriesSEXP, SEXP minDetectDrvCloneSzSEXP, SEXP extraTimeSEXP, SEXP keepPhylogSEXP, SEXP MMUEFSEXP, SEXP full2mutator_SEXP, SEXP n2SEXP, SEXP p2SEXP, SEXP PDBaselineSEXP, SEXP cPDetect_iSEXP, SEXP checkSizePEverySEXP, SEXP AND_DrvProbExitSEXP, SEXP fixation_listSEXP);

SEXP OncoSimulR_BNB_Algo5(SEXP restrictTableSEXP, SEXP numDriversSEXP, SEXP numGenesSEXP, SEXP typeCBN_SEXP,  SEXP sSEXP, SEXP deathSEXP, SEXP muSEXP, SEXP initSizeSEXP, SEXP sampleEverySEXP, SEXP detectionSizeSEXP, SEXP finalTimeSEXP, SEXP initSpSEXP, SEXP initItSEXP, SEXP seedSEXP, SEXP verbositySEXP, SEXP speciesFSSEXP, SEXP ratioForceSEXP, SEXP typeFitness_SEXP, SEXP maxramSEXP, SEXP mutationPropGrowthSEXP, SEXP initMutantSEXP, SEXP maxWallTimeSEXP, SEXP keepEverySEXP,  SEXP shSEXP, SEXP KSEXP, SEXP detectionDriversSEXP, SEXP onlyCancerSEXP, SEXP errorHitWallTimeSEXP, SEXP maxNumTriesSEXP, SEXP errorHitMaxTriesSEXP, SEXP minDetectDrvCloneSzSEXP, SEXP extraTimeSEXP);

SEXP OncoSimulR_evalRGenotype(SEXP rGSEXP, SEXP rFESEXP, SEXP verboseSEXP, SEXP prodNegSEXP, SEXP calledBy_SEXP, SEXP currentTimeSEXP);

SEXP OncoSimulR_evalRGenotypeAndMut(SEXP rGSEXP, SEXP rFESEXP, SEXP muEFSEXP, SEXP full2mutator_SEXP, SEXP verboseSEXP, SEXP prodNegSEXP, SEXP currentTimeSEXP);

// SEXP OncoSimulR_readFitnessEffects(SEXP rFESEXP, SEXP echoSEXP);
SEXP OncoSimulR_accessibleGenotypes(SEXP ySEXP, SEXP xSEXP, SEXP numMutSEXP, SEXP thSEXP);

SEXP OncoSimulR_genot2AdjMat(SEXP ySEXP, SEXP xSEXP, SEXP numMutSEXP);

SEXP OncoSimulR_peaksLandscape(SEXP ySEXP, SEXP xSEXP, SEXP numMutSEXP, SEXP thSEXP);

SEXP OncoSimulR_accessibleGenotypes_former(SEXP ySEXP, SEXP xSEXP, SEXP numMutSEXP, SEXP thSEXP);

// The number is the number of arguments
R_CallMethodDef callMethods[]  = {
  {"OncoSimulR_nr_BNB_Algo5", (DL_FUNC) &OncoSimulR_nr_BNB_Algo5, 37},
  {"OncoSimulR_BNB_Algo5", (DL_FUNC) &OncoSimulR_BNB_Algo5, 32},
  {"OncoSimulR_evalRGenotype", (DL_FUNC) &OncoSimulR_evalRGenotype, 6},
  {"OncoSimulR_evalRGenotypeAndMut", (DL_FUNC) &OncoSimulR_evalRGenotypeAndMut, 7},
  {"OncoSimulR_accessibleGenotypes", (DL_FUNC) &OncoSimulR_accessibleGenotypes, 4},
  {"OncoSimulR_genot2AdjMat", (DL_FUNC) &OncoSimulR_genot2AdjMat, 3},
  {"OncoSimulR_peaksLandscape", (DL_FUNC) &OncoSimulR_peaksLandscape, 4},
  {"OncoSimulR_accessibleGenotypes_former", (DL_FUNC) &OncoSimulR_accessibleGenotypes_former, 4},
  //  {"OncoSimulR_readFitnessEffects", (DL_FUNC) &OncoSimulR_readFitnessEffects, 2},
  {NULL, NULL, 0}
};

void R_init_OncoSimulR(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
