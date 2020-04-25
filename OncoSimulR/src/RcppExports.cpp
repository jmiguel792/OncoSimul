// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// nr_BNB_Algo5
Rcpp::List nr_BNB_Algo5(Rcpp::List rFE, Rcpp::NumericVector mu_, double death, double initSize, double sampleEvery, double detectionSize, double finalTime, int initSp, int initIt, double seed, int verbosity, int speciesFS, double ratioForce, Rcpp::CharacterVector typeFitness_, int maxram, int mutationPropGrowth, Rcpp::IntegerVector initMutant_, double maxWallTime, double keepEvery, double K, int detectionDrivers, bool onlyCancer, bool errorHitWallTime, int maxNumTries, bool errorHitMaxTries, double minDetectDrvCloneSz, double extraTime, bool keepPhylog, Rcpp::List MMUEF, Rcpp::IntegerVector full2mutator_, double n2, double p2, double PDBaseline, double cPDetect_i, double checkSizePEvery, bool AND_DrvProbExit, Rcpp::List fixation_list);
RcppExport SEXP OncoSimulR_nr_BNB_Algo5(SEXP rFESEXP, SEXP mu_SEXP, SEXP deathSEXP, SEXP initSizeSEXP, SEXP sampleEverySEXP, SEXP detectionSizeSEXP, SEXP finalTimeSEXP, SEXP initSpSEXP, SEXP initItSEXP, SEXP seedSEXP, SEXP verbositySEXP, SEXP speciesFSSEXP, SEXP ratioForceSEXP, SEXP typeFitness_SEXP, SEXP maxramSEXP, SEXP mutationPropGrowthSEXP, SEXP initMutant_SEXP, SEXP maxWallTimeSEXP, SEXP keepEverySEXP, SEXP KSEXP, SEXP detectionDriversSEXP, SEXP onlyCancerSEXP, SEXP errorHitWallTimeSEXP, SEXP maxNumTriesSEXP, SEXP errorHitMaxTriesSEXP, SEXP minDetectDrvCloneSzSEXP, SEXP extraTimeSEXP, SEXP keepPhylogSEXP, SEXP MMUEFSEXP, SEXP full2mutator_SEXP, SEXP n2SEXP, SEXP p2SEXP, SEXP PDBaselineSEXP, SEXP cPDetect_iSEXP, SEXP checkSizePEverySEXP, SEXP AND_DrvProbExitSEXP, SEXP fixation_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::List >::type rFE(rFESEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< double >::type death(deathSEXP);
    Rcpp::traits::input_parameter< double >::type initSize(initSizeSEXP);
    Rcpp::traits::input_parameter< double >::type sampleEvery(sampleEverySEXP);
    Rcpp::traits::input_parameter< double >::type detectionSize(detectionSizeSEXP);
    Rcpp::traits::input_parameter< double >::type finalTime(finalTimeSEXP);
    Rcpp::traits::input_parameter< int >::type initSp(initSpSEXP);
    Rcpp::traits::input_parameter< int >::type initIt(initItSEXP);
    Rcpp::traits::input_parameter< double >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    Rcpp::traits::input_parameter< int >::type speciesFS(speciesFSSEXP);
    Rcpp::traits::input_parameter< double >::type ratioForce(ratioForceSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type typeFitness_(typeFitness_SEXP);
    Rcpp::traits::input_parameter< int >::type maxram(maxramSEXP);
    Rcpp::traits::input_parameter< int >::type mutationPropGrowth(mutationPropGrowthSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type initMutant_(initMutant_SEXP);
    Rcpp::traits::input_parameter< double >::type maxWallTime(maxWallTimeSEXP);
    Rcpp::traits::input_parameter< double >::type keepEvery(keepEverySEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type detectionDrivers(detectionDriversSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyCancer(onlyCancerSEXP);
    Rcpp::traits::input_parameter< bool >::type errorHitWallTime(errorHitWallTimeSEXP);
    Rcpp::traits::input_parameter< int >::type maxNumTries(maxNumTriesSEXP);
    Rcpp::traits::input_parameter< bool >::type errorHitMaxTries(errorHitMaxTriesSEXP);
    Rcpp::traits::input_parameter< double >::type minDetectDrvCloneSz(minDetectDrvCloneSzSEXP);
    Rcpp::traits::input_parameter< double >::type extraTime(extraTimeSEXP);
    Rcpp::traits::input_parameter< bool >::type keepPhylog(keepPhylogSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MMUEF(MMUEFSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type full2mutator_(full2mutator_SEXP);
    Rcpp::traits::input_parameter< double >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< double >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< double >::type PDBaseline(PDBaselineSEXP);
    Rcpp::traits::input_parameter< double >::type cPDetect_i(cPDetect_iSEXP);
    Rcpp::traits::input_parameter< double >::type checkSizePEvery(checkSizePEverySEXP);
    Rcpp::traits::input_parameter< bool >::type AND_DrvProbExit(AND_DrvProbExitSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type fixation_list(fixation_listSEXP);
    __result = Rcpp::wrap(nr_BNB_Algo5(rFE, mu_, death, initSize, sampleEvery, detectionSize, finalTime, initSp, initIt, seed, verbosity, speciesFS, ratioForce, typeFitness_, maxram, mutationPropGrowth, initMutant_, maxWallTime, keepEvery, K, detectionDrivers, onlyCancer, errorHitWallTime, maxNumTries, errorHitMaxTries, minDetectDrvCloneSz, extraTime, keepPhylog, MMUEF, full2mutator_, n2, p2, PDBaseline, cPDetect_i, checkSizePEvery, AND_DrvProbExit, fixation_list));
    return __result;
END_RCPP
}
// BNB_Algo5
Rcpp::List BNB_Algo5(Rcpp::IntegerMatrix restrictTable, int numDrivers, int numGenes, Rcpp::CharacterVector typeCBN_, double s, double death, double mu, double initSize, double sampleEvery, double detectionSize, double finalTime, int initSp, int initIt, int seed, int verbosity, int speciesFS, double ratioForce, Rcpp::CharacterVector typeFitness_, int maxram, int mutationPropGrowth, int initMutant, double maxWallTime, double keepEvery, double sh, double K, int detectionDrivers, bool onlyCancer, bool errorHitWallTime, int maxNumTries, bool errorHitMaxTries, double minDetectDrvCloneSz, double extraTime);
RcppExport SEXP OncoSimulR_BNB_Algo5(SEXP restrictTableSEXP, SEXP numDriversSEXP, SEXP numGenesSEXP, SEXP typeCBN_SEXP, SEXP sSEXP, SEXP deathSEXP, SEXP muSEXP, SEXP initSizeSEXP, SEXP sampleEverySEXP, SEXP detectionSizeSEXP, SEXP finalTimeSEXP, SEXP initSpSEXP, SEXP initItSEXP, SEXP seedSEXP, SEXP verbositySEXP, SEXP speciesFSSEXP, SEXP ratioForceSEXP, SEXP typeFitness_SEXP, SEXP maxramSEXP, SEXP mutationPropGrowthSEXP, SEXP initMutantSEXP, SEXP maxWallTimeSEXP, SEXP keepEverySEXP, SEXP shSEXP, SEXP KSEXP, SEXP detectionDriversSEXP, SEXP onlyCancerSEXP, SEXP errorHitWallTimeSEXP, SEXP maxNumTriesSEXP, SEXP errorHitMaxTriesSEXP, SEXP minDetectDrvCloneSzSEXP, SEXP extraTimeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type restrictTable(restrictTableSEXP);
    Rcpp::traits::input_parameter< int >::type numDrivers(numDriversSEXP);
    Rcpp::traits::input_parameter< int >::type numGenes(numGenesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type typeCBN_(typeCBN_SEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    Rcpp::traits::input_parameter< double >::type death(deathSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type initSize(initSizeSEXP);
    Rcpp::traits::input_parameter< double >::type sampleEvery(sampleEverySEXP);
    Rcpp::traits::input_parameter< double >::type detectionSize(detectionSizeSEXP);
    Rcpp::traits::input_parameter< double >::type finalTime(finalTimeSEXP);
    Rcpp::traits::input_parameter< int >::type initSp(initSpSEXP);
    Rcpp::traits::input_parameter< int >::type initIt(initItSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    Rcpp::traits::input_parameter< int >::type speciesFS(speciesFSSEXP);
    Rcpp::traits::input_parameter< double >::type ratioForce(ratioForceSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type typeFitness_(typeFitness_SEXP);
    Rcpp::traits::input_parameter< int >::type maxram(maxramSEXP);
    Rcpp::traits::input_parameter< int >::type mutationPropGrowth(mutationPropGrowthSEXP);
    Rcpp::traits::input_parameter< int >::type initMutant(initMutantSEXP);
    Rcpp::traits::input_parameter< double >::type maxWallTime(maxWallTimeSEXP);
    Rcpp::traits::input_parameter< double >::type keepEvery(keepEverySEXP);
    Rcpp::traits::input_parameter< double >::type sh(shSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type detectionDrivers(detectionDriversSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyCancer(onlyCancerSEXP);
    Rcpp::traits::input_parameter< bool >::type errorHitWallTime(errorHitWallTimeSEXP);
    Rcpp::traits::input_parameter< int >::type maxNumTries(maxNumTriesSEXP);
    Rcpp::traits::input_parameter< bool >::type errorHitMaxTries(errorHitMaxTriesSEXP);
    Rcpp::traits::input_parameter< double >::type minDetectDrvCloneSz(minDetectDrvCloneSzSEXP);
    Rcpp::traits::input_parameter< double >::type extraTime(extraTimeSEXP);
    __result = Rcpp::wrap(BNB_Algo5(restrictTable, numDrivers, numGenes, typeCBN_, s, death, mu, initSize, sampleEvery, detectionSize, finalTime, initSp, initIt, seed, verbosity, speciesFS, ratioForce, typeFitness_, maxram, mutationPropGrowth, initMutant, maxWallTime, keepEvery, sh, K, detectionDrivers, onlyCancer, errorHitWallTime, maxNumTries, errorHitMaxTries, minDetectDrvCloneSz, extraTime));
    return __result;
END_RCPP
}
// evalRGenotype
double evalRGenotype(Rcpp::IntegerVector rG, Rcpp::List rFE, bool verbose, bool prodNeg, Rcpp::CharacterVector calledBy_, double currentTime);
RcppExport SEXP OncoSimulR_evalRGenotype(SEXP rGSEXP, SEXP rFESEXP, SEXP verboseSEXP, SEXP prodNegSEXP, SEXP calledBy_SEXP, SEXP currentTimeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type rG(rGSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type rFE(rFESEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type prodNeg(prodNegSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type calledBy_(calledBy_SEXP);
    Rcpp::traits::input_parameter< double >::type currentTime(currentTimeSEXP);
    __result = Rcpp::wrap(evalRGenotype(rG, rFE, verbose, prodNeg, calledBy_, currentTime));
    return __result;
END_RCPP
}
// evalRGenotypeAndMut
Rcpp::NumericVector evalRGenotypeAndMut(Rcpp::IntegerVector rG, Rcpp::List rFE, Rcpp::List muEF, Rcpp::IntegerVector full2mutator_, bool verbose, bool prodNeg, double currentTime);
RcppExport SEXP OncoSimulR_evalRGenotypeAndMut(SEXP rGSEXP, SEXP rFESEXP, SEXP muEFSEXP, SEXP full2mutator_SEXP, SEXP verboseSEXP, SEXP prodNegSEXP, SEXP currentTimeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type rG(rGSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type rFE(rFESEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type muEF(muEFSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type full2mutator_(full2mutator_SEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type prodNeg(prodNegSEXP);
    Rcpp::traits::input_parameter< double >::type currentTime(currentTimeSEXP);
    __result = Rcpp::wrap(evalRGenotypeAndMut(rG, rFE, muEF, full2mutator_, verbose, prodNeg, currentTime));
    return __result;
END_RCPP
}

// evalRGenotypeAndMut
Rcpp::IntegerVector accessibleGenotypes(Rcpp::IntegerMatrix y, Rcpp::NumericVector f, Rcpp::IntegerVector numMut, double th);
RcppExport SEXP OncoSimulR_accessibleGenotypes(SEXP ySEXP, SEXP fSEXP, SEXP numMutSEXP, SEXP thSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
// Rcpp::RNGScope __rngScope;
 Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type y(ySEXP);
 Rcpp::traits::input_parameter< Rcpp::NumericVector >::type f(fSEXP);
 Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numMut(numMutSEXP);
 Rcpp::traits::input_parameter< double >::type th(thSEXP);
 __result = Rcpp::wrap(accessibleGenotypes(y, f, numMut, th));
 return __result;
 END_RCPP
}

// genotype fitness matrix to adjacency matrix of genotypes
Rcpp::NumericMatrix genot2AdjMat(Rcpp::IntegerMatrix y, Rcpp::NumericVector f, Rcpp::IntegerVector numMut);
RcppExport SEXP OncoSimulR_genot2AdjMat(SEXP ySEXP, SEXP fSEXP, SEXP numMutSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
// Rcpp::RNGScope __rngScope;
 Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type y(ySEXP);
 Rcpp::traits::input_parameter< Rcpp::NumericVector >::type f(fSEXP);
 Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numMut(numMutSEXP);
 __result = Rcpp::wrap(genot2AdjMat(y, f, numMut));
 return __result;
 END_RCPP
}


// evalRGenotypeAndMut
Rcpp::IntegerVector peaksLandscape(Rcpp::IntegerMatrix y, Rcpp::NumericVector f, Rcpp::IntegerVector numMut, double th);
RcppExport SEXP OncoSimulR_peaksLandscape(SEXP ySEXP, SEXP fSEXP, SEXP numMutSEXP, SEXP thSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
// Rcpp::RNGScope __rngScope;
 Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type y(ySEXP);
 Rcpp::traits::input_parameter< Rcpp::NumericVector >::type f(fSEXP);
 Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numMut(numMutSEXP);
 Rcpp::traits::input_parameter< double >::type th(thSEXP);
 __result = Rcpp::wrap(peaksLandscape(y, f, numMut, th));
 return __result;
 END_RCPP
}

// just for testing. Eventually remove
Rcpp::IntegerVector accessibleGenotypes_former(Rcpp::IntegerMatrix y, Rcpp::NumericVector f, Rcpp::IntegerVector numMut, double th);
RcppExport SEXP OncoSimulR_accessibleGenotypes_former(SEXP ySEXP, SEXP fSEXP, SEXP numMutSEXP, SEXP thSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
// Rcpp::RNGScope __rngScope;
 Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type y(ySEXP);
 Rcpp::traits::input_parameter< Rcpp::NumericVector >::type f(fSEXP);
 Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type numMut(numMutSEXP);
 Rcpp::traits::input_parameter< double >::type th(thSEXP);
 __result = Rcpp::wrap(accessibleGenotypes_former(y, f, numMut, th));
 return __result;
 END_RCPP
}


// // readFitnessEffects
// void readFitnessEffects(Rcpp::List rFE, bool echo);
// RcppExport SEXP OncoSimulR_readFitnessEffects(SEXP rFESEXP, SEXP echoSEXP) {
// BEGIN_RCPP
//     Rcpp::RNGScope __rngScope;
//     Rcpp::traits::input_parameter< Rcpp::List >::type rFE(rFESEXP);
//     Rcpp::traits::input_parameter< bool >::type echo(echoSEXP);
//     readFitnessEffects(rFE, echo);
//     return R_NilValue;
// END_RCPP
// }
