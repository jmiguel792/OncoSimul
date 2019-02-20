//     Copyright 2013, 2014, 2015, 2016 Ramon Diaz-Uriarte

//     This program is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.

//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.

// This file is poorly named: these are functions in new_restric.cpp that
// are declared here so that new_restrict_1.cpp can compile.


#ifndef _NEW_RESTRICT_0_H__
#define _NEW_RESTRICT_0_H__

#include "new_restrict.h"
#include "debug_common.h"
#include "common_classes.h"
// #include "randutils.h" //Nope, until we have gcc-4.8 in Win; full C++11
#include <Rcpp.h>
#include <limits>
#include <random>

std::string concatIntsString(const std::vector<int>& ints,
			     const std::string sep = ", ");

void checkLegitGenotype(const Genotype& ge,
			const fitnessEffectsAll& F);


evalFVars_struct evalFVars(const fitnessEffectsAll& F,
			   const std::vector<Genotype>& Genotypes,
			   const std::vector<spParamsP>& popParams);

double totalPop(const std::vector<spParamsP>& popParams);

std::vector<int> GeneToModule(const std::vector<int>& Drv,
			     const
			      std::vector<Gene_Module_struct>& Gene_Module_tabl,
			      const bool sortout, const bool uniqueout);


std::vector<double> evalPosetConstraints(const std::vector<int>& mutatedModules,
					 const std::vector<Poset_struct>& Poset,
					 const std::vector<int>& allPosetG);


std::vector<double> evalEpistasis(const std::vector<int>& mutatedModules,
				  const std::vector<epistasis>& Epistasis);

std::vector<double> evalOrderEffects(const std::vector<int>& mutatedM,
				     const std::vector<epistasis>& OE);

// double evalGenotypeFDFitnessEcuation(const Genotype& ge,
// 	const fitnessEffectsAll& F,
// 	const std::vector<Genotype>& Genotypes,
// 				     const std::vector<spParamsP>& popParams);
#endif
