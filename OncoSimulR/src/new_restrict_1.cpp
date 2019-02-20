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


// #include "randutils.h" //Nope, until we have gcc-4.8 in Win; full C++11
#include "debug_common.h"
#include "common_classes.h"
#include "new_restrict.h"
#include "new_restrict_0.h"
#include "exprtk.hpp"
#include <Rcpp.h>
#include <iomanip>
#include <algorithm>
#include <random>
#include <string>
#include <sstream>


using namespace Rcpp;
using std::vector;
using std::back_inserter;




double evalGenotypeFDFitnessEcuation(const Genotype& ge,
	const fitnessEffectsAll& F,
	const std::vector<Genotype>& Genotypes,
	const std::vector<spParamsP>& popParams){

  double f;

  evalFVars_struct symbol_table_struct = evalFVars(F, Genotypes, popParams);

	std::map<std::string, double> EFVMap = symbol_table_struct.evalFVarsmap;

  std::string gs = concatIntsString(ge.flGenes);

  std::string expr_string = F.fitnessLandscape.flFDFmap.at(gs);

	double N = totalPop(popParams);

  typedef exprtk::symbol_table<double> symbol_table_t;
  typedef exprtk::expression<double> expression_t;
  typedef exprtk::parser<double> parser_t;

	symbol_table_t symbol_table;
  for(auto& iterator : EFVMap){
		symbol_table.add_variable(iterator.first, iterator.second);
  }
	symbol_table.add_constant("N", N);//We reserve N to total population size
  symbol_table.add_constants();

  expression_t expression;
  expression.register_symbol_table(symbol_table);

	parser_t parser;

	if (!parser.compile(expr_string, expression)){
				Rcpp::Rcout << "\nexprtk parser error: \n" << std::endl;
	      for (std::size_t i = 0; i < parser.error_count(); ++i){
	         typedef exprtk::parser_error::type error_t;
	         error_t error = parser.get_error(i);
		 // RDU: FIXME?
	         // Rcpp::Rcout <<
		 REprintf("Error[%02zu] Position: %02zu Type: [%14s] Msg: %s Expression: %s\n",
			 i,
			 error.token.position,
			 exprtk::parser_error::to_str(error.mode).c_str(),
			 error.diagnostic.c_str(),
			 expr_string.c_str());
		 // << std::endl;
	      }
				std::string errorMessage1 = "Wrong evalGenotypeFDFitnessEcuation evaluation, ";
				std::string errorMessage2 = "probably bad fitness columm especification.";
				std::string errorMessage = errorMessage1 + errorMessage2;
				throw std::invalid_argument(errorMessage);
	  }

  f = expression.value();

  return f;
}

std::vector<double> evalGenotypeFitness(const Genotype& ge,
	const fitnessEffectsAll& F,
	const std::vector<Genotype>& Genotypes,
	const std::vector<spParamsP>& popParams){

  // check_disable_later
  checkLegitGenotype(ge, F);

  std::vector<double> s;

	if( ((ge.orderEff.size() + ge.epistRtEff.size() +
       ge.rest.size() + ge.flGenes.size() ) == 0) && !F.frequencyDependentFitness ) {
    Rcpp::warning("WARNING: you have evaluated fitness of a genotype of length zero.");
    // s.push_back(1.0); //Eh??!! 1? or 0? FIXME It should be empty! and have prodFitness
    // deal with it.
    return s;
  }

  // If we are dealing with a fitness landscape, that is as far as we go here
  // at least for now. No other genes affect fitness.
  // But this can be easily fixed in the future; do not return
  // s below, but keep adding, maybe the noIntGenes.
  // Recall also  prodFitness uses, well, the prod of 1 + s
  // so we want an s s.t. 1 + s = birth rate passed,
  // which is the value in the fitness landscape as interpreted now.
  // i.e., s = birth rate - 1;

  if(F.fitnessLandscape.NumID.size()) {
		std::string gs = concatIntsString(ge.flGenes);
		if(F.frequencyDependentFitness){//possible also with Genotype.size()==0 and popParams.size==0 ?
			if(F.fitnessLandscape.flFDFmap.find(gs) == F.fitnessLandscape.flFDFmap.end()) {
	    	s.push_back(-1.0);
			} else {
	      s.push_back(evalGenotypeFDFitnessEcuation(ge, F, Genotypes, popParams) - 1);
	    	}
		}else{
			if(F.fitnessLandscape.flmap.find(gs) == F.fitnessLandscape.flmap.end()) {
	   		s.push_back(-1.0);
	    } else {
	      s.push_back(F.fitnessLandscape.flmap.at(gs) - 1);
	    }
		}
	}

  // Genes without any restriction or epistasis are just genes. No modules.
  // So simple we do it here.
  if(F.genesNoInt.shift > 0) {
    int shift = F.genesNoInt.shift;
    for(auto const & r : ge.rest ) {
      s.push_back(F.genesNoInt.s[r - shift]);
    }
  }

  // For the rest, there might be modules. Three different effects on
  // fitness possible: as encoded in Poset, general epistasis, order effects.

  // Epistatis and poset are checked against all mutations. Create single
  // sorted vector with all mutations and map to modules, if needed. Then
  // eval.

  // Why not use a modified genotypeSingleVector without the no ints? We
  // could, but not necessary. And you can place genes in any order you
  // want, since this is not for order restrictions. That goes below.
  // Why do I put the epist first? See previous answer.
  // Why do I sort if one to one? binary searches. Not done below for order.
  std::vector<int> mutG (ge.epistRtEff);
  // A gene can be involved in epistasis and order. This gene would only
  // be in the orderEff vector, as seen in "createNewGenotype" or
  // "convertGenotypeFromInts"
  mutG.insert( mutG.end(), ge.orderEff.begin(), ge.orderEff.end());
  std::vector<int> mutatedModules;
  if(F.gMOneToOne) {
    sort(mutG.begin(), mutG.end());
    mutatedModules = mutG;
  } else {
    mutatedModules = GeneToModule(mutG, F.Gene_Module_tabl, true, true);
  }
  std::vector<double> srt =
    evalPosetConstraints(mutatedModules, F.Poset, F.allPosetG);
  std::vector<double> se =
    evalEpistasis(mutatedModules, F.Epistasis);

  // For order effects we need a new vector of mutatedModules:
  if(F.gMOneToOne) {
    mutatedModules = ge.orderEff;
  } else {
    mutatedModules = GeneToModule(ge.orderEff, F.Gene_Module_tabl, false, true);
  }

  std::vector<double> so =
    evalOrderEffects(mutatedModules, F.orderE);

  // I keep s, srt, se, so separate for now for debugging.
  s.insert(s.end(), srt.begin(), srt.end());
  s.insert(s.end(), se.begin(), se.end());
  s.insert(s.end(), so.begin(), so.end());

  return s;
}
