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

