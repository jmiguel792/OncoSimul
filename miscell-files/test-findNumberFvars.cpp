#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string findVars(std::string& muExpression) {
  
  std::string fRelOrAbs;
  
  if(muExpression.find("f_") != std::string::npos){
    //std::cout << "f_ found" << std::endl;
    std::string s = "if(";
    unsigned first = muExpression.find_first_of(s);
    unsigned end_pos = first + s.length();
    unsigned last = muExpression.find_first_of(">");
    fRelOrAbs = muExpression.substr(end_pos, last-end_pos);
    std::string::iterator spaces = std::remove(fRelOrAbs.begin(), fRelOrAbs.end(), ' '); //remove whitespaces
    fRelOrAbs.erase(spaces, fRelOrAbs.end());
    
  } else if(muExpression.find("n_") != std::string::npos){
    //std::cout << "n_ found" << std::endl;
    std::string s = "if(";
    unsigned first = muExpression.find_first_of(s);
    unsigned end_pos = first + s.length();
    unsigned last = muExpression.find_first_of(">");
    fRelOrAbs = muExpression.substr(end_pos, last-end_pos);
    std::string::iterator spaces = std::remove(fRelOrAbs.begin(), fRelOrAbs.end(), ' '); //remove whitespaces
    fRelOrAbs.erase(spaces, fRelOrAbs.end());
    
  } else { std::cout << "nothing found" << std::endl; }
  
  //std::cout << "fRelOrAbs: " << fRelOrAbs << std::endl;
  return fRelOrAbs;
}

