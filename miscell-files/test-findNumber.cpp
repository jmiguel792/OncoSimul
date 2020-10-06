#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double findNumber(std::string& muExpression){
  
  std::vector<double> vtime;
  
  if(muExpression.find("T") != std::string::npos){
    //std::cout << "T found" << std::endl;
    std::string s = ">";
    //std::cout << "s found: " << s << std::endl;
    unsigned first = muExpression.find_first_of(s);
    unsigned end_pos = first + s.length();
    unsigned last = muExpression.find_first_of(")");
    double time = std::stod(muExpression.substr(end_pos, last-end_pos)); //string to double
    //std::cout << "time: " << time << std::endl;
    vtime.push_back(time);
  }
  
  //std::cout << "time: " << vtime[0] << std::endl;
  return vtime[0];
}
