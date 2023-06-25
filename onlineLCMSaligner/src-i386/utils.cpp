#include <Rcpp.h>
#include <map>
using namespace Rcpp;


// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector findDuplicatePosition(NumericVector x) {
  NumericVector pos(x.size(),-1.0);
  std::set<double> used_elems;
  for(int i = 0; i < x.size(); i++){
    if(NumericVector::is_na(x[i])) continue;
    double tt= x[i];
    const bool is_in = used_elems.find(tt) != used_elems.end();
    if(is_in){
      pos[i]=x[i];
    }else{

      used_elems.insert(tt);
    }
  }
  return pos;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
set.seed(512)

svec <- 100000

tvec <- sample(1:10000,size=svec,replace=TRUE)

tvec[sample(seq_along(tvec),size = 500)] <- NA_real_

findDuplicatePosition(tvec)



*/
