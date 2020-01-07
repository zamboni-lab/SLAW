#include <Rcpp.h>
#include <map>
using namespace Rcpp;

// A script to merge clique if a bigger clique containing the package is detected.
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
List resize( const List& x, int n ){
  int oldsize = x.size() ;
  List y(n) ;
  for( int i=0; i<oldsize; i++) y[i] = x[i] ;
  return y ;
}

//Test of the C++ code to check overlap
// [[Rcpp::export]]
List mergeCliques(List cliques,List new_cliques,IntegerVector assignment,IntegerVector size, IntegerVector maxCliques){
  int max_clique = maxCliques[0];
//We convert all the lcique to set because it is easier.

  int percent = 0;
  //We build a percentage of progress
  for(int i=0;i<new_cliques.size();i++){

    NumericVector cc = new_cliques[i];
    if(cc.size()<=1){
      continue;
    }
    int cval = (i*10)/new_cliques.size();
    
    if(cval!=percent){
      percent = cval;
      Rcout << percent*10 << " ";
    }
    bool to_change = true;
    for(auto it = cc.begin();it != cc.end(); it++){
      if(size[*it-1]>cc.size()){
        to_change = false;
        break;
      }
    }

    //If any of the feature is already in a clique we d ont change it.
    if(!to_change){
      continue;
    }
    //We update the list if needed.
    max_clique++;
    if(max_clique>=cliques.size()){
      cliques = resize(cliques,cliques.size() +
        std::min(1000,(int)cliques.size()+1));
    }
   // tWe change the annotation of the clique components
   for(auto it = cc.begin();it != cc.end(); it++){
     assignment[*it-1] = max_clique;
     size[*it-1] = cc.size();
   }
   //We update the size with the new clique size
  cliques[max_clique] = cc;
  //We add the lement to the lcique list if necessary.
  }

  //We create a new index vector.

  std::vector<std::vector<int> > res_cpp;
  for(int i = 0;i < cliques.size(); i++){
    std::vector<int> temp;
    res_cpp.push_back(temp);
  }
  int cliques_counter = 0;
  //We create a new_index list to be sure that it works.
  std::map<int,int> idx_list;

  for(int i = 0; i <assignment.size();i++){
    //In this case we push a new cluster directly
    if(IntegerVector::is_na(assignment[i])){
      //We skip this assignment.
      continue;

    }else{
      int a = assignment[i];
      //We chekc if the cluster already exist or not
      auto it = idx_list.find(a);
      if(it == idx_list.end()){
        idx_list.emplace(a,cliques_counter);
        cliques_counter++;
      }
      res_cpp[idx_list[a]].push_back(i+1);
    }
    //We add the element to the clique
  }

  res_cpp.resize(cliques_counter);
  //We add the new values to the data.
  List res(res_cpp.size());
  for(int i = 0;i < res.size(); i++){
    NumericVector temp(res_cpp[i].begin(),res_cpp[i].end());
    res[i] = temp;
  }

  for(auto it=assignment.begin();it!=assignment.end();it++){
    if(IntegerVector::is_na(*it)) continue;
    //+1 for R compatibilty
    *it = idx_list[*it]+1;
  }
  maxCliques[0]=max_clique+1;
  return(res);
}
