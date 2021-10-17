#include <RcppArmadillo.h>



/*
 * Normalize mutual information.
 */

// [[Rcpp::export]]
arma::sp_mat normalize_muti(arma::sp_mat muti,  // lower triangular
                            arma::vec entropy){
  for (int j = 0; j < muti.n_cols; j++){
    for (int i = j; i < muti.n_cols; i++){
      muti(i, j) = muti(i, j) / (entropy(i) + entropy(j) - muti(i, j));
      muti(i, j) = 1 - muti(i, j);
    }
  }
  return muti;
}



/*
 * Compute average distances between groups.
 */

arma::mat avg_dist(arma::mat& D,
                   arma::uvec& group,
                   arma::uvec& ugroup,
                   int linkage = 1){
  
  arma::mat D2(ugroup.size(), ugroup.size());
  D2.fill(2);
  
  for (int i = 0; i < ugroup.size()-1L; i++){
    
    arma::uvec I = arma::find(group == ugroup(i));
    
    for (int j = i+1L; j < ugroup.size(); j++){
      
      arma::uvec J = arma::find(group == ugroup(j));
      
      if (linkage == 0){
        D2(i, j) = (arma::accu(D(I, J)) / (I.size() * J.size()));
      } else if (linkage == 1){
        D2(i, j) = D(I, J).min();
      }
    }
  }
  
  return D2;
}



/*
 * Assign small groups into large groups or combine small groups.
 */

// [[Rcpp::export]]
arma::uvec assign_small(arma::mat D,
                        arma::uvec group,
                        arma::uvec ugroup,
                        arma::uvec sgroup,
                        double fraction = 0.05,
                        int linkage = 1){  // 0 for average linkage, 1 for single
  
  // fill_dist
  // arma::uvec ugroup = arma::unique(group);
  arma::mat D2(ugroup.size(), ugroup.size());
  if (ugroup.size() < D.n_cols){
    D2 = avg_dist(D, 
                  group,
                  ugroup,
                  linkage);
  } else{
    D2 = D;
    for (int j = 0; j < D.n_cols; j++){  // set upper tri of D2 to 2
      for (int i = j; i < D.n_rows; i++){
        D2(i, j) = 2;
      }
    }
  }
  
  int large = std::max(std::ceil(fraction * group.size()), 3.0);
  int best;
  arma::uvec sub(2);
  
  int i = 0, j = 0, sij = 0;
  int row = 0, col = 0;
  while (true){
    
    best = D2.index_min();
    if (D2(best) > 1) break;
    sub = arma::sort(arma::ind2sub(arma::size(D2), best));
    i = sub(0);  // smaller index
    j = sub(1);  // larger index
    
    // skip if both groups are large or already merged
    if (//(sgroup(i) + sgroup(j) > sgroup.max()) ||
        (sgroup(i) >= large && sgroup(j) >= large) ||
          (ugroup(i) == ugroup(j))){
      
      D2(i, j) = 2;
      continue;
    }
    
    // merge j into i (larger into smaller)
    
    // update groups
    for (int l = 0; l < group.size(); l++){
      
      if (group(l) == ugroup(j)){
        group(l) = ugroup(i);
      } else if (group(l) > ugroup(j)){
        group(l)--;
      }
    }
    
    // update distances
    for (int l = 0; l < D2.n_cols; l++){
      
      if (l == i || l == j) continue;
      
      if (l < i){
        row = l;
        col = i;
      } else {
        row = i;
        col = l;
      }
      
      if (linkage == 0){
        D2(row, col) = (D2(row, col) * sgroup(i) + std::min(D2(j, l), D2(l, j)) * sgroup(j)) / (sgroup(i) + sgroup(j));
      } else if (linkage == 1){
        D2(row, col) = std::min(D2(row, col), std::min(D2(j, l), D2(l, j)));
      }
    }
    D2.shed_row(j);
    D2.shed_col(j);
    
    // update labels
    ugroup.shed_row(j);
    for (int l = j; l < ugroup.size(); l++){
      ugroup(l)--;
    }
    
    // update sizes
    sij = sgroup(i) + sgroup(j);
    sgroup.shed_row(j);
    sgroup(i) = sij;
  }
  
  return group;
}