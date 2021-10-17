#include <RcppArmadillo.h>

using namespace Rcpp;

extern "C" {
#include "bnlearn_4.7/rcore.h"
#include "bnlearn_4.7/matrix.h"
#include "bnlearn_4.7/sets.h"
}

extern "C" {
  
  // ALLOCATIONS.C
  void *Calloc1D(size_t R, size_t size) {
    
    void *p = NULL;
    
    if (R == 0)
      return NULL;
    
    p = calloc(R, size);
    
    if (!p)
      Rf_error("unable to allocate a %d array.", R);
    
    return(p);
    
  }///*CALLOC1D*/
  
  void BN_Free1D(void *p) {
    
    free(p);
    
  }///*FREE1D*/
  
  void **Calloc2D(size_t R, size_t C, size_t size) {
    
    void **p = NULL;
    
    /* no corner cases, both dimensions required to be positive. */
    if ((R == 0) || (C == 0))
      Rf_error("trying to allocate a %dx%d two-dimensional array.", R, C);
    
    p = (void**) Calloc1D(R, sizeof(void *));
    
    for (int i = 0; i < R; i++)
      p[i] = Calloc1D(C, size);
    
    return p;
    
  }///*CALLOC2D*/
  
  void BN_Free2D(void **p, size_t R) {
    
    int i = 0;
    
    for (i = 0; i < R; i++)
      free(p[i]);
    free(p);
    
  }///*FREE2D*/
  
  /* initialize a one-dimensional contingency table. */
  int fill_1d_table(int *xx, int **n, int llx, int num) {
    
    int i = 0, ncomplete = 0;
    
    *n = (int*) Calloc1D(llx, sizeof(int));
    
    /* first fill the counts into the table... */
    for (i = 0; i < num; i++)
      if (xx[i] != NA_INTEGER)
        (*n)[xx[i] - 1]++;
      
      /* ... then add them up to count the number of complete observations. */
      for (i = 0; i < llx; i++)
        ncomplete += (*n)[i];
      
      return ncomplete;
      
  }///*FILL_1D_TABLE*/
  
  /* variadic version of mkString(). */
  SEXP mkStringVec(int n, ...) {
    
    va_list strings;
    int i = 0;
    SEXP vec;
    
    PROTECT(vec = Rf_allocVector(STRSXP, n));
    va_start(strings, n);
    for (i = 0; i < n; i++)
      SET_STRING_ELT(vec, i, Rf_mkChar(va_arg(strings, char *)));
    va_end(strings);
    UNPROTECT(1);
    
    return vec;
    
  }///*MKSTRINGVEC*/
  
  int c_has_path(int start, int stop, int *amat, int n, SEXP nodes,
                 int ugraph, int notdirect, int *path, int *counter, int debuglevel) {
    
    int i = 0, a1 = 0, a2 = 0, path_pos = 0, cur = start;
    
    /* remove any arc between start and stop if asked to. */
    if (notdirect) {
      
      a1 = amat[CMC(start, stop, n)];
      a2 = amat[CMC(stop, start, n)];
      amat[CMC(start, stop, n)] =  amat[CMC(stop, start, n)] = 0;
      
    }///*THEN*/
    
    /* initialize the position counters for the rows of the adjacency matrix. */
    memset(counter, '\0', n * sizeof(int));
    /* initialize the path array. */
    memset(path, '\0', n * sizeof(int));
    
    /* iterate until the other node is found. */
    while (cur != stop) {
      
      if (debuglevel > 0) {
        
        Rprintf("* currently at '%s'.\n", NODE(cur));
        Rprintf("  > current path is:\n");
        for (int i = 0; i < path_pos ; i ++)
          Rprintf("'%s' ", NODE(path[i]));
        Rprintf("'%s' \n", NODE(cur));
        
      }///*THEN*/
      
      there:
        
        /* find the next child of the 'cur' node. */
        for (i = 0; (i < n) && (counter[cur] < n); i++) {
          
          if (!ugraph) {
            
            if (amat[CMC(cur, counter[cur], n)] != 0)
              break;
            
          }///*THEN*/
          else {
            
            /* if we are looking for a path in the underlying (undirected)
             * graph, check also the symmetric entry of the adjacency matrix. */
            if ((amat[CMC(cur, counter[cur], n)] != 0) ||
            (amat[CMC(counter[cur], cur, n)] != 0))
              break;
            
          }///*ELSE*/
          
          counter[cur]++;
          
        }///*FOR*/
        
        /* the column indexes range from 0 to n - 1;  counter value of n means
         * that all the children of that node has beeen visited. */
        if (counter[cur] == n) {
          
          /* if this node is the first one in the path, there search is finished
           * and the 'stop' node was not found; return FALSE. */
          if  (path_pos == 0) {
            
            /* remove any arc between start and stop if asked to. */
            if (notdirect) {
              
              amat[CMC(start, stop, n)] = a1;
              amat[CMC(stop, start, n)] = a2;
              
            }///*THEN*/
            
            return FALSE;
            
          }///*THEN*/
          
          if (debuglevel > 0)
            Rprintf("  > node '%s' has no more children, going back to '%s'.\n",
                    NODE(cur), NODE(path[path_pos - 1]));
          
          /* this node has no more children, skip back to the previous one. */
          cur = path[--path_pos];
          path[path_pos + 1] = 0;
          
          goto there;
          
        }///*THEN*/
        else {
          
          /* increment the counter to get to the next node to check, unless
           * that one was the last. */
          if (counter[cur] < n)
            counter[cur]++;
          
          /* do not visit an already visited node */
          for (i = path_pos - 1; i >= 0; i--) {
            
            if ((counter[cur] - 1) == path[i]) {
              
              if (debuglevel > 0)
                Rprintf("  @ node '%s' already visited, skipping.\n", NODE(path[i]));
              
              goto there;
              
            }///*THEN*/
            
          }///*FOR*/
          
          /* update the path. */
          path[path_pos++] = cur;
          /* the current node is now the children we have just found. */
          cur = counter[cur] - 1;
          
          if (debuglevel > 0)
            Rprintf("  > jumping to '%s'.\n", NODE(cur));
          
        }///*ELSE*/
        
    }///*WHILE*/
    
    
    /* remove any arc between start and stop if asked to. */
    if (notdirect) {
      
      amat[CMC(start, stop, n)] = a1;
      amat[CMC(stop, start, n)] = a2;
      
    }///*THEN*/
    
    /* node 'stop' has been found, return TRUE. */
    return TRUE;
    
  }///*C_HAS_PATH*/
  
  bool sexp_has_path(int i, int j, SEXP amat, SEXP nodes) {
    
    bool indicate = false;
    int nnodes = Rf_length(nodes);
    int *am = NULL;
    int *path = NULL, *scratch = NULL;
    
    /* allocate buffers for c_has_path(). */
    path = (int*)Calloc1D(nnodes, sizeof(int));
    scratch = (int*)Calloc1D(nnodes, sizeof(int));
    
    /* save pointers to the numeric/integer matrices. */
    am = INTEGER(amat);
    
    indicate = c_has_path(i, j, am, nnodes, nodes, false, false, path, scratch,
                          false);
    
    Free1D(path);
    Free1D(scratch);
    
    return indicate;
  }///*SEXP_HAS_PATH*/
  
  double dlik(SEXP x, double *nparams) {
    
    int i = 0;
    int *n = NULL, *xx = INTEGER(x), llx = NLEVELS(x), num = Rf_length(x);
    double res = 0;
    
    /* initialize the contingency table. */
    fill_1d_table(xx, &n, llx, num);
    
    /* compute the entropy from the marginal frequencies. */
    for (i = 0; i < llx; i++)
      if (n[i] != 0)
        res += (double)n[i] * log((double)n[i] / num);
      
      /* we may want to store the number of parameters. */
      if (nparams)
        *nparams = llx - 1;
      
      Free1D(n);
      
      return res;
      
  }///*DLIK*/
  
  double cdlik(SEXP x, SEXP y, double *nparams) {
    
    int i = 0, j = 0, k = 0;
    int **n = NULL, *nj = NULL;
    int llx = NLEVELS(x), lly = NLEVELS(y), num = Rf_length(x);
    int *xx = INTEGER(x), *yy = INTEGER(y);
    double res = 0;
    
    /* initialize the contingency table and the marginal frequencies. */
    n = (int **) Calloc2D(llx, lly, sizeof(int));
    nj = (int*) Calloc1D(lly, sizeof(int));
    
    /* compute the joint frequency of x and y. */
    for (k = 0; k < num; k++)
      n[xx[k] - 1][yy[k] - 1]++;
    
    /* compute the marginals. */
    for (i = 0; i < llx; i++)
      for (j = 0; j < lly; j++)
        nj[j] += n[i][j];
    
    /* compute the conditional entropy from the joint and marginal
     frequencies. */
    for (i = 0; i < llx; i++)
      for (j = 0; j < lly; j++)
        if (n[i][j] != 0)
          res += (double)n[i][j] * log((double)n[i][j] / (double)nj[j]);
        
        /* we may want to store the number of parameters. */
        if (nparams)
          *nparams = (llx - 1) * lly;
        
        Free1D(nj);
        Free2D(n, llx);
        
        return res;
        
  }///*CDLIK*/
  
  /* return the unique elements from an input vector.*/
  SEXP unique(SEXP array) {
    
    int *d = NULL, i = 0, k = 0, dup_counter = 0, n = Rf_length(array);
    int *res = NULL, *a = NULL;
    SEXP dup, result = R_NilValue;
    
    PROTECT(dup = Rf_duplicated(array, FALSE));
    d = LOGICAL(dup);
    
    switch(TYPEOF(array)) {
    
    case INTSXP:
      
      a = INTEGER(array);
      
      for (i = 0; i < n; i++)
        if ((d[i] == 0) && (a[i] != NA_INTEGER))
          dup_counter++;
        
        PROTECT(result = Rf_allocVector(INTSXP, dup_counter));
        res = INTEGER(result);
        
        for (i = 0; i < n; i++)
          if ((d[i] == 0) && (a[i] != NA_INTEGER))
            res[k++] = a[i];
          
          break;
          
    case STRSXP:
      
      for (i = 0; i < n; i++)
        if (d[i] == 0)
          dup_counter++;
        
        PROTECT(result = Rf_allocVector(STRSXP, dup_counter));
        
        for (i = 0; i < n; i++)
          if (d[i] == 0)
            SET_STRING_ELT(result, k++, STRING_ELT(array, i));
          
          break;
          
    default:
      
      Rf_error("this SEXP type is not handled in unique().");
    
    }///*SWITCH*/
    
    UNPROTECT(2);
    
    return result;
    
  }///*UNIQUE*/
  
  /* get the list element named str, or return NULL. */
  SEXP getListElement(SEXP list, char *str) {
    
    SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);
    int i = 0;
    
    for (i = 0; i < Rf_length(list); i++) {
      
      if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
        
        elmt = VECTOR_ELT(list, i);
        break;
        
      }///*THEN*/
      
    }///*FOR*/
    
    return elmt;
    
  }/*GETLISTELEMENT*/
  
  /* transform an integer vector into a factor. */
  SEXP int2fac(SEXP vector, int *nlevels) {
    
    int i = 0, *l = NULL, *r = NULL, *v = INTEGER(vector);
    SEXP result, levels, lvls;
    
    if (!nlevels) {
      
      PROTECT(levels = unique(vector));
      
    }///*THEN*/
    else {
      
      PROTECT(levels = Rf_allocVector(INTSXP, *nlevels));
      l = INTEGER(levels);
      
      for (i = 0; i < *nlevels; i++)
        l[i] = i;
      
    }///*ELSE*/
    
    /* match the elements of the vector against the levels, preserving NAs and
     * setting to NA those elements that cannot be matched to a level. */
    PROTECT(result = Rf_match(levels, vector, 0));
    r = INTEGER(result);
    
    for (i = 0; i < Rf_length(result); i++)
      if ((r[i] == 0) || (v[i] == NA_INTEGER))
        r[i] = NA_INTEGER;
      
      /* set the levels of the factor. */
      PROTECT(lvls = Rf_coerceVector(levels, STRSXP));
      Rf_setAttrib(result, R_LevelsSymbol, lvls);
      
      /* set the class of the return value. */
      Rf_setAttrib(result, R_ClassSymbol, Rf_mkString("factor"));
      
      UNPROTECT(3);
      
      return result;
      
  }///*INT2FAC*/
  
  /* identify different configurations of factors and assign them unique integer
   * codes, computed using the general formula for the column-mayor indexing. */
  void cfg(SEXP parents, int *configurations, int *nlevels) {
    
    int i = 0, **columns = NULL, *levels = NULL;
    int ncol = Rf_length(parents), nrow = Rf_length(VECTOR_ELT(parents, 0));
    SEXP temp;
    
    /* dereference the columns of the data frame. */
    columns = (int **) Calloc1D(ncol, sizeof(int *));
    levels = (int*)Calloc1D(ncol, sizeof(int));
    for (i = 0; i < ncol; i++) {
      
      temp = VECTOR_ELT(parents, i);
      columns[i] = INTEGER(temp);
      levels[i] = NLEVELS(temp);
      
    }///*FOR*/
    
    c_fast_config(columns, nrow, ncol, levels, configurations, nlevels, 0);
    
    Free1D(columns);
    Free1D(levels);
    
  }///*CFG*/
  
  void c_fast_config(int **columns, int nrow, int ncol, int *levels, int *configurations,
                     int *nlevels, int offset) {
    
    int i = 0, j = 0, cfgmap = 0;
    long long *cumlevels = NULL, nl = 0;
    
    /* create the cumulative products of the number of levels. */
    cumlevels = (long long*) Calloc1D(ncol, sizeof(long long));
    
    /* set the first one to 1 ... */
    cumlevels[0] = 1;
    
    /* ... then compute the following ones. */
    for (j = 1; j < ncol; j++)
      cumlevels[j] = cumlevels[j - 1] * levels[j - 1];
    
    /* compute the number of possible configurations. */
    nl = cumlevels[ncol - 1] * levels[ncol - 1];
    
    if (nl >= INT_MAX)
      Rf_error("attempting to create a factor with more than INT_MAX levels.");
    
    /* if nlevels is not a NULL pointer, save the number of possible
     * configurations. */
    if (nlevels)
      *nlevels = nl;
    
    for (i = 0; i < nrow; i++) {
      
      /* reset the configuration mapping of the new row. */
      cfgmap = 0;
      
      for (j = 0; j < ncol; j++) {
        
        if (columns[j][i] == NA_INTEGER) {
          
          cfgmap = NA_INTEGER;
          break;
          
        }///*THEN*/
        else {
          
          cfgmap += (columns[j][i] - 1) * cumlevels[j];
          
        }///*ELSE*/
        
      }///*FOR*/
      
      /* save the configuration, applying the offset to non-NA values. */
      configurations[i] = cfgmap + offset * (cfgmap != NA_INTEGER);
      
    }///*FOR*/
    
    Free1D(cumlevels);
    
  }///*C_FAST_CONFIG*/
  
  /* wrapper around the cfg() function for use in C code. */
  SEXP c_configurations(SEXP parents, int factor, int all_levels) {
    
    int i = 0, *res = NULL, nlevels = 0;
    SEXP temp, result;
    
    /* compute the configurations. */
    PROTECT(temp = Rf_allocVector(INTSXP, Rf_length(VECTOR_ELT(parents, 0))));
    cfg(parents, INTEGER(temp), &nlevels);
    
    if (factor) {
      
      /* convert the configurations from an integer array to to a factor. */
      if (all_levels)
        result = int2fac(temp, &nlevels);
      else
        result = int2fac(temp, NULL);
      
    }///*THEN*/
    else {
      
      /* keep the configurations as they are, but add 1 to get indexing right
       * in R. */
      result = temp;
      res = INTEGER(result);
      
      for (i = 0; i < Rf_length(result); i++)
        if (res[i] != NA_INTEGER)
          res[i]++;
        
    }///*ELSE*/
    
    UNPROTECT(1);
    
    return result;
    
  }///*C_CONFIGURATIONS*/
  
  SEXP c_dataframe_column(SEXP dataframe, SEXP name, int drop, int keep_names) {
    
    SEXP try2, result, colnames = Rf_getAttrib(dataframe, R_NamesSymbol);
    int *idx = NULL, nnames = Rf_length(name), name_type = TYPEOF(name);
    
    if (dataframe == R_NilValue)
      return R_NilValue;
    
    switch(name_type) {
    
    case STRSXP:
      
      /* column names passed as strings; match the corresponding indexes. */
      PROTECT(try2 = Rf_match(colnames, name, 0));
      idx = INTEGER(try2);
      break;
      
    case REALSXP:
      
      /* these are almost good enough, coerce them to integers. */
      PROTECT(try2 = Rf_coerceVector(name, INTSXP));
      idx = INTEGER(try2);
      break;
      
    case INTSXP:
      
      /* these are already indexes, nothing to do. */
      idx = INTEGER(name);
      break;
      
    default:
      
      Rf_error("this SEXP type is not handled in minimal.data.frame.column().");
    
    }///*SWITCH*/
    
    if ((nnames > 1) || (drop == 0)) {
      
      PROTECT(result = Rf_allocVector(VECSXP, nnames));
      
      for (int i = 0; i < nnames; i++)
        SET_VECTOR_ELT(result, i, VECTOR_ELT(dataframe, idx[i] - 1));
      
      if (keep_names)
        Rf_setAttrib(result, R_NamesSymbol, name);
      
      UNPROTECT(1);
      
    }///*THEN*/
    else {
      
      if (*idx != 0)
        result = VECTOR_ELT(dataframe, *idx - 1);
      else
        result = R_NilValue;
      
    }///*ELSE*/
    
    if (name_type != INTSXP)
      UNPROTECT(1);
    
    return result;
    
  }///*C_DATAFRAME_COLUMN*/
  
  double c_loglik_dnode(SEXP target, SEXP parents, SEXP data, double *nparams, int debuglevel) {
    
    // double *nparams = 0;
    double loglik = 0;
    char *t = (char *)CHAR(STRING_ELT(target, 0));
    SEXP nodes, node_t, data_t, parent_vars, config;
    
    /* get the node cached information. */
    // nodes = getListElement(x, "nodes");
    // node_t = getListElement(nodes, t);
    /* get the parents of the node. */
    // parents = getListElement(node_t, "parents");
    /* extract the node's column from the data frame. */
    PROTECT(data_t = c_dataframe_column(data, target, TRUE, FALSE));
    
    if (Rf_length(parents) == 0) {
      
      loglik = dlik(data_t, nparams);
      
    }///*THEN*/
    else {
      
      /* generate the configurations of the parents. */
      PROTECT(parent_vars = c_dataframe_column(data, parents, FALSE, FALSE));
      PROTECT(config = c_configurations(parent_vars, TRUE, TRUE));
      /* compute the log-likelihood. */
      loglik = cdlik(data_t, config, nparams);
      
      UNPROTECT(2);
      
    }///*ELSE*/
    
    if (debuglevel > 0)
      Rprintf("  > loglikelihood is %lf.\n", loglik);
    
    UNPROTECT(1);
    
    return loglik;
    
  }///*LOGLIK_DNODE*/
  
}



/*
 * Check whether there is a path from i to j in amat.
 */

// [[Rcpp::export]]
bool cpp_has_path(int i, int j, SEXP amat, SEXP nodes){
  return sexp_has_path(i, j, amat, nodes);
}



/*
 * Check if x is in X.
 */

bool is_in_ivec(int x, arma::ivec& X){
  return(std::find(X.begin(), X.end(), x) != X.end());
}




/*
 * Check if any element in row i of vs is in X.
 */

bool vs_is_in_ivec(int i, arma::mat& vs, arma::ivec& X){
  int x;
  for (int j = 1; j < vs.n_cols; j++){
    x = vs(i, j);
    if (std::find(X.begin(), X.end(), x) != X.end())
      return true;
  }
  return false;
}



/*
 * Reference replacement of original into copy, allowing for
 * deep copy by creating copy and replacing with original.
 */

// [[Rcpp::export]]
void deep_copy_NumericVector(Rcpp::NumericVector& original,
                             Rcpp::NumericVector& copy){
  std::copy(original.begin(), original.end(), copy.begin());
}



/*
 * Rcpp implementation of bnlearn:::vstruct.apply() for speed.
 */

// [[Rcpp::export]]
Rcpp::IntegerMatrix vstruct_apply_cpp(Rcpp::IntegerMatrix undirMat,
                                      arma::mat& vs,
                                      Rcpp::CharacterVector nodes,
                                      int maxp = 8, 
                                      bool debug = false){
  
  // arma::umat dirMat(undirMat.n_rows, undirMat.n_cols, arma::fill::zeros);
  Rcpp::IntegerMatrix dirMat(undirMat.nrow(), undirMat.ncol());
  
  double max_a;
  int y, x, z;
  
  if (debug)
    Rcpp::Rcout << "----------------------------------------------------------------\n";
  for (int i = 0; i < vs.n_rows; i++){
    
    max_a = vs(i, 0);
    y = vs(i, 1);
    x = vs(i, 2);
    z = vs(i, 3);
    
    if (!((undirMat(y, x) || dirMat(y, x)) && 
        (undirMat(z, x) || dirMat(z, x)))){
      
      if (debug)
        Rcpp::Rcout << "* not applying v-structure " << y << " -> " << x <<
          " <- " << z << " (" << max_a << ")\n";
      Rf_warning("vstructure is not applicable, because one or both arcs are oriented in the opposite direction.",
                 y, x, z);
      continue;
    }
    
    // Rcpp::IntegerMatrix dirIntMat = as<Rcpp::IntegerMatrix>(Rcpp::wrap(dirMat));
    if (cpp_has_path(x, y, dirMat, nodes) ||
        cpp_has_path(x, z, dirMat, nodes)){
      
      if (debug)
        Rcpp::Rcout << "* not applying v-structure " << y << " -> " << x <<
          " <- " << z << " (" << max_a << ")\n";
      Rf_warning("vstructure is not applicable, because one or both arcs introduce cycles in the graph.", 
                 y, x, z);
      continue;
    }
    
    dirMat(y, x) = 1L;
    dirMat(z, x) = 1L;
    
    if (Rcpp::sum(dirMat(_, x)) > maxp){
      
      if (undirMat(y, x) > 0)
        dirMat(y, x) = 0L;
      if (undirMat(z, x) > 0)
        dirMat(z, x) = 0L;
      
      if (debug)
        Rcpp::Rcout << "* not applying v-structure " << y << " -> " << x <<
          " <- " << z << " (" << max_a << ")\n";
      Rf_warning("vstructure is not applicable, because exceeds maximum number of parents.",
                 y, x, z);
      continue;
    }
    
    undirMat(x, y) = 0L;
    undirMat(y, x) = 0L;
    undirMat(x, z) = 0L;
    undirMat(z, x) = 0L;
    
    if (debug)
      Rcpp::Rcout << "* applying v-structure " << y << " -> " << x << " <- " << z <<
        " (" << max_a << ")\n";
    
    vs(i, 0) = -1;
  }
  
  dirMat += undirMat;
  
  return dirMat;
}



/*
 * Score functions
 */

// [[Rcpp::export]]
double R_loglik_dnode(SEXP target, SEXP parents, SEXP data, double k, int debuglevel = 0){
  
  double nparams = 0;
  double loglik = c_loglik_dnode(target, 
                                 parents,
                                 data,
                                 &nparams,
                                 debuglevel);
  
  return loglik - (k * nparams);
}

double cpp_loglik_dnode(int target, 
                        Rcpp::List parentList, 
                        Rcpp::CharacterVector nodes, 
                        SEXP data, 
                        double k, 
                        int debuglevel = 0){
  
  Rcpp::IntegerVector t = Rcpp::IntegerVector::create(target);
  Rcpp::IntegerVector parents = parentList[target];
  return R_loglik_dnode(nodes[t], 
                        nodes[parents],
                             data, 
                             k, 
                             debuglevel);
}

double local_bic_cpp(Rcpp::IntegerVector t,  // int target, 
                     Rcpp::IntegerVector parents,  // Rcpp::List parentList, 
                     Rcpp::CharacterVector nodes, 
                     SEXP data, 
                     double k,
                     int is_discrete = 1,
                     int debuglevel = 0){
  
  // originally intended to switch between
  // discrete and Gaussian
  if (is_discrete){
    return R_loglik_dnode(nodes[t], 
                          nodes[parents],
                               data, 
                               k, 
                               debuglevel);
  } else{
    return -1;
  }
}




/*
 * Indicate conflicts induced by applying v-structure i.
 */

void indicate_conflicts(int i,
                        arma::mat& vs,
                        arma::uvec& ignore,
                        arma::imat& conflicts){
  // Rcpp::Rcout << "* indicating conflict against " << vs(i, 1) << " -> " << 
  //   vs(i, 2) << " <- " << vs(i, 3) << "\n";
  for (int j = 0; j < vs.n_rows; j++){
    if (ignore(j) == 0 &&  // j is not oriented
        (vs(j, 2) == vs(i, 1) || vs(j, 2) == vs(i, 3)) &&  // j is centered at y or z
        (vs(j, 1) == vs(i, 2) || vs(j, 3) == vs(i, 2))){  // j is oriented from x
      conflicts(i, j) = 1;
      // Rcpp::Rcout << "* indicating conflict with " << vs(j, 1) << " -> " << 
      //   vs(j, 2) << " <- " << vs(j, 3) << "\n";
    }
  }
}




/*
 * Update conflicts, for use if reversals are allowed (not checked).
 */

void update_conflicts(Rcpp::IntegerMatrix& dirMat,
                      arma::mat& vs,
                      arma::uvec& conflicts_vec,
                      arma::uvec& ignore,
                      arma::ivec& updated,
                      bool remove = true,
                      bool update = false,
                      bool debug = false){
  
  // arma::uvec conflicts_vec = arma::find(conflicts.col(i) > 0);
  for (int k = 0; k < conflicts_vec.size(); k++){
    
    if (remove){
      dirMat(vs(conflicts_vec(k), 1), vs(conflicts_vec(k), 2)) = 0L;
      dirMat(vs(conflicts_vec(k), 3), vs(conflicts_vec(k), 2)) = 0L;
      ignore(conflicts_vec(k)) = 0;
      if (debug && update)
        Rcpp::Rcout << "* removing v-structure " << vs(conflicts_vec(k), 1) << " -> " << 
          vs(conflicts_vec(k), 2) << " <- " << vs(conflicts_vec(k), 3) << "\n";
    } else{
      dirMat(vs(conflicts_vec(k), 1), vs(conflicts_vec(k), 2)) = 1L;
      dirMat(vs(conflicts_vec(k), 3), vs(conflicts_vec(k), 2)) = 1L;
      ignore(conflicts_vec(k)) = 2;
    }
    
    if (update){
      updated(k+1) = vs(conflicts_vec(k), 2);
    }
  }
}




/*
 * Get score difference, for use if reversals are allowed (not checked). 
 */

double get_delta(arma::vec& delta,
                 int i, 
                 arma::uvec& conflicts_vec){
  double temp_delta = delta(i);
  for (int k = 0; k < conflicts_vec.size(); k++){
    temp_delta += delta(conflicts_vec(k));
  }
  return temp_delta;
}



/*
 * Cache v-structure score differences.
 */

void vs_cache_fill(Rcpp::IntegerMatrix& dirMat,
                   arma::vec& reference,
                   Rcpp::CharacterVector& nodes,
                   arma::mat& vs,
                   arma::vec& delta,
                   arma::uvec& ignore,
                   SEXP data,
                   int& nscores,
                   arma::ivec& updated,  // int updated = -1,
                   arma::imat& conflicts,
                   int reverse = 0,
                   int is_discrete = 1,
                   double k = 0,
                   bool debug = false){
  
  constexpr double DOUBLE_MIN = 0;  // std::numeric_limits<double>::lowest();
  // constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();
  
  double temp_score = 0, max_a = 0;
  int x, y, z;
  Rcpp::IntegerVector t(1);
  Rcpp::IntegerVector indices = Rcpp::seq(0, nodes.size()-1);
  Rcpp::LogicalVector bool_parents;
  arma::uvec conflicts_vec;
  
  if (updated(0) == -2) return;
  
  for (int i = 0; i < vs.n_rows; i++){
    
    max_a = vs(i, 0);
    x = vs(i, 2);
    y = vs(i, 1);
    z = vs(i, 3);
    
    if (ignore(i) || 
        (reverse == 0 && updated(0) >= 0 && !is_in_ivec(x, updated)) ||
        (reverse > 0 && updated(0) >= 0 && !vs_is_in_ivec(i, vs, updated))
    ) continue;
    
    if (debug)
      Rcpp::Rcout << "* scoring v-structure " << y << " -> " << x << " <- " << z <<
        " (" << delta(i) << ")\n";
    
    t(0) = x;
    bool_parents = (dirMat(_, t(0)) > 0);
    bool_parents(z) = true;  // direct z -> x
    bool_parents(y) = true;  // direct z -> y
    // temp_score = cpp_loglik_dnode2(t, indices[bool_parents],
    //                                nodes, data, k, 1 * debug);
    temp_score = local_bic_cpp(t, indices[bool_parents],
                               nodes, data, k, is_discrete, 1 * debug);
    nscores++;
    delta(i) = temp_score - reference(t(0));
    
    if (reverse > 0 && arma::any(conflicts.col(i) > 0)){  // are any conflicts
      conflicts_vec = arma::find(conflicts.col(i) > 0);
      
      // drop v-structures centered at y and then z
      int j = 1;
      while (j <= 3){
        
        if (dirMat(x, vs(i, j))){  // oriented x to y or x to z
          
          if (debug)
            Rcpp::Rcout << "* additionally scoring node " << vs(i, j) << "\n";
          
          t(0) = vs(i, j);  // y (j==1) and z (j==3)
          bool_parents = (dirMat(_, t(0)) > 0);  // parents of y or z
          for (int k = 0; k < conflicts_vec.size(); k++){
            
            // if centered at y or z, remove v-structure (should always be one of y or z)
            if (vs(conflicts_vec(k), 2) == t(0)){
              bool_parents(vs(conflicts_vec(k), 1)) = false;
              bool_parents(vs(conflicts_vec(k), 3)) = false;
            }
          }
          // temp_score = cpp_loglik_dnode2(t, indices[bool_parents], 
          //                                nodes, data, k, 1 * debug);
          temp_score = local_bic_cpp(t, indices[bool_parents],
                                     nodes, data, k, is_discrete, 1 * debug);
          nscores++;
          for (int k = 0; k < conflicts_vec.size(); k++){
            if (vs(conflicts_vec(k), 2) == t(0)){
              delta(conflicts_vec(k)) = temp_score - reference(t(0));  // change in score of y or z
              break;
            }
          }
        }
        j += 2;  // move from y to z
      }
    }
  }
}



/*
 * Operation step for greedily applying v-structures, somewhat
 * modeled after bnlearn:::call_hc_opt_step.
 */

void vs_opt_step(arma::ivec& updated,
                 arma::imat& conflicts,
                 Rcpp::IntegerMatrix& undirMat,
                 Rcpp::IntegerMatrix& dirMat,
                 arma::vec& reference,
                 Rcpp::CharacterVector& nodes,
                 arma::mat& vs,
                 arma::vec& delta,
                 arma::uvec& ignore,
                 int reverse = 0,
                 int maxp = 8,
                 bool debug = false){
  
  constexpr double DOUBLE_MIN = 0;  // std::numeric_limits<double>::lowest();
  // constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();
  const double tol = arma::datum::eps * 2.0;
  
  double temp_delta = 0, best_delta = DOUBLE_MIN, max_a = 0;
  int x, y, z, best_i = -1;
  
  arma::uvec conflicts_vec;
  for (int i = 0; i < vs.n_rows; i++){
    
    if (ignore(i)) continue;
    
    max_a = vs(i, 0);
    y = vs(i, 1);
    x = vs(i, 2);
    z = vs(i, 3);
    
    if (reverse == 0 && 
        (dirMat(x, y) ||
        dirMat(x, z))){
      
      if (debug)
        Rcpp::Rcout << "* not applying v-structure " << y << " -> " << x <<
          " <- " << z << " (" << delta(i) << ")\n";
      Rf_warning("vstructure is not applicable, because one or both arcs are oriented in the opposite direction.",
                 y, x, z);
      
      delta(i) = DOUBLE_MIN;
      ignore(i) = 1;
      continue;
    }
    
    if (reverse == 0 && 
        (delta(i) - best_delta) < tol) continue;
    
    if (maxp < reference.size() &&
        (Rcpp::sum(dirMat(_, x)) + (dirMat(y, x) == 0) + (dirMat(z, x) == 0)) > maxp){
      
      if (debug)
        Rcpp::Rcout << "* not applying v-structure " << y << " -> " << x <<
          " <- " << z << " (" << delta(i) << ")\n";
      Rf_warning("vstructure is not applicable, because exceeds maximum number of parents.",
                 y, x, z);
      
      delta(i) = DOUBLE_MIN;
      ignore(i) = 1;
      continue;
    }
    
    if (reverse > 0){
      
      // get conflicts
      conflicts_vec = arma::find(conflicts.col(i) > 0);
      temp_delta = get_delta(delta, i, conflicts_vec);
      
      if ((temp_delta - best_delta) < tol) continue;
      
      // if too many conflicts
      if (conflicts_vec.size() > reverse){
        
        if (debug)
          Rcpp::Rcout << "* not applying v-structure " << y << " -> " << x <<
            " <- " << z << " (" << delta(i) << ")\n";
        Rf_warning("vstructure is not applicable, because too many conflicting v-structures.",
                   y, x, z);
        
        // delta(i) = DOUBLE_MIN;
        // ignore(i) = 1;
        continue;
      }
    }
    
    // if induces a cycle
    if (reverse > 0)
      update_conflicts(dirMat, vs, conflicts_vec, ignore, updated, true, false);  // remove conflicts
    if (cpp_has_path(x, y, dirMat, nodes) ||
        cpp_has_path(x, z, dirMat, nodes)){
      
      if (debug)
        Rcpp::Rcout << "* not applying v-structure " << y << " -> " << x <<
          " <- " << z << " (" << delta(i) << ")\n";
      Rf_warning("vstructure is not applicable, because one or both arcs introduce cycles in the graph.",
                 y, x, z);
      
      if (reverse > 0){
        update_conflicts(dirMat, vs, conflicts_vec, ignore, updated, false, false);  // restore conflicts
      }
      else{
        delta(i) = DOUBLE_MIN;
        ignore(i) = 1;
      }
      continue;
    }
    if (reverse > 0)
      update_conflicts(dirMat, vs, conflicts_vec, ignore, updated, false, false);  // restore conflicts
    
    best_delta = delta(i);
    best_i = i;
  }
  
  updated.fill(-1);
  if (best_i < 0)
    return;
  
  int i = best_i;
  max_a = vs(i, 0);
  y = vs(i, 1);
  x = vs(i, 2);
  z = vs(i, 3);
  
  if (debug)
    Rcpp::Rcout << "* applying v-structure " << y << " -> " << x << " <- " << z <<
      " (" << delta(i) << ")\n";
  
  updated(0) = x;  // indicate updated
  reference(x) += delta(i);  // update node score for x
  
  if (reverse > 0){
    
    // remove conflicts to v-structure i
    conflicts_vec = arma::find(conflicts.col(i) > 0);
    update_conflicts(dirMat, vs, conflicts_vec, ignore, updated, true, true);  // remove in graph and update updated
    // conflicts.col(i).zeros();  // remove indication of conflicts on v-structure i
    for (int k = 0; k < conflicts_vec.size(); k++){
      conflicts.row(conflicts_vec(k)).zeros();  // remove conflicts by removed v-structures
    }
    
    // indicate conflicts induced by the applied v-structure i
    indicate_conflicts(i, vs, ignore, conflicts);
    
    // update node scores for y and z
    int j = 1;
    while (j <= 3){
      for (int k = 0; k < conflicts_vec.size(); k++){
        if (vs(conflicts_vec(k), 2) == vs(i, j)){  // first conflicting vs centered at y or z
          reference(vs(conflicts_vec(k), 2)) += delta(conflicts_vec(k));
          break;
        }
      }
      j += 2;
    }
    
    if (debug)
      Rcpp::Rcout << "reference = " << arma::sum(reference) << " (+" << 
        get_delta(delta, i, conflicts_vec) << ")" << std::endl;
  }
  
  // orient v-structure i
  dirMat(y, x) = dirMat(z, x) = 1L;
  // undirMat(x, y) = undirMat(y, x) = undirMat(x, z) = undirMat(z, x) = 0L;
  
  if (reverse == 0 && debug)
    Rcpp::Rcout << "reference = " << arma::sum(reference) << " (+" << 
      delta(i) << ")" << std::endl;
  
  delta(i) = DOUBLE_MIN;
  ignore(i) = 2;
}



/*
 * Greedily apply v-structures for the hybrid greedy
 * initialization (HGI) algorithm. Modeled somewhat after
 * bnlearn:::vstruct.apply(). 
 */

// [[Rcpp::export]]
List vstruct_apply_hgi(Rcpp::IntegerMatrix undirMat,
                       arma::vec& reference,
                       Rcpp::CharacterVector& nodes,
                       arma::mat& vs,
                       SEXP data,
                       arma::vec& delta,
                       int reverse = 0,
                       int is_discrete = 1,
                       double k = 0,
                       int maxp = 8, 
                       bool just_delta = false,
                       bool debug = false){
  
  Rcpp::IntegerMatrix dirMat(undirMat.nrow(), undirMat.ncol());
  arma::uvec ignore(vs.n_rows, arma::fill::zeros);
  // Rcpp::IntegerVector updated(reverse + 1);
  arma::ivec updated(reverse + 1);
  // int updated;
  arma::imat conflicts(vs.n_rows, vs.n_rows, arma::fill::zeros);
  // conflicts.fill(-1);
  int nscores = 0;
  if (delta.size() == vs.n_rows){  // supplied a valid delta
    // updated(0) = -2;  // skip first vs_cache_fill()
    updated.fill(-2);
  } else{
    delta = arma::vec(vs.n_rows, arma::fill::zeros);
    // updated(0) = -1;  // keep first vs_cache_fill()
    updated.fill(-1);
  }
  
  while (true){
    
    // Rcpp::Rcout << "New iteration ==================================================" << std::endl;
    
    vs_cache_fill(dirMat,
                  reference,
                  nodes,
                  vs,
                  delta,
                  ignore,
                  data,
                  nscores,
                  updated,
                  conflicts,
                  reverse,
                  is_discrete,
                  k, 
                  debug);
    
    if (just_delta) break;
    
    // Rcpp::Rcout << "Filled ==================================================" << std::endl;
    
    vs_opt_step(updated,
                conflicts,
                undirMat,
                dirMat,
                reference,
                nodes,
                vs,
                delta,
                ignore,
                reverse,
                maxp, 
                debug);
    
    // Rcpp::Rcout << "Updated ==================================================" << std::endl;
    
    // Rcpp::Rcout << "undirMat = " << Rcpp::sum(undirMat) << std::endl;
    // Rcpp::Rcout << "updated = " << updated << std::endl;
    // Rcpp::Rcout << "reference = " << arma::sum(reference) << std::endl;
    
    if (updated(0) < 0) break;
  }
  
  // dirMat += undirMat;
  
  return Rcpp::List::create(_["graph"] = dirMat, _["nscores"] = nscores, _["delta"] = delta, _["ignore"] = ignore);
}



/*
 * Convert PDAG g to a DAG. Modeled after pcalg::pdag2dag().
 */

// [[Rcpp::export]]
Rcpp::List pdag2dag_cpp(arma::umat g, 
                        SEXP nodes, 
                        bool direct_all = false){
  
  arma::umat graph;
  int npairs = 0;
  bool not_yet = false;
  if (arma::accu(g) == 0){
    graph = g;
  } else{
    
    arma::umat gm = g;
    arma::umat a = g;
    
    bool go_on = true, is_sink = true, adj_check = true;
    int y, z, i;
    
    arma::uvec removed(a.n_rows, arma::fill::zeros);
    arma::uvec x_vec(1, arma::fill::zeros); 
    
    while (go_on && 
           // (a.size() > 0) && 
           (arma::accu(a) > 0)){
      
      not_yet = true;
      
      for (int x = 0; x < a.n_rows; x++){
        
        if (removed(x) > 0)
          continue;
        
        // check if is a sink
        is_sink = true;
        for (int y = 0; y < a.n_cols; y++){
          
          // if directed edge outwards, not a sink
          if ((a(x, y) == 1) && (a(y, x) == 0)){
            is_sink = false;
            break;
          }
        }
        if (!is_sink)
          continue;
        
        // Rcpp::Rcout << x << " is a sink" << std::endl;
        
        // check if all adjacent vertices are adjacent to each other
        arma::uvec undirected = (a.row(x).t() == 1);
        arma::uvec adj = arma::find((a.col(x) == 1) + undirected);
        
        // if (x == 47) 
        //   adj.t().print("adj = ");
        
        adj_check = true;
        i = 0;
        while (i < adj.size() && adj_check){
          y = adj(i);
          
          // if (x == 47) 
          //   Rcpp::Rcout << "y = " << y << std::endl;
          
          // if undirected
          if ((a(x, y) == 1) && (a(y, x) == 1)){
            
            // check if adjacent to others
            for (int j = 0; j < adj.size(); j++){
              z = adj(j);
              
              if (y != z &&
                  (a(y, z) == 0) && 
                  (a(z, y) == 0)){
                adj_check = false;
                break;
              }
            }
          }
          i++;
        }
        if (!adj_check)
          continue;
        not_yet = false;
        
        // Rcpp::Rcout << x << " passes adjacency check" << std::endl;
        
        // orient edges
        arma::uvec into_x = arma::find(undirected);
        if (into_x.size() > 0){
          x_vec(0) = x;
          gm(x_vec, into_x).zeros();
          gm(into_x, x_vec).ones();
          // Rcpp::Rcout << "arma::accu(gm) = " << arma::accu(gm) << std::endl;
          // Rcpp::Rcout << "Oriented " << into_x.size() << " edges into " << x << std::endl;
        }
        a.row(x).zeros();
        a.col(x).zeros();
        removed(x) = 1L;
      }
      go_on = !not_yet;
    }
    
    arma::uvec undir = arma::find((gm == gm.t()) % (gm == 1));
    gm(undir).zeros();
    npairs = arma::accu(gm) + undir.size() / 2;
    
    if (direct_all){
      arma::uvec ij(2);
      int n = gm.n_cols;
      for (int i = 0; i < undir.size(); i++){
        ij(0) = undir(i) % n;
        ij(1) = (undir(i) - ij(0)) / n;
        // if (arma::sum(gm.col(ij(0)) > gm.col(ij(1))))
        //   ij = arma::reverse(ij);
        // ij.t().print("ij = ");
        if (gm(ij(1), ij(0)) > 0)
          continue;
        bool path_ji = cpp_has_path(ij(1), ij(0), as<Rcpp::IntegerMatrix>(Rcpp::wrap(gm)), nodes);
        if (path_ji){
          // Rcpp::Rcout << path_ji << " path from " << ij(1) << " -> " << ij(0) << std::endl;
          continue;
        }
        gm(ij(0), ij(1)) = 1L;  // direct the edge
        // Rcpp::Rcout << "Directing edge " << ij(0) << " -> " << ij(1) << std::endl;
      }
      
      // Rcpp::Rcout << "arma::accu(gm2) = " << arma::accu(gm2) << std::endl;
      // Rcpp::Rcout << "npairs = " << npairs << std::endl;
      // int ndelete = npairs - arma::accu(gm2);
      // Rcpp::Rcout << arma::accu(gm) << " edges remaining (deleted " << npairs - arma::accu(gm) << " undirected edges)\n";
    }
    npairs = undir.size();
    
    graph = gm;
  }
  return Rcpp::List::create(_["graph"] = graph, _["success"] = !not_yet, _["undirected"] = npairs);
}



/*
 * Check which of Meek's rules R1, R2, R3, and R4 \
 * compels y -> x in a, if any.
 */

int check_R(int y,
            int x,
            arma::umat& a,
            Rcpp::IntegerMatrix& d,
            Rcpp::CharacterVector& nodes,
            arma::mat& cache,
            int R1 = true, int R2 = true, int R3 = true, int R4 = true,
            bool verbose = false){
  
  int R = 0;
  arma::uvec Z3;  // for R3
  for (int z = 0; z < a.n_cols; z++){
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // R1: search for a structure z -> y -- x with z -/- x
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    if (R1){
      // if (verbose)
      //   Rcpp::Rcout << "R1: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
      
      // if z -> y and z -/- x
      if (a(z, y) == 1 && a(y, z) == 0 &&  // if z -> y
          a(z, x) == 0 && a(x, z) == 0){  // and z -/- x
        
        R = 1;
        break;  // out of for (z)
      }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // R2: orientation compelled by cycle
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    if (R2){
      // if (verbose)
      //   Rcpp::Rcout << "R2: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
      
      // if x -> y induces a cycle but y -> x doesn't
      if (cache(x, y) == -1 ||  // detected cycle induced if x -> y
          cpp_has_path(y, x, d, nodes)){  // path y ---> x means cycle induced if x -> y
        // (a(y, z) == 1 && a(z, y) == 0 &&  // y -> z
        // a(z, x) == 1 && a(x, z) == 0)){  // z -> x
        
        R = 2;
        break;  // out of for (z)
      } else R2 = false;  // don't need to check this for any other z
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // R3
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    if (R3){
      // if (verbose)
      //   Rcpp::Rcout << "R1: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
      
      // if y -- z -> x
      if (a(z, x) == 1 && a(x, z) == 0 &&  // z -> x
          a(z, y) == 1 && a(y, z) == 1){  // y -- z
        
        // search if adjacent to any previous Z3
        if (Z3.size() >= 1){
          for (int z3 = 0; z3 < Z3.size(); z3++){  // for each previous y -- z3 -> x
            if (a(z, Z3(z3)) == 0 && a(Z3(z3), z) == 0){  // unshielded z -/- z3
              
              R = 3;
              break;  // out of for (z3)
            }
          }
          if (R) break;  // out of for (z)
        }
        Z3.insert_rows(Z3.n_rows, z);
      }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // R4
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    if (R4){
      // if (verbose)
      //   Rcpp::Rcout << "R4: (x, y, z) = (" << x << ", " << y << ", " << z << ")" << std::endl;
      
      // if y -- z -> x  (same as R3)
      if (a(z, x) == 1 && a(x, z) == 0 &&  // z -> x
          a(z, y) == 1 && a(y, z) == 1){  // z -- y
        
        for (int z4 = 0; z4 < a.n_cols; z4++){
          
          // y -- z4 -> z
          if (a(z4, z) == 1 && a(z, z4) == 0 &&  // z4 -> z
              a(y, z4) == 0 && a(z4, y) == 0){  // y -- z4
            
            R = 4;
            break;  // out of for (z4)
          }
        }
        if (R) break;  // out of for (z)
      }
    }
  }  // end for (z)
  
  // return result
  return R;
}



/*
 * Greedy extension to a DAG in HGI.
 */

// [[Rcpp::export]]
arma::umat pdag2dag_greedy(arma::umat& a,
                           Rcpp::IntegerMatrix& d,
                           arma::vec& reference,
                           Rcpp::CharacterVector& nodes,
                           SEXP data,
                           arma::vec& nscores,
                           int is_discrete = 1,
                           double k = 0,
                           int maxp = 8,
                           bool verbose = false){
  
  if (arma::accu(a) == 0)
    return a;
  
  constexpr double DOUBLE_MIN = 0;  // std::numeric_limits<double>::lowest();
  // constexpr double DOUBLE_MIN = std::numeric_limits<double>::lowest();
  const double tol = arma::datum::eps * 2.0;
  
  bool go_on = true, is_sink = true, adj_check = true, is_sink_adj = true;
  bool reset = false, restart = false;
  int y, z, i, nu, upos, best, from, to, iter = 1;
  double temp_score, best_delta;
  bool extend = true, R1 = true, R2 = true, R3 = true, R4 = true, detected_R3 = false, del = true;
  int R = 0;
  bool rules = false;
  
  arma::ivec del_xy(2);
  // arma::uvec meek_xy(2);
  // meek_xy.fill(-1);
  
  Rcpp::IntegerVector t(1);
  t(0) = -1;
  Rcpp::IntegerVector indices = Rcpp::seq(0, nodes.size()-1);
  Rcpp::LogicalVector bool_parents;
  
  arma::uvec removed(a.n_rows, arma::fill::zeros);
  arma::uvec adj, Z3;
  arma::uvec sinks(a.n_cols, arma::fill::zeros);
  // arma::ivec adj(a.n_rows);
  // adj.fill(-1);
  // arma::ivec Z3(a.n_rows);
  // Z3.fill(-1);
  
  // greedy
  arma::mat cache(a.n_rows, a.n_cols, arma::fill::zeros);
  arma::vec best_op(3, arma::fill::zeros);
  bool bool_cache = true;
  
  // while
  while (go_on && 
         arma::accu(a) > 0){
    
    go_on = false;
    if (verbose){
      Rcpp::Rcout << "==================================================" << std::endl;
      Rcpp::Rcout << "Starting iteration " << iter << ": " << arma::sum(reference) << std::endl;
      Rcpp::Rcout << "==================================================" << std::endl;
    }
    
    // if (iter == 55) break;
    
    bool_cache = true;
    del = false;
    rules = false;
    while (bool_cache || del){
      
      bool_cache = false;
      if (verbose){
        if (del){
          Rcpp::Rcout << "Deleting" << std::endl;
        } else{
          Rcpp::Rcout << "Scoring" << std::endl;
        }
      }
      
      // for each node x
      for (int x = 0; x < a.n_rows; x++){
        
        if (removed(x) > 0) continue;
        if (sinks(x) > 0){
          // Rcpp::Rcout << x << " is already a sink" << std::endl;
          continue;
        }
        
        is_sink_adj = extend;
        if (adj.size() > 0){
          adj.reset();
        }
        nu = 0;
        
        // for each node y
        for (int y = 0; y < a.n_cols; y++){
          
          if (removed(y) > 0 ||
              y == x) continue;
          
          ////////////////////////////////////////////////////////////////////////////////////////////////////
          // sink and adjacencies
          ////////////////////////////////////////////////////////////////////////////////////////////////////
          if (is_sink_adj){
            
            // if directed edge outwards, not a sink
            if (a(x, y) == 1 && a(y, x) == 0){  // x -> y outward
              is_sink_adj = false;
            } else if (a(y, x) == 1){  // y -> x or y -- x
              
              // Rcpp::Rcout << "checking " << x << ": adjacent to " << y << std::endl;
              
              if (a(x, y) == 1)  // y -- x
                nu++;  // increment number of undirected edges
              
              // regardless of orientation, y must be adjacent to all other undirected edges z
              for (const auto &z:adj){
                
                // if y -/- z and y -- x OR z -- x
                if (((a(y, x) == 1 && a(x, y) == 1) ||  // y -- x
                    (a(z, x) == 1 && a(x, z) == 1)) &&  // z -- x (either edges undirected)
                    a(z, y) == 0 && a(y, z) == 0){  // y -/- z
                  is_sink_adj = false;
                  // Rcpp::Rcout << "checking " << x << ": " << y << " not adjacent to " << z << std::endl;
                  break;
                } else{
                  // Rcpp::Rcout << "checking " << x << ": " << y << " is adjacent to " << z << std::endl;
                }
              }
              
              // if still potentially a sink satisfying adjacencies, add as an adjacency
              if (is_sink_adj){
                // Rcpp::Rcout << "checking " << x << ": adding " << y << " as an adjacency of " << x << std::endl;
                adj.insert_rows(adj.n_rows, 1);
                adj(adj.n_rows-1) = y;  // y is an adjacency
              }
            }
          }
          
          if (rules){
            
            // if not undirected or already scored and not better
            // if not scored, score it
            if (!(a(y, x) == 1 && a(x, y) == 1)) continue;  // not undirected y -- x
            //   ||  // not undirected y -- x
            // (cache(y, x) != 0 &&  // already scored and
            // (cache(y, x) - best_op(0)) < tol)) continue;  // not better
            
            // flag if orienting y -> x induces a cycle
            if (cpp_has_path(x, y, d, nodes)){  // path x ---> y means y -> x induces a cycle
              cache(y, x) = -1;  // don't delete from a yet because x -> y might be a thing, which will be checked by R2
              continue;
            }
            
            if (verbose)
              Rcpp::Rcout << "Checking rules for (x, y) = (" << x << ", " << y << 
                ") with score (" << cache(y, x) << ")" << std::endl;
            
            R = check_R(y, x, a, d, nodes, cache, 
                        true, true, true, true, verbose);
            
            if (R > 0){
              if (cache(y, x) == 0){
                
                // cache score
                t(0) = x;
                bool_parents = (d(_, x) > 0);
                bool_parents(y) = true;
                // temp_score = cpp_loglik_dnode2(t, indices[bool_parents],
                //                                nodes, data, k, 0);
                temp_score = local_bic_cpp(t, indices[bool_parents],
                                           nodes, data, k, is_discrete, 0);
                nscores(0)++;
                bool_parents(y) = false;  // to preserve for sink scoring
                cache(y, x) = temp_score - reference(x);
              }
              
              if (verbose)
                Rcpp::Rcout << "R" << R << ": (x, y) = (" << x << ", " << y << 
                  ") with score (" << cache(y, x) << ")" << std::endl;
              
              // determine improvement or not
              if (cache(y, x) < tol){
                cache(y, x) = cache(y, x) - 2;
              }
              if ((cache(y, x) - best_op(0)) > tol){
                best_op(0) = cache(y, x);
                best_op(1) = y;
                best_op(2) = x;
              } else if (best_op(0) <= 0 &&  // non-positive improvement
                //          (cache(y, x) + cache(x, y) - best_op(0)) < tol){  // worse
                // best_op(0) = cache(y, x) + cache(x, y);
                (cache(y, x) - best_op(0)) < tol){  // worse
                best_op(0) = cache(y, x);
                best_op(1) = y;
                best_op(2) = x;
              }
            }
            R = 0;
          }  // end if (rules)
        }  // end for (y)
        
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        // continue investigating sink x
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        if (!is_sink_adj || del)
          continue;
        // Rcpp::Rcout << x << " is a sink" << std::endl;
        
        // remove sink with no undirected edges
        if (nu == 0){
          if (verbose)
            Rcpp::Rcout << "removing node " << x << std::endl;
          a.row(x).zeros();
          a.col(x).zeros();
          removed(x) = 1;  // indicate removed
          bool_cache = true;  // keep scoring because might find another sink
          go_on = true;  // keep going because might find another sink
          continue;
        }
        
        if (t(0) != x){
          t(0) = x;
          bool_parents = (d(_, x) > 0);
        }
        for (const auto &y:adj){
          // for (int y = 0; y < undirected.size(); y++){
          
          if (a(x, y) == 1){  // y -- x since y in adj already has y -> x
            if (cache(y, x) == 0){
              
              // score for adding y -> x
              bool_parents(y) = true;
              // temp_score = cpp_loglik_dnode2(t, indices[bool_parents],
              //                                nodes, data, k, 0);
              temp_score = local_bic_cpp(t, indices[bool_parents],
                                         nodes, data, k, is_discrete, 0);
              nscores(0)++;
              bool_parents(y) = false;
              cache(y, x) = temp_score - reference(x);
            }
            
            // determine improvement or not
            if (cache(y, x) < tol){
              cache(y, x) = cache(y, x) - 2;
            }
            if ((cache(y, x) - best_op(0)) > tol){
              best_op(0) = cache(y, x);
              best_op(1) = y;
              best_op(2) = x;
            } else if (best_op(0) <= 0 &&  // non-positive improvement
              //          (cache(y, x) + cache(x, y) - best_op(0)) < tol){  // worse
              // best_op(0) = cache(y, x) + cache(x, y);
              (cache(y, x) - best_op(0)) < tol){  // worse
              best_op(0) = cache(y, x);
              best_op(1) = y;
              best_op(2) = x;
            }
          }
        }
      }  // end for (x)
      
      // turn on rules if fail to find any sink
      if (!bool_cache && best_op(0) == 0){
        if (!rules){
          bool_cache = true;
          rules = true;
          // test
          // go_on = false;
          // break;
        } else{
          break;
        }
      }
    }  // end while (bool_cache)
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // apply best_op
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    if (best_op(0) > tol){
      if (verbose){
        Rcpp::Rcout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
        Rcpp::Rcout << "orienting " << best_op(1) << " -> " << best_op(2) << " (" << best_op(0) << ")" << std::endl;
        Rcpp::Rcout << "++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      }
      d(best_op(1), best_op(2)) = 1L;
      a(best_op(1), best_op(2)) = 1L;
      a(best_op(2), best_op(1)) = 0L;
      reference(best_op(2)) += best_op(0);
      cache.col(best_op(2)).zeros();  // all nodes -> best_op(2) need to be re-scored
      sinks(best_op(1)) = 0L;  // outgoing edge, so not a sink
      best_op.zeros();
      go_on = true;
    } else if (best_op(0) < 0){
      if (verbose){
        Rcpp::Rcout << "--------------------------------------------------" << std::endl;
        Rcpp::Rcout << "deleting " << best_op(1) << " -> " << best_op(2) << " (" << best_op(0) << ")" << std::endl;
        Rcpp::Rcout << "--------------------------------------------------" << std::endl;
      }
      a(best_op(1), best_op(2)) = 0L;
      a(best_op(2), best_op(1)) = 0L;
      best_op.zeros();
      go_on = true;
    }
    iter++;
  }  // end while (go_on)
  if (verbose)
    Rcpp::Rcout << "arma::accu(a % a.t()) / 2 = " << (arma::accu(a % a.t()) / 2) << std::endl;
  return a;
}



// [[Rcpp::export]]
Rcpp::IntegerMatrix apply_cpdag_rules(arma::umat& pdag,
                                      Rcpp::CharacterVector& nodes,  // bool orient_conflict = true,
                                      bool remove_invalid = true,
                                      bool debug = false){
  
  arma::mat cache(pdag.n_rows, pdag.n_rows);
  Rcpp::IntegerMatrix dag = as<Rcpp::IntegerMatrix>(wrap(pdag - pdag % pdag.t()));
  
  int Rxy = 0;
  int Ryx = 0;
  bool go_on = true;
  
  while (go_on){
    
    if (debug)
      Rcpp::Rcout << "* scanning node pairs\n";
    
    go_on = false;
    
    // for all distinct node pairs
    for (int x = 0; x < pdag.n_rows-1; x++){
      
      for (int y = x+1; y < pdag.n_rows; y++){
        
        // investigate undirected edges between distinct nodes
        if (x == y || cache(x, y) == 1 ||
            pdag(x, y) == 0 || pdag(y, x) == 0) continue;
        
        // check x -> y
        Rxy = check_R(x, y, pdag, dag, nodes, cache, 
                      true, true, true, true, debug);
        
        // check y -> x
        Ryx = check_R(y, x, pdag, dag, nodes, cache, 
                      true, true, true, true, debug);
        
        if (Rxy > 0 && Ryx > 0){
          
          // remove invalid arc
          if (remove_invalid && Rxy == 2 && Ryx == 2){
            
            if (debug)
              Rcpp::Rcout << "  > removing invalid arc " << 
                nodes(x) << " -- " << nodes(y) << "\n";
            
            pdag(x, y) = pdag(y, x) = 0L;
          }
          else{
          
            if (debug)
              Rcpp::Rcout << "  > ignoring " << nodes(x) << " -- " << nodes(y) << 
                " since " << nodes(x) << " -> " << nodes(y) << " compelled by R" << Rxy << 
                  " and " << nodes(y) << " -> " << nodes(x) << " compelled by R" << Ryx << "\n";
            
            cache(x, y) = 1;
          }
        } 
        else if (Rxy > 0){
          
          if (debug)
            Rcpp::Rcout << "  > orienting " << nodes(x) << " -> " << nodes(y) << 
              " compelled by R" << Rxy << "\n";
          
          // apply x -> y
          dag(x, y) = pdag(x, y) = 1L;
          pdag(y, x) = 0L;
          go_on = true;
        }
        else if (Ryx > 0){
          
          if (debug)
            Rcpp::Rcout << "  > orienting " << nodes(y) << " -> " << nodes(x) << 
              " compelled by R" << Ryx << "\n";
          
          // apply y -> x
          dag(y, x) = pdag(y, x) = 1L;
          pdag(x, y) = 0L;
          go_on = true;
        }
        // reset for next node pair
        Rxy = Ryx = 0;
      }
    }
  }
  if (debug)
    Rcpp::Rcout << "* the graph is unchanged, stopping.\n";
  
  return as<Rcpp::IntegerMatrix>(wrap(pdag));
}