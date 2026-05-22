// #include <Rcpp.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double sum_nan_cpp(const arma::mat& M){
  arma::mat clone = M;
  clone.replace(NA_REAL, 0);
  return  arma::accu(clone);
}


// [[Rcpp::export]]
arma::mat mat_prod(arma::mat& A, arma::mat& B){
  // % corresponds to element-wise product
  return A%B;
}


// https://stats.stackexchange.com/questions/494701/how-to-convert-a-list-of-integers-into-a-probability-distribution-such-that-the
// to keep this in mind, might be useful
// function that given the convolved neighbourhoods (that tells you the density per cell of the matrix)
// sample one cell, I want sample both x and y, but just the index
// as matrices are actually column vectors
// [[Rcpp::export]]
int sample_next_xy(arma::mat& M){
  arma::uvec v = arma::find_finite(M);
  // M.elem(v) corresponds to the weights used to sample
  // as the spawn matrix contains these
  arma::uvec i = RcppArmadillo::sample(v, 1, false, M.elem(v));
  return i[0];
}

// [[Rcpp::export]]
arma::mat select_row(arma::mat& M, arma::uvec& v){
  arma::mat r = M.rows(v);
  // NumericMatrix res = wrap(r);
  // return res;
  return r;
}

// verison of extract neighbourhood in cpp
// slightly different and simplified than its R version
// [[Rcpp::export]]
arma::mat extract_neighbourhood_cpp(arma::mat& map_mat, const int& x, const int& y, const int& radius){
  // -1 because 0-indexed
  int xl = x-radius-1;
  int xr = x+radius-1;
  int yl = y-radius-1;
  int yr = y+radius-1;
  // check whether we go out of the border
  if(xl<0) xl = 0 ;
  if(yl<0) yl = 0 ;
  if(xr>(map_mat.n_rows-1)) xr = map_mat.n_rows-1;
  if(yr>(map_mat.n_cols-1)) yr = map_mat.n_cols-1;
  return map_mat.submat(xl, yl, xr, yr);
}

// function to obtain a spawning map of the area of interest given the livability map
// and the convolved map
// [[Rcpp::export]]
arma::mat define_spawn_map(arma::mat& conv_map_mat, arma::mat& live_map_mat, const int& x, const int& y, const int& radius){
  // -1 because 0-indexed
  int xl = x-radius-1;
  int xr = x+radius-1;
  int yl = y-radius-1;
  int yr = y+radius-1;
  // check whether we go out of the border
  if(xl<0) xl = 0 ;
  if(yl<0) yl = 0 ;
  if(xr>(conv_map_mat.n_rows-1)) xr = conv_map_mat.n_rows-1;
  if(yr>(conv_map_mat.n_cols-1)) yr = conv_map_mat.n_cols-1;
  // return the product of the two
  // % operator in armadillo is the element-wise product
  arma::mat spawn_map = 1-(conv_map_mat.submat(xl, yl, xr, yr)%live_map_mat.submat(xl, yl, xr, yr));
  // clamp it betwwen 0 and 1
  spawn_map.clamp(0, 1);
  return spawn_map;
}

// check if neighbourhood has space
// slightly different and simplified than its R version
// [[Rcpp::export]]
bool check_neighbourhood_cpp(arma::mat& M, const double& ratio){
  double num = arma::accu(M==-1);
  arma::uvec n_nan = arma::find_nonfinite(M);
  double den = M.n_rows*M.n_cols-n_nan.n_elem;
  return num/den>ratio;
}

// THIS HAS TO BE THOUGHT IN A BETTER WAY
// check if the spawning map ha free spave to let other populations to spawn
// [[Rcpp::export]]
bool check_spawn_map_cpp(arma::mat& M){
  // check how many elements of M are non 0
  // which means probability of spawning
  // (NA are by default non livable)
  // those are spaces where a new population could spawn
  double count = arma::accu(arma::find(M>0));
  // arma::uvec n_nan = arma::find_nonfinite(M);
  // check if the ratio of available places and total places is above a
  // certain number
  return count>0;
}

// [[Rcpp::export]]
arma::mat make_binary(arma::mat& M, const bool& keep_na){
  arma::mat bin = arma::mat(M.n_rows, M.n_cols, arma::fill::zeros);
  // arma::uvec pop_vec = arma::find(Ma>0);
  bin.elem(arma::find(M>=0)).ones();
  // if we want to keep the NA value (not for convolution as na gives problems)
  if(keep_na) bin.elem(arma::find_nonfinite(M)).operator/=(NA_REAL);
  return bin;
}


// [[Rcpp::export]]
arma::mat convolve_mat_full_cpp(arma::mat& M, arma::mat& kernel){
  // get binary
  arma::mat bin=make_binary(M, false);
  arma::uvec is_na = arma::find_nonfinite(M);
  arma::uvec is_pop = arma::find(bin==1);
  // convolve it
  arma::mat conv_mat = arma::conv2(bin, kernel.t(), "same");
  // set the elements that were NA back to NA
  conv_mat.elem(is_na).operator/=(NA_REAL);
  // set the elements that were occupied by a pop NA
  conv_mat.elem(is_pop).operator/=(NA_REAL);
  return conv_mat;
}

// version where we do not need to convert into binary
// guess maybe is faster
// // [[Rcpp::export]]
// arma::mat convolve_bin_mat_full_cpp(arma::mat& M, arma::mat& bin, arma::mat& kernel){
//   // get binary
//   // arma::mat bin=make_binary(M, false);
//   // convolve it
//   arma::mat conv_mat = arma::conv2(bin, kernel.t(), "same");
//   conv_mat.elem(arma::find_nonfinite(M)).operator/=(NA_REAL);
//   return conv_mat;
// }


// [[Rcpp::export]]
IntegerVector convert_xy_new(const int& x_par, const int& y_par, const int& x_new, const int& y_new, const int& radius){
  // -1 because we want the corner out of the area
  int xl = x_par-radius-1;
  int yl = y_par-radius-1;
  // check whether we go out of the border
  if(xl<0) xl = 0 ;
  if(yl<0) yl = 0 ;
  // now we add the top left corner we obtained to the x_new, y_new values
  // that we obtained
  return IntegerVector::create(xl+x_new, yl+y_new);
}

// [[Rcpp::export]]
void update_conv_map_mat_cpp(arma::mat& conv_map_mat, arma::mat& map_mat, arma::mat& kernel, const int& x, const int& y, const int& radius){
  int radius_db = radius*2;
  arma::mat map_mat_neighb = extract_neighbourhood_cpp(map_mat, x, y, radius_db);
  arma::mat conv_map_mat_neighb = convolve_mat_full_cpp(map_mat_neighb, kernel);
  // get the limit of the map_mat_neigh
  // as it is double the size of what it should be sue to correct convolution
  // NEVER SURE, BUT IT SHOULD BE RIGHT
  int xl_n;
  int xr_n;
  int yl_n;
  int yr_n;
  // already thinking in 0-indexed
  if(x-radius_db>=1){
    xl_n = radius;
    xr_n = radius_db+radius;
    if(xr_n>=map_mat_neighb.n_rows) xr_n=map_mat_neighb.n_rows-1;
  }
  else if(x-radius<1){
    xl_n = 0;
    xr_n = x+radius-1;
    if(xr_n>=map_mat_neighb.n_rows) xr_n=map_mat_neighb.n_rows-1;
  }
  else if(x-radius>=1){
    xl_n = x-radius-1;
    xr_n = x+radius-1;
    if(xr_n>=map_mat_neighb.n_rows) xr_n=map_mat_neighb.n_rows-1;
  }
  if(y-radius_db>=1){
    yl_n = radius;
    yr_n = radius_db+radius;
    if(yr_n>=map_mat_neighb.n_cols) yr_n=map_mat_neighb.n_cols-1;
  }
  else if(y-radius<1){
    yl_n = 0;
    yr_n = y+radius-1;
    if(yr_n>=map_mat_neighb.n_cols) yr_n=map_mat_neighb.n_cols-1;
  }
  else if(y-radius>=1){
    yl_n = y-radius-1;
    yr_n = y+radius-1;
    if(yr_n>=map_mat_neighb.n_cols) yr_n=map_mat_neighb.n_cols-1;
  }

  // get the right indices of the conv_map_mat
  // to update the area
  // -1 because 0-indexed
  int xl = x-radius-1;
  int xr = x+radius-1;
  int yl = y-radius-1;
  int yr = y+radius-1;
  //  check whether we go out of the border
  if(xl<0) xl = 0;
  if(yl<0) yl = 0;
  if(xr>(conv_map_mat.n_rows-1)) xr = conv_map_mat.n_rows-1;
  if(yr>(conv_map_mat.n_cols-1)) yr = conv_map_mat.n_cols-1;
  // update matrix
  conv_map_mat.submat(xl, yl, xr, yr) = conv_map_mat_neighb.submat(xl_n, yl_n, xr_n, yr_n);
  // return conv_map_mat;
}

// if need to be added in the argument, maybe having a bin version fo the matrix can
// save computation
// arma::mat& bin_map_mat,
//
//
// [[Rcpp::export]]
void bd_loop(NumericMatrix& pop_mat, arma::mat& map_mat, arma::mat& conv_map_mat, arma::mat& live_map_mat, arma::mat& kernel, const NumericVector& pop_ids, const NumericVector& sampled_events, const NumericVector& N,const int& radius, NumericVector& next_pop_v, const int& offset, NumericVector& n_v, const int& n_max, const int& t, NumericVector& i_mat_v) {

  int next_pop = next_pop_v[0];
  int i_mat = i_mat_v[0];
  int n = n_v[0];

  for (int i = 0; i < pop_ids.length(); i++){
    int pop_chosen = pop_ids[i];
    // event birth=1, death=2
    int event = sampled_events[i];
    // convert the pop chosen into row index
    // of the dataframe
    int row_id = pop_chosen+offset-1;
    if (n == n_max && event==1) {
      continue;
    }else if(event==1){
      // 2nd col is parent -- 0-indexed
      int next_parent = pop_mat(row_id, 1);
      // 6th col is x -- 0-indexed
      int x_par = pop_mat(row_id, 5);
      // 7th col is y -- 0-indexed
      int y_par = pop_mat(row_id, 6);
      // now we need to check whether the neighbourhood of this population
      // allows for a new population to spawn
      // to do so we use both the convolved map
      // and the liveability map, the resulting map will be a product between the two
      arma::mat spawn_map = define_spawn_map(conv_map_mat, live_map_mat, x_par, y_par, radius);
      if(!check_spawn_map_cpp(spawn_map)) continue;
      // now we need to decide where the next pop will spawn
      // sample one cell as single index, spawn_map is
      // a map of probabilities
      int sampled_ind =  sample_next_xy(spawn_map);
      int nrow =  spawn_map.n_rows;
      // note that the new coordinates are in the 1 based coordinates system
      int x_new_nc = sampled_ind%nrow+1;
      int y_new_nc = sampled_ind/nrow+1;
      // conveert them to their position in the original map
      IntegerVector xy_new = convert_xy_new(x_par, y_par, x_new_nc, y_new_nc, radius);
      int x_new = xy_new[0];
      int y_new = xy_new[1];
      // set values of new matrix and binary matrix
      map_mat(x_new-1, y_new-1) = next_pop;
      // maybe it has some advantage to have a bin mat
      // bin_map_mat(x_new, y_new) = 1;
      // update the convolved matrix
      update_conv_map_mat_cpp(conv_map_mat, map_mat, kernel, x_new, y_new, radius);

      // TO FIND A WAY HOW TO INTRODUCE THE FUNCTION FOR THE POPULATION
      // SIZE
      //
      // create and update row in the matrix
      NumericVector new_row = NumericVector::create(next_parent,
                                                    next_pop,
                                                    t,
                                                    0,
                                                    1,
                                                    x_new,
                                                    y_new,
                                                    N[i],
                                                     1);
      pop_mat(i_mat-1, _) = new_row;
      // increase all the counters
      next_pop+=1;
      n+=1;
      i_mat+=1;
      // dead condition
    }else if (n > 1 && event == 2){
      // 5th col is x -- 0-indexed
      int x_kill = pop_mat(row_id, 5);
      // 6th col is y -- 0-indexed
      int y_kill = pop_mat(row_id, 6);
      // set -1 in map_mat and 0 in binary
      map_mat(x_kill-1, y_kill-1) = -1;
      // maybe it has some advantage to have a bin mat
      // bin_map_mat(x_kill-1, y_kill-1) = 0;
      // update the convolved matrix
      update_conv_map_mat_cpp(conv_map_mat, map_mat, kernel, x_kill, y_kill, radius);
      // update t_end (4th) -- 0-indexed
      pop_mat(row_id, 3) = t;
      // set alive to dead (5th) -- 0-indexed
      pop_mat(row_id, 4) = 0;
      // decrease counter
      n-=1;
    }
  }

  next_pop_v[0] = next_pop;
  i_mat_v[0] = i_mat;
  n_v[0] = n;
  // return ????????;
}


///////////////////////////////////////////////////////////////////////
//                EXTRA                                              //
///////////////////////////////////////////////////////////////////////


// just a wrapper
// [[Rcpp::export]]
arma::mat conv_2d_cpp(const arma::mat& pad, const arma::mat& kernel) {

  arma::mat conv_mat = arma::conv2(pad, kernel, "same");

  return conv_mat;
}


