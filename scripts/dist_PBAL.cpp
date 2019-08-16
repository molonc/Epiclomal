#include <Rcpp.h>
using namespace Rcpp;

double dist_pair(NumericVector v1, NumericVector v2){
    int len = v1.size();
    double l = 0;
    double diff_sum = 0;
    for(int i = 0; i < len; i++){
        int v1_i = v1[i];
        int v2_i = v2[i];
        if(IntegerVector::is_na(v1_i) || IntegerVector::is_na(v2_i)) continue;
        diff_sum += std::abs(v1_i - v2_i);
        l++;
    }
    return(diff_sum / l);
}

// [[Rcpp::export]]
NumericMatrix dist_PBAL(NumericMatrix d){
    int rows = d.nrow();
    NumericMatrix dist_data(rows, rows);
    rownames(dist_data) = rownames(d);
    colnames(dist_data) = rownames(d);
    for(int i = 0; i < rows; i++){
        for(int j = i; j < rows; j++){
            double dist = dist_pair(d.row(i), d.row(j));
            dist_data(j,i) = dist;
        }
    }
    return(dist_data);
}