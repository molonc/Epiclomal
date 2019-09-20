#include <Rcpp.h>
using namespace Rcpp;

// Average distance between methylation of two cells
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

// create distance matrix with PBAL distance
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

// Fill in NA values with average (mean) for the loci or cell
// [[Rcpp::export]]
NumericMatrix impute_means(NumericMatrix input_data){
    NumericMatrix imputed_data(clone(input_data));
    for (int i = 0; i < input_data.ncol(); i++){
        std::vector<int> na_idx;
        float count = 0;
        float total = 0;
        for (int j = 0; j < input_data.nrow(); j++){
            int val = input_data(j, i);
            if (IntegerVector::is_na(val)) {
                na_idx.push_back(j);
                continue;
            }
            total += val;
            count ++;
        }
        float mean = total/count;
        for(std::vector<int>::iterator idx = na_idx.begin(); idx != na_idx.end(); idx ++){
            if (isnan(mean)) imputed_data(*idx, i) = NA_REAL;
            else imputed_data(*idx, i) = mean;
        }
    }
    return(imputed_data);
}

// Fill in NA values with average (median) for the loci or cell
// [[Rcpp::export]]
IntegerMatrix impute_medians(IntegerMatrix input_data){
    IntegerMatrix imputed_data(clone(input_data));
    for (int i = 0; i < input_data.ncol(); i++){
        std::vector<int> na_idx;
        int n_0 = 0;
        int n_1 = 0;
        for (int j = 0; j < input_data.nrow(); j ++){
            int val = input_data(j, i);
            if (IntegerVector::is_na(val)) {
                na_idx.push_back(j);
                continue;
            }
            if (val == 0) n_0++;
            else if (val == 1) n_1++;
        }
        int median;
        if (n_1 > n_0) median = 1;
        else median = 0;
        for(std::vector<int>::iterator idx = na_idx.begin(); idx != na_idx.end(); idx ++){
            if (n_0 + n_1 == 0) imputed_data(*idx, i) = NA_INTEGER;
            else imputed_data(*idx, i) = median;
        }
    }
    return(imputed_data);
}