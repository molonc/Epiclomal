#include <Rcpp.h>
using namespace Rcpp;

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