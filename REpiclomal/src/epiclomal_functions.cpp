#include <Rcpp.h>
using namespace Rcpp;

// Average distance between methylation of two cells
double dist_pair(NumericVector v1, NumericVector v2) {
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
NumericMatrix dist_PBAL(NumericMatrix d) {
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
NumericMatrix impute_means(NumericMatrix input_data) {
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
        ## MA: added std:: below, without it it may give errors on some platforms
            if (std::isnan(mean)) imputed_data(*idx, i) = NA_REAL;
            else imputed_data(*idx, i) = mean;
        }
    }
    return(imputed_data);
}

// Fill in NA values with average (median) for the loci or cell
// [[Rcpp::export]]
IntegerMatrix impute_medians(IntegerMatrix input_data) {
    IntegerMatrix imputed_data(clone(input_data));
    for (int i = 0; i < input_data.ncol(); i++){
        std::vector<int> na_idx;
        int n_0 = 0;
        int n_1 = 0;
        for (int j = 0; j < input_data.nrow(); j++){
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

NumericVector find_cluster_vector(IntegerVector cluster_set, DataFrame estimate_epi, int rstart, int rend, int cid) {
    NumericVector vector;
    int cluster_id = cluster_set[cid];
    int cluster_idx;
    NumericVector epi_cluster = estimate_epi["cluster_id"];
    for (int idx = 0; idx < estimate_epi.nrows(); idx++) {
        if (epi_cluster[idx] == cluster_id) {
            cluster_idx = idx;
            break;
        }
    }
    for (int col = rstart; col < rend + 1; col ++) {
        NumericVector loci_col = estimate_epi[col];
        vector.push_back(loci_col[cluster_idx]);
    }

    return(vector);
}

// [[Rcpp::export]]
bool find_variable_region(int r, IntegerVector cluster_set, DataFrame estimate_epi, int rstart, int rend) {
    bool variable_region = false;
    for (int cid1 = 0; cid1 < cluster_set.size() - 1; cid1++) {
        for (int cid2 = cid1 + 1; cid2 < cluster_set.size(); cid2++) {
            NumericVector vector1 = find_cluster_vector(cluster_set, estimate_epi, rstart, rend, cid1);
            NumericVector vector2 = find_cluster_vector(cluster_set, estimate_epi, rstart, rend, cid2);

            float diff_sum = 0;
            for (int i = 0; i < vector1.length(); i++) {
                diff_sum += std::abs(vector1[i]-vector2[i]);
            }
            float di = (diff_sum)/vector1.length();

            std::cout << "Region " << r << " cluster1 " << cluster_set[cid1] << " cluster2 " << cluster_set[cid2] << " di " << di << "\n";

            if (di > 0.5) { // large distance
                variable_region = true;
                return(variable_region);
            }
        }
    }
    return(variable_region);
}

// [[Rcpp::export]]
void impute_missing_data(IntegerMatrix meth_data, IntegerMatrix data_estimate_corr, DataFrame estimate_Z, int rstart, int rend, bool variable_region) {
    NumericVector estimate_cluster = estimate_Z[1];
    NumericVector clusters = unique(estimate_cluster);
    std::sort(clusters.begin(), clusters.end());
    for (int j = rstart - 1; j < rend; j++) {
        // traverse by cluster
        for (int c_idx = 0; c_idx < clusters.size(); c_idx++) {
            int c = clusters[c_idx];
            std::vector<int> na_idx;
            int n_0 = 0;
            int n_1 = 0;
            int median = NA_INTEGER;
            for (int cell = 0; cell < estimate_Z.nrow(); cell++) {
                if (estimate_cluster[cell] != c) continue;
                int val = meth_data(cell, j);
                if (IntegerVector::is_na(val)) {
                    na_idx.push_back(cell);
                    continue;
                }
                if (val == 0) n_0++;
                else if (val == 1) n_1++;
            }
            if (n_0 == 0 && n_1 == 0) { // if none of the cells in this cluster had an observed value at this CpG
                // replace the naive matrix with the median of all cells
                for (int cell = 0; cell < meth_data.nrow(); cell++) {
                    int val = meth_data(cell, j);
                    if (val == 0) {
                        n_0++;
                    }
                    else if (val == 1) n_1++;
                }
                if (n_1 > n_0) median = 1;
                else median = 0;
                // also correct the estimated matrix only in this case and if the region is not variable
                // may have to check that this is a region that is mostly similar across clusters
                if(!variable_region) {
                    for (std::vector<int>::iterator idx = na_idx.begin(); idx != na_idx.end(); idx++) {
                        data_estimate_corr(*idx, j) = median;
                    }
                }
            }
            else {
                if (n_1 > n_0) median = 1;
                else median = 0;
            }
            for (std::vector<int>::iterator idx = na_idx.begin(); idx != na_idx.end(); idx++) {
                meth_data(*idx, j) = median;
            }
        }
    }
}
