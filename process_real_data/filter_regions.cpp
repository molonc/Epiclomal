#include <Rcpp.h>

using namespace Rcpp;

IntegerVector non_na_intersection(NumericVector x, NumericVector y) {
    IntegerVector indices;
    for (int i = 0; i < x.length(); i++) {
        if (!NumericVector::is_na(x[i]) && !NumericVector::is_na(y[i])){
            indices.push_back(i);
        }
    }
    return(indices);
}

double cpp_mean(NumericVector arr, IntegerVector indices) {
    double sum = 0;
    double n = indices.length();

    for (int i = 0; i < n; i++) {
        int idx = indices[i];
        sum += arr[idx];
    }

    if (n == 0) {
        return(NA_REAL);
    }
    else {
        // std::cout << sum << "\t" << n << "\n";
        return(sum / n);
    }
}

double stdDev(NumericVector arr, IntegerVector indices) {
    double sum = 0;
    double n = indices.length();
    double mean = cpp_mean(arr, indices);

    for (int i = 0; i < n; i++) {
        int idx = indices[i];
        sum += pow((arr[idx] - mean), 2);
    }

    if (n == 0) {
        return (NA_REAL);
    }
    else {
        return(sqrt(sum / n));
    }
}

double find_pearson_corr(NumericVector x, NumericVector y) {
    IntegerVector indices = non_na_intersection(x, y);

    double xMean = cpp_mean(x, indices);
    double yMean = cpp_mean(y, indices);

    double xStdDev = stdDev(x, indices);
    double yStdDev = stdDev(y, indices);

    if(NumericVector::is_na(xMean) || NumericVector::is_na(yMean) || NumericVector::is_na(xStdDev) || NumericVector::is_na(yStdDev)) {
        return NA_REAL;
        // std::cout << "something is NA";
    }

    double sum = 0;
    double n = indices.length();

    for (int i = 0; i < n; i++) {
        int idx = indices[i];
        sum += (x[idx] - xMean) * (y[idx] - yMean);
    }

    if(n == 0) {
        return (NA_REAL);
    }
    else {
        // std::cout << " sum " << sum << " n " << n << " xstdDev " << xStdDev << " yStdDev " << yStdDev << "\n";
        return(sum / (n * xStdDev * yStdDev));
    }
}

bool region_already_removed(CharacterVector regions_to_remove, std::string region_name){
    for (int i = 0; i < regions_to_remove.length(); i++) {
        if (as<std::string>(regions_to_remove[i]) == region_name) {
            return true;
        }
    }
    return false;
}

int number_of_na(NumericVector region) {
    int no_na = 0;
    for (int i = 0; i < region.length(); i++) {
        if (NumericVector::is_na(region[i])) {
            no_na++;
        }
    }
    return(no_na);
}

// find regions with high correlation to remove and region weights
// [[Rcpp::export]]
List find_correlated_regions_cpp(NumericMatrix mean_meth_matrix, double coef_t) {
    int cols = mean_meth_matrix.ncol();
    CharacterVector regions = colnames(mean_meth_matrix);
    CharacterVector regions_to_remove;
    IntegerVector region_weights(cols, 1);

    for (int r1 = 0; r1 < (cols - 1); r1++) {
        if (region_already_removed(regions_to_remove, as<std::string>(regions[r1]))) continue;

        for (int r2 = (r1 + 1); r2 < cols; r2++) {
            if (region_already_removed(regions_to_remove, as<std::string>(regions[r2]))) continue;

            double pearson_corr = find_pearson_corr(mean_meth_matrix.column(r1), mean_meth_matrix.column(r2));

            if (pearson_corr == NA_REAL) continue;
            else {
                if((pearson_corr > coef_t) || (pearson_corr < (-1 * coef_t))) {
                    if (number_of_na(mean_meth_matrix.column(r2)) >= number_of_na(mean_meth_matrix.column(r1))) {
                        regions_to_remove.push_back(regions[r2]);
                        region_weights[r1] = region_weights[r1] + region_weights[r2];
                    }
                    else {
                        regions_to_remove.push_back(regions[r1]);
                        region_weights[r2] = region_weights[r1] + region_weights[r2];
                        break;
                    }
                }
            }
        }
    }
    DataFrame region_weights_df = DataFrame::create(Named("region") = colnames(mean_meth_matrix), Named("weight") = region_weights);

    return(List::create(Named("regions_to_remove") = regions_to_remove, Named("region_weights") = region_weights_df));
}

// find regions with high correlation to remove and region weights
// [[Rcpp::export]]
DataFrame find_correlated_regions_LV(NumericMatrix mean_meth_matrix, double coef_t) {
    int cols = mean_meth_matrix.ncol();
    CharacterVector regions = colnames(mean_meth_matrix);
    LogicalVector keep_region(cols, true);
    IntegerVector region_weights(cols, 1);

    for (int r1 = 0; r1 < (cols - 1); r1++) {
        if (!keep_region[r1]) continue;

        for (int r2 = (r1 + 1); r2 < cols; r2++) {
            if (!keep_region[r2]) continue;

            double pearson_corr = find_pearson_corr(mean_meth_matrix.column(r1), mean_meth_matrix.column(r2));

            if (pearson_corr == NA_REAL) continue;
            else {
                if((pearson_corr > coef_t) || (pearson_corr < (-1 * coef_t))) {
                    if (number_of_na(mean_meth_matrix.column(r2)) >= number_of_na(mean_meth_matrix.column(r1))) {
                        // std::cout << "removing region 2: " << regions[r2] << " pearson_corr: " << pearson_corr << " region 1: " << regions[r1];
                        // std::cout << " r1 na: " << number_of_na(mean_meth_matrix.column(r1)) << " r2 na: " << number_of_na(mean_meth_matrix.column(r2)) << "\n";
                        keep_region[r2] = false;
                        region_weights[r1] = region_weights[r1] + region_weights[r2];
                    }
                    else {
                        // std::cout << "removing region 1: " << regions[r1] << " pearson_corr: " << pearson_corr << " region 2: " << regions[r2];
                        // std::cout << " r1 na: " << number_of_na(mean_meth_matrix.column(r1)) << " r2 na: " << number_of_na(mean_meth_matrix.column(r2)) << "\n";
                        keep_region[r1] = false;
                        region_weights[r2] = region_weights[r1] + region_weights[r2];
                        break;
                    }
                }
            }
        }
    }

    DataFrame correlated_regions = DataFrame::create(Named("region") = colnames(mean_meth_matrix), Named("keep_region") = keep_region, Named("weight") = region_weights);

    return(correlated_regions);
}