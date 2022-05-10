#include <Rcpp.h>

using namespace Rcpp;

// finds the indices that have values in both Vectors (non NA).
LogicalVector non_na_intersection(NumericVector x, NumericVector y) {
    int n = x.length();
    LogicalVector complete(n, TRUE);
    for (int i = 0; i < n; i++) {
        if (NumericVector::is_na(x[i]) || NumericVector::is_na(y[i])){
            complete[i] = FALSE;
        }
    }
    return(complete);
}

// Calculates the mean methylation value across cells
// Only includes cells that have values in both regions in the calculation
double cpp_mean(NumericVector arr, LogicalVector complete) {
    double sum = 0;
    double n = 0;

    for (int i = 0; i < arr.length(); i++) {
        if(complete[i]) {
            sum += arr[i];
            n++;
        }
    }

    if (n == 0) {
        return(NA_REAL);
    }

    // std::cout << sum << "\t" << n << "\n";
    return(sum / n);
}

// Calculates the standard deviation of methylation across cells
// Only includes cells that have values in both regions in the calculation
double stdDev(NumericVector arr, double mean, LogicalVector complete) {
    double sum = 0;
    double n = 0;

    for (int i = 0; i < arr.length(); i++) {
        if(complete[i]) {
            sum += pow((arr[i] - mean), 2);
            n++;
        }
    }

    if (n == 0) {
        return (NA_REAL);
    }

    return(sqrt(sum / n));
}

// Calculates the Pearson Correlation Coefficient between two Numeric Vectors
// is equivalent to cor(x, y, method = "pearson", use = "na.or.complete") in R
double find_pearson_corr(NumericVector x, NumericVector y) {
    LogicalVector complete = non_na_intersection(x, y);

    double xMean = cpp_mean(x, complete);
    double yMean = cpp_mean(y, complete);

    double xStdDev = stdDev(x, xMean, complete);
    double yStdDev = stdDev(y, yMean, complete);

    // double sum_X = 0, sum_Y = 0, sum_XY = 0;
    // double squareSum_X = 0, squareSum_Y = 0;
    double sum = 0;
    double n = 0;
    for (int i = 0; i < x.length(); i++) {
        if(complete[i]) {
            sum += (x[i] - xMean) * (y[i] - yMean);
            // sum_X += x[i];
            // sum_Y += y[i];
            // sum_XY += x[i] * y[i];
            // squareSum_X += pow(x[i], 2);
            // squareSum_Y += pow(y[i], 2);
            n++;
        }
    }

    if(n == 0) {
        // There are no cells where both regions have data
        return (NA_REAL);
    }
    // std::cout << " sum " << sum << " n " << n << " xstdDev " << xStdDev << " yStdDev " << yStdDev << "\n";
    return(sum / (n * xStdDev * yStdDev));

    // double corr = ((n * sum_XY) - (sum_X * sum_Y)) / sqrt(((n * squareSum_X) - (sum_X * sum_X)) * ((n * squareSum_Y) - (sum_Y * sum_Y)));
    // return(corr);
}

// Counts the number of NA values in a Numeric Vector
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
// If the absolute value of the Pearson Correlation Coefficient between a pair of regions is greater than the threshold,
// the region that has more missing data (NAs) is removed.
// [[Rcpp::export]]
DataFrame find_correlated_regions_cpp(NumericMatrix mean_meth_matrix, double coef_t) {
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