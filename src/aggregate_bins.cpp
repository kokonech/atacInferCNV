#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix aggregate_bins_cpp_sparse(IntegerVector targ_bins,
                                        IntegerVector query_hits,
                                        IntegerVector subject_hits,
                                        const Eigen::SparseMatrix<double> &signal_matrix,
                                        int num_bins) {
  int ncol_signal = signal_matrix.cols();  // Number of columns in the signal matrix
  NumericMatrix binned_signal(num_bins, ncol_signal);  // Initialize result matrix

  // Iterate over each target bin
  for (int bin : targ_bins) {
    // Find the regions that overlap with this bin
    std::vector<int> overlapping_regions;
    for (int i = 0; i < subject_hits.size(); i++) {
      if (subject_hits[i] == bin) {
        overlapping_regions.push_back(query_hits[i] - 1);  // Convert to 0-based indexing
      }
    }

    // Sum the signal for overlapping regions
    if (!overlapping_regions.empty()) {
      for (int col = 0; col < ncol_signal; col++) {
        double col_sum = 0;
        for (int region : overlapping_regions) {
          col_sum += signal_matrix.coeff(region, col);  // Access sparse matrix elements
        }
        binned_signal(bin - 1, col) = col_sum;  // Store in 0-based bin
      }
    }
  }

  return binned_signal;
}
