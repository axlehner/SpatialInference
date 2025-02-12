#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Haversine distance function
double haversine_cpp(double lat1, double lon1, double lat2, double lon2) {
  double R = 6371.0; // Radius of Earth in km
  double dlat = (lat2 - lat1) * M_PI / 180.0;
  double dlon = (lon2 - lon1) * M_PI / 180.0;
  double a = sin(dlat / 2) * sin(dlat / 2) +
             cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
             sin(dlon / 2) * sin(dlon / 2);
  double c = 2 * atan2(sqrt(a), sqrt(1 - a));
  return R * c;
}

// Compute distance matrix with kernel weights
// [[Rcpp::export]]
arma::mat DistMat(arma::mat M, double cutoff, std::string kernel="Bartlett") {
  long long int nrow = M.n_rows;
  arma::mat dmat(nrow, nrow, fill::zeros);

  for (long long int i = 0; i < nrow; i++) {
    dmat(i, i) = 1;
    for (long long int j = i+1; j < nrow; j++) {
      double d = haversine_cpp(M(i,0), M(i,1), M(j,0), M(j,1));
      int v = d <= cutoff;

      if (kernel == "Bartlett") {
        dmat(i, j) = dmat(j, i) = (1 - d / cutoff) * v;
      } else if (kernel == "Epanechnikov") {
        dmat(i, j) = dmat(j, i) = (1 - pow(d / cutoff, 2)) * v;
      } else if (kernel == "Gaussian") {
        dmat(i, j) = dmat(j, i) = exp(-0.5 * pow(d / cutoff, 2)) * v;
      } else if (kernel == "Parzen") {
        if (d <= 0.5 * cutoff) {
          dmat(i, j) = dmat(j, i) = (1 - 6 * pow(d / cutoff, 2) + 6 * pow(d / cutoff, 3)) * v;
        } else if (d <= cutoff) {
          dmat(i, j) = dmat(j, i) = 2 * pow(1 - d / cutoff, 3) * v;
        } else {
          dmat(i, j) = dmat(j, i) = 0;
        }
      } else if (kernel == "Biweight") {
        if (d <= cutoff) {
          dmat(i, j) = dmat(j, i) = (15.0 / 16.0) * pow(1 - pow(d / cutoff, 2), 2) * v;
        } else {
          dmat(i, j) = dmat(j, i) = 0;
        }
      } else {
        dmat(i, j) = dmat(j, i) = v;
      }
    }
  }
  return dmat;
}
