#include "avgabscor_criterion.h"

#include <cmath>
#include <algorithm>

namespace LHD {
AvgAbsCorCriterion::AvgAbsCorCriterion(const std::vector<std::vector<int>>* design)
   : Criterion(design) {
  const auto& a = *design;
  kCorDenominator = 1.0 * n_run_ * (n_run_ + 1) * (2 * n_run_ + 1) / 6 -
    1.0 * n_run_ * (n_run_ + 1) * (n_run_ + 1) / 4;
	kCorPairNum = 1.0 * k_var_ * (k_var_ - 1) / 2;
  cor_.resize(k_var_, std::vector<double>(k_var_));
  rho_avg_ = 0;
  for (int i = 0; i < k_var_; ++i) {
    cor_[i][i] = 1;
    for (int j = 0; j < i; ++j) {
      double rho = 0;
      for (int o = 0; o < n_run_; ++o) {
        rho += 1.0 * a[o][i] * a[o][j];
      }
      rho -= 1.0 * n_run_ * (n_run_ + 1) * (n_run_ + 1) / 4;
      rho /= kCorDenominator;
      cor_[i][j] = cor_[j][i] = rho;
			rho_avg_ += std::fabs(rho);
    }
  }
	rho_avg_ /= kCorPairNum;
}

void AvgAbsCorCriterion::SwapInCol(int col, int r1, int r2) {
  const auto& a = *design_;
	rho_avg_ *= kCorPairNum;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    int deta = (a[r1][i] - a[r2][i]) *
      (a[r2][col] - a[r1][col]);
		rho_avg_ -= std::fabs(cor_[i][col]);
    cor_[i][col] = cor_[col][i] = 
      (cor_[i][col] * kCorDenominator + deta) / kCorDenominator;
		rho_avg_ += std::fabs(cor_[i][col]);
  }
	rho_avg_ /= kCorPairNum;
}

double AvgAbsCorCriterion::PreSwapInCol(int col, int r1, int r2) const {
  const auto& a = *design_;
  double ret = rho_avg_ * kCorPairNum;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
		ret -= std::fabs(cor_[i][col]);
    int deta = (a[r1][i] - a[r2][i]) *
      (a[r2][col] - a[r1][col]);
    double tmp_cor = (cor_[i][col] * kCorDenominator + deta) / kCorDenominator;
		ret += std::fabs(tmp_cor);
  }
	return ret / kCorPairNum;
}
} // namespace LHD