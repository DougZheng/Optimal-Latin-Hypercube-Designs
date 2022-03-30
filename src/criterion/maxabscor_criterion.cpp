#include "maxabscor_criterion.h"

#include <cmath>
#include <algorithm>

namespace LHD {
MaxAbsCorCriterion::MaxAbsCorCriterion(const std::vector<std::vector<int>>* design)
   : Criterion(design) {
  const auto& a = *design;
  kCorDenominator = 1.0 * n_run_ * (n_run_ + 1) * (2 * n_run_ + 1) / 6 -
    1.0 * n_run_ * (n_run_ + 1) * (n_run_ + 1) / 4;
  cor_.resize(k_var_, std::vector<double>(k_var_));
  rho_max_ = -1;
  rho_max_i_ = rho_max_j_ = 0;
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
      if (std::fabs(rho) > rho_max_) {
        rho_max_ = std::fabs(rho);
        rho_max_i_ = i;
        rho_max_j_ = j;
      }
    }
  }
}

void MaxAbsCorCriterion::SwapInCol(int col, int r1, int r2) {
  const auto& a = *design_;
  double old_rho_max = rho_max_;
  int old_rho_max_i = rho_max_i_;
  int old_rho_max_j = rho_max_j_;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    int deta = (a[r1][i] - a[r2][i]) *
      (a[r2][col] - a[r1][col]);
    cor_[i][col] = cor_[col][i] = 
      (cor_[i][col] * kCorDenominator + deta) / kCorDenominator;
    if (std::fabs(cor_[i][col]) > rho_max_) {
      rho_max_ = std::fabs(cor_[i][col]);
      rho_max_i_ = i;
      rho_max_j_ = col;
    }
  }
  if ((old_rho_max_i == col || old_rho_max_j == col) &&
       rho_max_ == old_rho_max) {
    rho_max_ = -1;
    for (int i = 0; i < k_var_; ++i) {
      for (int j = 0; j < i; ++j) {
        if (std::fabs(cor_[i][j]) > rho_max_) {
          rho_max_ = std::fabs(cor_[i][j]);
          rho_max_i_ = i;
          rho_max_j_ = j;
        }
      }
    }
  }
}

double MaxAbsCorCriterion::PreSwapInCol(int col, int r1, int r2) const {
  const auto& a = *design_;
  double old_rho_max = rho_max_;
  int old_rho_max_i = rho_max_i_;
  int old_rho_max_j = rho_max_j_;
  double new_rho_max = -1;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    int deta = (a[r1][i] - a[r2][i]) *
      (a[r2][col] - a[r1][col]);
    double tmp_cor = (cor_[i][col] * kCorDenominator + deta) / kCorDenominator;
    if (std::fabs(tmp_cor) > new_rho_max) {
      new_rho_max = std::fabs(tmp_cor);
    }
  }
  if ((old_rho_max_i == col || old_rho_max_j == col) &&
       new_rho_max < old_rho_max) {
    for (int i = 0; i < k_var_; ++i) {
      if (i == col) continue;
      for (int j = 0; j < i; ++j) {
        if (j == col) continue;
        new_rho_max = std::max(new_rho_max, std::fabs(cor_[i][j]));
      }
    }
    return new_rho_max;
  } else {
    return std::max(new_rho_max, rho_max_);
  }
}
} // namespace LHD