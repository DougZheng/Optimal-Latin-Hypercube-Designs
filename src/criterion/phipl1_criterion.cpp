#include "phipl1_criterion.h"

#include <cmath>
#include <algorithm>

#include "comm/utils.h"

namespace LHD {
PhiPL1Criterion::PhiPL1Criterion(const std::vector<std::vector<int>>* design, int power)
  : Criterion(design),
    kPower(power) {
  const auto& a = *design;
  dis_.resize(n_run_, std::vector<int>(n_run_));
  phi_p_ = 0;
  for (int i = 0; i < n_run_; ++i) {
    dis_[i][i] = 0;
    for (int j = 0; j < i; ++j) {
      int l1_dis = 0;
      for (int o = 0; o < k_var_; ++o) {
        l1_dis += std::abs(a[i][o] - a[j][o]);
      }
      dis_[i][j] = dis_[j][i] = l1_dis;
      phi_p_ += 1.0 / Utils::QuickPow(dis_[i][j], power);
    }
  }
  phi_p_ = std::pow(phi_p_, 1.0 / power);
}

std::pair<double, double> PhiPL1Criterion::GetCriterionBound() const {
  static const auto s_bound = [&]() -> std::pair<double, double> {
    double d_avg = 1.0 * (n_run_ + 1) * k_var_ / 3;
    double d_floor = std::floor(d_avg);
    double d_ceil = std::ceil(d_avg);
    double lower_bound = std::pow(n_run_ * (n_run_ - 1) / 2 * 
      ((d_ceil - d_avg) / std::pow(d_floor, kPower) + 
      (d_avg - d_floor) / std::pow(d_ceil, kPower)), 1.0 / kPower);
    double upper_bound = 0;
    for (int i = 1; i < n_run_; ++i) {
      upper_bound += (n_run_ - i) / std::pow(i * k_var_, kPower);
    }
    upper_bound = std::pow(upper_bound, 1.0 / kPower);
    return {lower_bound, upper_bound};
  }();
  return s_bound;
}

void PhiPL1Criterion::SwapInCol(int col, int r1, int r2) {
  const auto& a = *design_;
  double phi_p_num = Utils::QuickPow(phi_p_, kPower);
  auto Update = [this, &a, &phi_p_num, col](int r1, int r2, int deta) {
    phi_p_num -= 1.0 / Utils::QuickPow(dis_[r1][r2], kPower);
    dis_[r1][r2] -= std::abs(a[r2][col] - a[r1][col]);
    dis_[r1][r2] += std::abs(a[r2][col] - a[r1][col] + deta);
    phi_p_num += 1.0 / Utils::QuickPow(dis_[r1][r2], kPower);
  };
  int deta = a[r2][col] - a[r1][col];
  for (int i = 0; i < n_run_; ++i) {
    if (i == r1 || i == r2) continue;
    i > r1 ? Update(i, r1, deta) : Update(r1, i, -deta);
    i > r2 ? Update(i, r2, -deta) : Update(r2, i, deta);
  }
  phi_p_ = std::pow(phi_p_num, 1.0 / kPower);
}

double PhiPL1Criterion::PreSwapInCol(int col, int r1, int r2) const {
  const auto& a = *design_;
  double phi_p_num = Utils::QuickPow(phi_p_, kPower);
  auto Update = [this, &a, &phi_p_num, col](int r1, int r2, int deta) {
    phi_p_num -= 1.0 / Utils::QuickPow(dis_[r1][r2], kPower);
    int tmp_dis = dis_[r1][r2] - std::abs(a[r2][col] - a[r1][col]) +
      std::abs(a[r2][col] - a[r1][col] + deta);
    phi_p_num += 1.0 / Utils::QuickPow(tmp_dis, kPower);
  };
  int deta = a[r2][col] - a[r1][col];
  for (int i = 0; i < n_run_; ++i) {
    if (i == r1 || i == r2) continue;
    i > r1 ? Update(i, r1, deta) : Update(r1, i, -deta);
    i > r2 ? Update(i, r2, -deta) : Update(r2, i, deta);
  }
  return std::pow(phi_p_num, 1.0 / kPower);
}
} // namespace LHD