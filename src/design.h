#pragma once

#include <cmath>
#include <cassert>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <random>

namespace SearchingAlgorithm {
class Design {
 public:
  using VecInt2D = std::vector<std::vector<int>>;
  using VecDouble2D = std::vector<std::vector<double>>;
  static VecInt2D ReadDesign(const std::string& file);
  Design(const VecInt2D& design);
  Design(int n, int k, int seed);
  inline int GetN() const { return n_run_; }
  inline int GetK() const { return k_var_; }
  inline VecInt2D GetDesign() const { return design_; }
  inline double GetPhiP() const { return update_dis_ ? phi_p_ : -1; };
  inline double GetCritVal(double w, double rho_max, double phi_p) const {
    return w * rho_max + 
      (1 - w) * (phi_p - phi_p_low_) / (phi_p_up_ - phi_p_low_);
  }
  inline double GetPhiPLow() const { return phi_p_low_; }
  inline double GetPhiPUp() const { return phi_p_up_; }
  inline void DisableCorr() { update_corr_ = false; }
  inline void DisableDis() { update_dis_ = false; }
  void EnableCorr();
  void EnableDis();
  void SwapInCol(int col, int r1, int r2);
  std::pair<double, double> PreSwapInCol(int col, int r1, int r2) const;
  double GetMaxAbsCorr(int col) const;
  double GetMaxAbsCorr() const;
  double GetMaxAbsCorrExcept(int col) const;
  void Display();
 private:
  inline double QuickPow15(double x) const {
    double y = x;
    y *= y, y *= y, y *= y, y *= y;
    return y / x;
  }
  void InitCorr();
  void InitDis();
  void InitPhiPBound();
  void MaintainCorr(int col, int r1, int r2);
  void MaintainDis(int col, int r1, int r2);
  double GetPreSwapMaxAbsCorr(int col, int r1, int r2) const;
  double GetPreSwapPhiP(int col, int r1, int r2) const;
 private:
  const int kPInPhi = 15;
  double kCorrDenominator;
  int n_run_;
  int k_var_;
  VecInt2D design_;
  VecDouble2D corr_;
  VecInt2D dis_;
  double phi_p_;
  double phi_p_low_;
  double phi_p_up_;
  bool update_corr_;
  bool update_dis_;
};

Design::VecInt2D Design::ReadDesign(const std::string& file) {
  std::ifstream in(file);
  if (in.fail()) {
    std::cerr << "open " << file << " failed." << std::endl;
    return VecInt2D();
  }
  int n, k;
  in >> n >> k;
  VecInt2D a(n, std::vector<int>(k));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      in >> a[i][j];
    }
  }
  return a;
}

Design::Design(const VecInt2D& design) : design_(design) {
  assert(design.size() > 0);
  n_run_ = design.size();
  k_var_ = design[0].size();
  InitCorr();
  InitDis();
}

Design::Design(int n, int k, int seed) : n_run_(n), k_var_(k) {
  std::mt19937 rng_(seed);
  design_.resize(n, std::vector<int>(k));
  std::vector<int> rnd_list(n);
  std::iota(rnd_list.begin(), rnd_list.end(), 1);
  for (int i = 0; i < k; ++i) {
    std::shuffle(rnd_list.begin(), rnd_list.end(), rng_);
    for (int j = 0; j < n; ++j) {
      design_[j][i] = rnd_list[j];
    }
  }
  InitCorr();
  InitDis();
  update_corr_ = true;
  update_dis_ = true;
}

void Design::InitCorr() {
  kCorrDenominator = 1.0 * n_run_ * (n_run_ + 1) * (2 * n_run_ + 1) / 6 -
    1.0 * n_run_ * (n_run_ + 1) * (n_run_ + 1) / 4;
  corr_.resize(k_var_, std::vector<double>(k_var_));
  for (int i = 0; i < k_var_; ++i) {
    corr_[i][i] = 1;
    for (int j = 0; j < i; ++j) {
      double rho = 0;
      for (int o = 0; o < n_run_; ++o) {
        rho += 1.0 * design_[o][i] * design_[o][j];
      }
      rho -= 1.0 * n_run_ * (n_run_ + 1) * (n_run_ + 1) / 4;
      rho /= kCorrDenominator;
      corr_[i][j] = corr_[j][i] = rho;
    }
  }
}

void Design::InitDis() {
  dis_.resize(n_run_, std::vector<int>(n_run_));
  phi_p_ = 0;
  for (int i = 0; i < n_run_; ++i) {
    dis_[i][i] = 0;
    for (int j = 0; j < i; ++j) {
      int l1_dis = 0;
      for (int o = 0; o < k_var_; ++o) {
        l1_dis += std::abs(design_[i][o] - design_[j][o]);
      }
      dis_[i][j] = dis_[j][i] = l1_dis;
      phi_p_ += 1.0 / QuickPow15(dis_[i][j]);
    }
  }
  phi_p_ = std::pow(phi_p_, 1.0 / kPInPhi);
  InitPhiPBound();
}

void Design::InitPhiPBound() {
  double d_avg = (n_run_ + 1) * k_var_ / 3.0;
  double d_floor = std::floor(d_avg);
  double d_ceil = std::ceil(d_avg);
  phi_p_low_ = std::pow(n_run_ * (n_run_ - 1) / 2 * 
    ((d_ceil - d_avg) / std::pow(d_floor, kPInPhi) + 
    (d_avg - d_floor) / std::pow(d_ceil, kPInPhi)), 1.0 / kPInPhi);
  phi_p_up_ = 0;
  for (int i = 1; i < n_run_; ++i) {
    phi_p_up_ += (n_run_ - i) / QuickPow15(i * k_var_);
  }
  phi_p_up_ = std::pow(phi_p_up_, 1.0 / kPInPhi);
}

void Design::EnableCorr() {
  if (!update_corr_) {
    update_corr_ = true;
    InitCorr();
  }
}

void Design::EnableDis() {
  if (!update_dis_) {
    update_dis_ = true;
    InitDis();
  }
}

void Design::SwapInCol(int col, int r1, int r2) {
  MaintainCorr(col, r1, r2);
  MaintainDis(col, r1, r2);
  std::swap(design_[r1][col], design_[r2][col]);
}

std::pair<double, double> Design::PreSwapInCol(int col, int r1, int r2) const {
  return std::make_pair(GetPreSwapMaxAbsCorr(col, r1, r2), 
    GetPreSwapPhiP(col, r1, r2));
}

double Design::GetMaxAbsCorr(int col) const {
  if (!update_corr_) return -1;
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    max_corr = std::max(max_corr, std::fabs(corr_[i][col]));
  }
  return max_corr;
}

double Design::GetMaxAbsCorr() const {
  if (!update_corr_) return -1;
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    for (int j = 0; j < i; ++j) {
      max_corr = std::max(max_corr, std::fabs(corr_[i][j]));
    }
  }
  return max_corr;
}

double Design::GetMaxAbsCorrExcept(int col) const {
  if (!update_corr_) return -1;
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    for (int j = 0; j < i; ++j) {
      if (j == col) continue;
      max_corr = std::max(max_corr, std::fabs(corr_[i][j]));
    }
  }
  return max_corr;
}

void Design::MaintainCorr(int col, int r1, int r2) {
  if (!update_corr_) return;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    int deta = (design_[r1][i] - design_[r2][i]) *
      (design_[r2][col] - design_[r1][col]);
    corr_[i][col] = corr_[col][i] = 
      (corr_[i][col] * kCorrDenominator + deta) / kCorrDenominator;
  }
}

void Design::MaintainDis(int col, int r1, int r2) {
  if (!update_dis_) return;
  double phi_p_num = QuickPow15(phi_p_);
  auto UpdateDis = [this, &phi_p_num, col](int r1, int r2, int deta) {
    phi_p_num -= 1.0 / QuickPow15(dis_[r1][r2]);
    dis_[r1][r2] -= std::abs(design_[r2][col] - design_[r1][col]);
    dis_[r1][r2] += std::abs(design_[r2][col] - design_[r1][col] + deta);
    phi_p_num += 1.0 / QuickPow15(dis_[r1][r2]);
  };
  int deta = design_[r2][col] - design_[r1][col];
  for (int i = 0; i < n_run_; ++i) {
    if (i == r1 || i == r2) continue;
    i > r1 ? UpdateDis(i, r1, deta) : UpdateDis(r1, i, -deta);
    i > r2 ? UpdateDis(i, r2, -deta) : UpdateDis(r2, i, deta);
  }
  phi_p_ = std::pow(phi_p_num, 1.0 / kPInPhi);
}

double Design::GetPreSwapMaxAbsCorr(int col, int r1, int r2) const {
  if (!update_corr_) return -1;
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    int deta = (design_[r1][i] - design_[r2][i]) *
      (design_[r2][col] - design_[r1][col]);
    double tmp_corr = (corr_[i][col] * kCorrDenominator + deta) / kCorrDenominator;
    max_corr = std::max(max_corr, std::fabs(tmp_corr));
  }
  return max_corr;
}

double Design::GetPreSwapPhiP(int col, int r1, int r2) const {
  if (!update_dis_) return -1;
  double phi_p_num = QuickPow15(phi_p_);
  auto UpdateDis = [this, &phi_p_num, col](int r1, int r2, int deta) {
    phi_p_num -= 1.0 / QuickPow15(dis_[r1][r2]);
    int tmp_dis = dis_[r1][r2] - std::abs(design_[r2][col] - design_[r1][col]) +
      std::abs(design_[r2][col] - design_[r1][col] + deta);
    phi_p_num += 1.0 / QuickPow15(tmp_dis);
  };
  int deta = design_[r2][col] - design_[r1][col];
  for (int i = 0; i < n_run_; ++i) {
    if (i == r1 || i == r2) continue;
    i > r1 ? UpdateDis(i, r1, deta) : UpdateDis(r1, i, -deta);
    i > r2 ? UpdateDis(i, r2, -deta) : UpdateDis(r2, i, deta);
  }
  return std::pow(phi_p_num, 1.0 / kPInPhi);
}

void Design::Display() {
  std::cout << n_run_ << " " << k_var_ << "\n";
  int w = std::floor(std::log10(n_run_)) + 1;
  for (int i = 0; i < n_run_; ++i) {
    for (int j = 0; j < k_var_; ++j) {
      std::cout << std::setw(w) << design_[i][j] << " \n"[j == k_var_ - 1];
    }
  }
  if (!update_dis_) InitDis();
  if (!update_corr_) InitCorr();
  std::cout << "PhiP: " << GetPhiP() << "\n";
  std::cout << "RhoMax: " << GetMaxAbsCorr() << std::endl;
}
} // namespace SearchingAlgorithm