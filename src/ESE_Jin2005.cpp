#include <cmath>
#include <ctime>
#include <cassert>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>
#include <random>

namespace ESEJin2005 {
class Design {
 public:
  using VecInt2D = std::vector<std::vector<int>>;
  using VecDouble2D = std::vector<std::vector<double>>;
  Design(const VecInt2D& design);
  Design(int n, int k, int seed);
  inline int GetN() { return n_run_; }
  inline int GetK() { return k_var_; }
  inline VecInt2D GetDesign() { return design_; }
  inline double GetPhiP() { return phi_p_; };
  inline double GetCritVal(double w, double rho_max) {
    static auto phi_p_bound = [this]() -> std::pair<double, double> {
      double d_avg = (n_run_ + 1) * k_var_ / 3.0;
      double d_floor = std::floor(d_avg);
      double d_ceil = std::ceil(d_avg);
      double phi_p_low = std::pow(n_run_ * (n_run_ - 1) / 2 * 
        ((d_ceil - d_avg) / std::pow(d_floor, kPInPhi) + 
        (d_avg - d_floor) / std::pow(d_ceil, kPInPhi)), 1.0 / kPInPhi);
      double phi_p_up = 0;
      for (int i = 1; i < n_run_; ++i) {
        phi_p_up += (n_run_ - i) * std::pow(i * k_var_, kPInPhi);
      }
      phi_p_up = std::pow(phi_p_up, 1.0 / kPInPhi);
      return {phi_p_low, phi_p_up};
    }();
    return w * rho_max + 
      (1 - w) * (phi_p_ - phi_p_bound.first) / (phi_p_bound.second - phi_p_bound.first);
  }
  void SwapInCol(int col, int r1, int r2);
  double GetMaxAbsCorr(int col);
  double GetMaxAbsCorr();
  double GetMaxAbsCorrExcept(int col);
  void Display();
 private:
  void InitCorr();
  void InitDis();
  void MaintainCorr(int col, int r1, int r2);
  void MaintainDis(int col, int r1, int r2);
  double QuickPow(double x, int p);
 private:
  const int kPInPhi = 15;
  double kCorrDenominator;
  int n_run_;
  int k_var_;
  VecInt2D design_;
  VecDouble2D corr_;
  VecInt2D dis_;
  double phi_p_;
};

Design::Design(const VecInt2D& design) 
    : design_(design) {
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
      phi_p_ += QuickPow(dis_[i][j], -kPInPhi);
    }
  }
  phi_p_ = std::pow(phi_p_, 1.0 / kPInPhi);
}

void Design::SwapInCol(int col, int r1, int r2) {
  MaintainCorr(col, r1, r2);
  MaintainDis(col, r1, r2);
  std::swap(design_[r1][col], design_[r2][col]);
}

double Design::GetMaxAbsCorr(int col) {
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    max_corr = std::max(max_corr, std::fabs(corr_[i][col]));
  }
  return max_corr;
}

double Design::GetMaxAbsCorr() {
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    for (int j = 0; j < i; ++j) {
      max_corr = std::max(max_corr, std::fabs(corr_[i][j]));
    }
  }
  return max_corr;
}

double Design::GetMaxAbsCorrExcept(int col) {
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
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    int deta = (design_[r1][i] - design_[r2][i]) *
      (design_[r2][col] - design_[r1][col]);
    corr_[i][col] = corr_[col][i] = (corr_[i][col] * kCorrDenominator + deta) / kCorrDenominator;
  }
}

void Design::MaintainDis(int col, int r1, int r2) {
  double phi_p_num = QuickPow(phi_p_, kPInPhi);
  auto UpdateDis = [this, &phi_p_num, col](int r1, int r2, int deta) {
    if (r1 <= r2) return;
    phi_p_num -= QuickPow(dis_[r1][r2], -kPInPhi);
    dis_[r1][r2] -= std::abs(design_[r2][col] - design_[r1][col]);
    dis_[r1][r2] += std::abs(design_[r2][col] - design_[r1][col] + deta);
    phi_p_num += QuickPow(dis_[r1][r2], -kPInPhi);
  };
  int deta = design_[r2][col] - design_[r1][col];
  for (int i = 0; i < n_run_; ++i) {
    if (i == r1 || i == r2) continue;
    UpdateDis(i, r1, deta);
    UpdateDis(i, r2, -deta);
    UpdateDis(r1, i, -deta);
    UpdateDis(r2, i, deta);
  }
  phi_p_ = std::pow(phi_p_num, 1.0 / kPInPhi);
}

void Design::Display() {
  std::cout << n_run_ << " " << k_var_ << "\n";
  int w = std::floor(std::log10(n_run_)) + 1;
  for (int i = 0; i < n_run_; ++i) {
    for (int j = 0; j < k_var_; ++j) {
      std::cout << std::setw(w) << design_[i][j] << " \n"[j == k_var_ - 1];
    }
  }
  std::cout << "PhiP: " << GetPhiP() << "\n";
  std::cout << "RhoMax: " << GetMaxAbsCorr() << std::endl;
}

double Design::QuickPow(double x, int p) {
  bool is_neg = p < 0;
  p = std::abs(p);
  double ret = 1;
  while (p > 0) {
    if (p & 1) {
      ret = ret * x;
    }
    x = x * x;
    p >>= 1;
  }
  return is_neg ? 1.0 / ret : ret;
}

class Solver {
 public:
  Solver(int seed = 19937);
  Design Solve(int n, int k, double w = 0.5, int limit_cnt = 10000);
 private:
  void ShuffleM(std::vector<std::pair<int, int>>& pair_list, int m);
 private:
  std::mt19937 rng_;
};

Solver::Solver(int seed) {
  rng_.seed(seed);
}

Design Solver::Solve(int n, int k, double w, int limit_cnt) {
  Design design(n, k, rng_());
  const int kPrintRatio = 100;
  int n_e = n * (n - 1) / 2;
  int J = std::min(50, n_e / 5);
  int M = std::min(100, 2 * n_e * k / J);
  double a1 = 0.8;
  double a2 = 0.9;
  double a3 = 0.7;
  int explore_dir = 0;
  std::uniform_real_distribution<double> uniform_dis(0.0, std::nextafter(1.0, 1.1));
  int iterator_cnt = 0;
  double bst_corr = design.GetMaxAbsCorr();
  double bst_val = design.GetCritVal(w, bst_corr);
  Design::VecInt2D bst_design = design.GetDesign();
  double T_h = 0.005 * bst_val;
  std::vector<std::pair<int, int>> pair_list;
  pair_list.reserve(n * (n - 1) / 2);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      pair_list.emplace_back(i, j);
    }
  }
  double cur_val = bst_val;
  double cur_corr = bst_corr;
  while (iterator_cnt < limit_cnt) {
    ++iterator_cnt;
    if (iterator_cnt % kPrintRatio == 0) {
      std::cerr << "iterator: " << iterator_cnt << 
        ", current val: " << cur_val <<
        ", rho max: " << cur_corr << "\n";
      std::cerr << "T_h: " << T_h << std::endl;
    }
    double old_val = bst_val;
    int n_acpt = 0;
    int n_imp = 0;
    for (int i = 0; i < M; ++i) {
      int col = i % k;
      double inner_bst_val = __DBL_MAX__;
      double inner_bst_corr = __DBL_MAX__;
      ShuffleM(pair_list, J);
      int bst_r1 = 0;
      int bst_r2 = 0;
      double corr_except = design.GetMaxAbsCorrExcept(col);
      for (int j = 0; j < J; ++j) {
        int r1 = pair_list[j].first;
        int r2 = pair_list[j].second;
        design.SwapInCol(col, r1, r2);
        double tmp_corr = std::max(corr_except, design.GetMaxAbsCorr(col));
        double tmp_val = design.GetCritVal(w, tmp_corr);
        if (tmp_val < inner_bst_val) {
          inner_bst_val = tmp_val;
          inner_bst_corr = tmp_corr;
          bst_r1 = r1;
          bst_r2 = r2;
        }
        design.SwapInCol(col, r1, r2);
      }
      design.SwapInCol(col, bst_r1, bst_r2);
      if (inner_bst_val - cur_val <= T_h * uniform_dis(rng_)) {
        ++n_acpt;
        if (inner_bst_val < bst_val) {
          bst_val = inner_bst_val;
          bst_design = design.GetDesign();
          ++n_imp;
        }
        cur_val = inner_bst_val;
        cur_corr = inner_bst_corr;
      }
      else {
        design.SwapInCol(col, bst_r1, bst_r2);
      }
    }
    bool flag_imp = old_val - bst_val > 0;
    double acpt_ratio = static_cast<double>(n_acpt) / M;
    double imp_ratio = static_cast<double>(n_imp) / M;
    if (flag_imp) {
      if (acpt_ratio > 0.1 && n_acpt > n_imp) {
        T_h *= a1;
      }
      else if (acpt_ratio > 0.1 && n_acpt == n_imp) {
        ;
      }
      else {
        T_h /= a1;
      }
    }
    else {
      if (acpt_ratio < 0.1) {
        T_h /= a3;
        explore_dir = 1;
      }
      else if (acpt_ratio > 0.8) {
        T_h *= a2;
        explore_dir = -1;
      }
      else if (explore_dir == 1) {
        T_h /= a3;
      }
      else if (explore_dir == -1) {
        T_h *= a2;
      }
    }
  }
  return Design(bst_design);
}

void Solver::ShuffleM(std::vector<std::pair<int, int>>& pair_list, int m) {
  int n = pair_list.size();
  assert(m <= n);
  for (int i = 0; i < m; ++i) {
    int j = rng_() % (n - i);
    std::swap(pair_list[i], pair_list[j]);
  }
}
} // namespace ESEJin2005

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cout << "Usage: ${exe} ${n} ${k} [$w ${cnt}]" << std::endl;
    return 1;
  }
  int n = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  double w = argc >= 4 ? std::stod(argv[3]) : 0.5;
  int cnt = argc >= 5 ? std::stoi(argv[4]) : 10000;
  ESEJin2005::Solver().Solve(n, k, w, cnt).Display();
  return 0;
}