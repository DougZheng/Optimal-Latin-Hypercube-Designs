#include <cmath>
#include <ctime>
#include <cassert>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
#include <random>

namespace ESEJin2005 {
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
  inline double GetPhiP() const { return phi_p_; };
  inline double GetCritVal(double w, double rho_max, double phi_p) {
    return w * rho_max + 
      (1 - w) * (phi_p - phi_p_low_) / (phi_p_up_ - phi_p_low_);
  }
  void SwapInCol(int col, int r1, int r2);
  std::pair<double, double> PreSwapInCol(int col, int r1, int r2) const;
  double GetMaxAbsCorr(int col) const;
  double GetMaxAbsCorr() const;
  double GetMaxAbsCorrExcept(int col) const;
  void Display() const;
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
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    if (i == col) continue;
    max_corr = std::max(max_corr, std::fabs(corr_[i][col]));
  }
  return max_corr;
}

double Design::GetMaxAbsCorr() const {
  double max_corr = 0;
  for (int i = 0; i < k_var_; ++i) {
    for (int j = 0; j < i; ++j) {
      max_corr = std::max(max_corr, std::fabs(corr_[i][j]));
    }
  }
  return max_corr;
}

double Design::GetMaxAbsCorrExcept(int col) const {
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
    corr_[i][col] = corr_[col][i] = 
      (corr_[i][col] * kCorrDenominator + deta) / kCorrDenominator;
  }
}

void Design::MaintainDis(int col, int r1, int r2) {
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

void Design::Display() const {
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

class Searcher {
 public:
  Searcher(int n, int k);
  inline void SetMCol(int m_col) { m_col_ = m_col; }
  inline void SetJColPair(int j_col_pair) { j_col_pair_ = j_col_pair; }
  inline void SetAlpha1(double alpha1) { alpha1_ = alpha1; }
  inline void SetAlpha2(double alpha2) { alpha2_ = alpha2; }
  inline void SetAlpha3(double alpha3) { alpha3_ = alpha3; }
  inline void SetW(double w) { w_ = w; }
  inline void SetIterateCnt(int iterate_cnt) { iterate_cnt_ = iterate_cnt; }
  inline void SetSeed(int seed) { rng_.seed(seed); }
  inline void SetPrintFrequence(int print_frequence) { print_frequence_ = print_frequence; }
  Design Search();
  Design IncrementalSearch(Design design);
 private:
  void InitDefaultParam();
  void ShuffleM(std::vector<std::pair<int, int>>& pair_list, int m);
 private:
  int n_;
  int k_;
  int m_col_;
  int j_col_pair_;
  double alpha1_;
  double alpha2_;
  double alpha3_;
  double w_;
  int iterate_cnt_;
  std::mt19937 rng_;
  int print_frequence_;
};

Searcher::Searcher(int n, int k) : n_(n), k_(k) {
  InitDefaultParam();
}

void Searcher::InitDefaultParam() {
  int n_e = n_ * (n_ - 1) / 2;
  j_col_pair_ = std::min(50, (n_e + 4) / 5);
  m_col_ = k_;
  alpha1_ = 0.8;
  alpha2_ = 0.9;
  alpha3_ = 0.7;
  rng_.seed(0);
  w_ = 0.5;
  iterate_cnt_ = 1000;
  print_frequence_ = 100;
}

Design Searcher::Search() {
  Design design(n_, k_, rng_());
  return IncrementalSearch(design);
}

Design Searcher::IncrementalSearch(Design design) {
  assert(design.GetN() == n_ && design.GetK() == k_);
  std::uniform_real_distribution<double> uniform_dis(0.0, std::nextafter(1.0, 1.1));
  int cnt = 0;
  double opt_corr = design.GetMaxAbsCorr();
  double opt_val = design.GetCritVal(w_, opt_corr, design.GetPhiP());
  Design::VecInt2D opt_design = design.GetDesign();
  double T_h = 0.005 * opt_val;
  std::vector<std::pair<int, int>> pair_list;
  pair_list.reserve(n_ * (n_ - 1) / 2);
  for (int i = 0; i < n_; ++i) {
    for (int j = i + 1; j < n_; ++j) {
      pair_list.emplace_back(i, j);
    }
  }
  double cur_val = opt_val;
  double cur_corr = opt_corr;
  int max_no_imp_cnt = 0;
  int explore_dir = 0;
  auto StopSearch = [this, &cnt, &max_no_imp_cnt]() -> bool {
    return iterate_cnt_ >= 0 && cnt >= iterate_cnt_ || max_no_imp_cnt >= 1000;
  };
  auto PrintLog = [this, &cnt, &cur_val, &cur_corr, &design]() -> void {
    std::cerr << "Iteration: " << cnt
      << ", Val: " << cur_val 
      << ", PhiP: " << design.GetPhiP()
      << ", RhoMax: " << cur_corr << std::endl;
  };
  while (true) {
    bool stop = StopSearch();
    if (stop || cnt % print_frequence_ == 0) {
      PrintLog();
      if (stop) break;
    }
    ++cnt;

    double old_val = opt_val;
    int n_acpt = 0;
    int n_imp = 0;
    for (int i = 0; i < m_col_; ++i) {
      int col = i % k_;
      double inner_opt_val = __DBL_MAX__;
      double inner_opt_corr = __DBL_MAX__;
      ShuffleM(pair_list, j_col_pair_);
      std::pair<int, int> opt_pair;
      double corr_except = design.GetMaxAbsCorrExcept(col);
      for (int j = 0; j < j_col_pair_; ++j) {
        auto pair = pair_list[j];
        auto tmp_ret = design.PreSwapInCol(col, pair.first, pair.second);
        double tmp_corr = std::max(corr_except, tmp_ret.first);
        double tmp_val = design.GetCritVal(w_, tmp_corr, tmp_ret.second);
        if (tmp_val < inner_opt_val) {
          inner_opt_val = tmp_val;
          inner_opt_corr = tmp_corr;
          opt_pair = pair;
        }
      }
      if (inner_opt_val - cur_val <= T_h * uniform_dis(rng_)) {
        design.SwapInCol(col, opt_pair.first, opt_pair.second);
        ++n_acpt;
        if (inner_opt_val < opt_val) {
          opt_val = inner_opt_val;
          opt_design = design.GetDesign();
          ++n_imp;
        }
        cur_val = inner_opt_val;
        cur_corr = inner_opt_corr;
      }
    }

    bool flag_imp = old_val - opt_val > 0;
    max_no_imp_cnt = flag_imp ? 0 : max_no_imp_cnt + 1;
    double acpt_ratio = 1.0 * n_acpt / m_col_;
    double imp_ratio = 1.0 * n_imp / m_col_;
    if (flag_imp) {
      if (acpt_ratio > 0.1 && n_acpt > n_imp) {
        T_h *= alpha1_;
      } else if (acpt_ratio > 0.1 && n_acpt == n_imp) {
        ;
      } else {
        T_h /= alpha1_;
      }
    } else {
      if (acpt_ratio < 0.1) {
        T_h /= alpha3_;
        explore_dir = 1;
      } else if (acpt_ratio > 0.8) {
        T_h *= alpha2_;
        explore_dir = -1;
      } else if (explore_dir == 1) {
        T_h /= alpha3_;
      } else if (explore_dir == -1) {
        T_h *= alpha2_;
      }
    }
  }
  return Design(opt_design);
}

void Searcher::ShuffleM(std::vector<std::pair<int, int>>& pair_list, int m) {
  int n = pair_list.size();
  assert(m <= n);
  for (int i = 0; i < m; ++i) {
    int j = rng_() % (n - i);
    std::swap(pair_list[i], pair_list[j]);
  }
}
} // namespace ESEJin2005

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cout << "Usage: ${exe} ${n} ${k} ${w} ${cnt} [${file}]" << std::endl;
    return 1;
  }
  // ./ESE_Jin2005 20 8 0.5 1000 >LHD
  // ./ESE_Jin2005 20 8 0.5 1000 LHD
  int n = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  double w = std::stod(argv[3]);
  int cnt = std::stoi(argv[4]);
  auto searcher = ESEJin2005::Searcher(n, k);
  searcher.SetW(w);
  searcher.SetIterateCnt(cnt);
  if (argc > 5) {
    auto design = ESEJin2005::Design::ReadDesign(argv[5]);
    searcher.IncrementalSearch(design).Display();
  } else {
    searcher.Search().Display();
  }
  return 0;
}