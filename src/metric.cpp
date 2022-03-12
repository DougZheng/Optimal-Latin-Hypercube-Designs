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

class Design {
 public:
  using VecInt2D = std::vector<std::vector<int>>;
  using VecDouble2D = std::vector<std::vector<double>>;
  Design(const VecInt2D& design);
  Design(int n, int k, int seed);
  inline int GetN() { return n_run_; }
  inline int GetK() { return k_var_; }
  inline VecInt2D GetDesign() { return design_; }
  double CalcPhiP(int p = 15, int t = 1);
  int CalcMinL1();
  double CalcMinL2();
  double CalcRhoMax();
  double CalcRhoAvg();
  void Display();
 private:
  int n_run_;
  int k_var_;
  VecInt2D design_;
  double phi_p_;
  double min_l1_;
  double min_l2_;
  double rho_max_;
  double rho_avg_;
};

Design::Design(const VecInt2D& design) : design_(design) {
  assert(design.size() > 0);
  n_run_ = design.size();
  k_var_ = design[0].size();
}

Design::Design(int n, int k, int seed) : n_run_(n), k_var_(k) {
  std::mt19937 rng_(seed);
  design_.resize(n);
  for (int i = 0; i < n; ++i) {
    design_.resize(k);
  }
  std::vector<int> rnd_list(n);
  std::iota(rnd_list.begin(), rnd_list.end(), 1);
  for (int i = 0; i < k; ++i) {
    std::shuffle(rnd_list.begin(), rnd_list.end(), rng_);
    for (int j = 0; j < n; ++j) {
      design_[j][i] = rnd_list[j];
    }
  }
}

double Design::CalcPhiP(int p, int t) {
  if (n_run_ <= 1) {
    return 0;
  }
  assert(t == 1 || t == 2);
  int min_num = design_[0][0];
  int max_num = design_[0][0];
  for (int i = 1; i < n_run_; ++i) {
    min_num = std::min(min_num, design_[i][0]);
    max_num = std::max(max_num, design_[i][0]);
  }
  auto GetNum = [min_num, max_num](int x) -> double {
    return x;
    return static_cast<double>(x - min_num) / (max_num - min_num);
  };
  double phi_p = 0;
  for (int i = 0; i < n_run_; ++i) {
    for (int j = i + 1; j < n_run_; ++j) {
      double dij = 0;
      for (int o = 0; o < k_var_; ++o) {
        double aio = GetNum(design_[i][o]);
        double ajo = GetNum(design_[j][o]);
        if (t == 1) dij += std::fabs(aio - ajo);
        else dij += (aio - ajo) * (aio - ajo);
      }
      double d = (t == 1 ? dij : std::sqrt(dij));
      phi_p += std::pow(d, -p);
    }
  }
  phi_p = std::pow(phi_p, 1.0 / p);
  return phi_p;
}

int Design::CalcMinL1() {
  if (n_run_ <= 1) {
    return 0;
  }
  int min_l1 = __INT_MAX__;
  for (int i = 0; i < n_run_; ++i) {
    for (int j = i + 1; j < n_run_; ++j) {
      int dij = 0;
      for (int o = 0; o < k_var_; ++o) {
        dij += std::abs(design_[i][o] - design_[j][o]);
      }
      min_l1 = std::min(min_l1, dij);
    }
  }
  return min_l1;
}

double Design::CalcMinL2() {
  if (n_run_ <= 1) {
    return 0;
  }
  double min_l2 = __DBL_MAX__;
  for (int i = 0; i < n_run_; ++i) {
    for (int j = i + 1; j < n_run_; ++j) {
      double dij = 0;
      for (int o = 0; o < k_var_; ++o) {
        dij += (design_[i][o] - design_[j][o]) * (design_[i][o] - design_[j][o]);
      }
      dij = std::sqrt(dij);
      min_l2 = std::min(min_l2, dij);
    }
  }
  return min_l2;
}

double Design::CalcRhoMax() {
  if (k_var_ <= 1) {
    return 0;
  }
  double rho_max = 0;
  std::vector<double> avg(k_var_);
  for (int i = 0; i < k_var_; ++i) {
    for (int j = 0; j < n_run_; ++j) {
      avg[i] += design_[j][i];
    }
    avg[i] /= static_cast<double>(n_run_);
  }
  std::vector<double> var(k_var_);
  for (int i = 0; i < n_run_; ++i) {
    for (int j = 0; j < k_var_; ++j) {
      var[j] += (design_[i][j] - avg[j]) * (design_[i][j] - avg[j]);
    }
  }
  for (int i = 0; i < k_var_; ++i) {
    for (int j = i + 1; j < k_var_; ++j) {
      double rho = 0;
      for (int o = 0; o < n_run_; ++o) {
        rho += (design_[o][i] - avg[i]) * (design_[o][j] - avg[j]);
      }
      rho /= std::sqrt(var[i] * var[j]);
      rho_max = std::max(rho_max, std::fabs(rho));
    }
  }
  return rho_max;
}

double Design::CalcRhoAvg() {
  if (k_var_ <= 1) {
    return 0;
  }
  double rho_avg = 0;
  std::vector<double> avg(k_var_);
  for (int i = 0; i < k_var_; ++i) {
    for (int j = 0; j < n_run_; ++j) {
      avg[i] += design_[j][i];
    }
    avg[i] /= static_cast<double>(n_run_);
  }
  std::vector<double> var(k_var_);
  for (int i = 0; i < n_run_; ++i) {
    for (int j = 0; j < k_var_; ++j) {
      var[j] += (design_[i][j] - avg[j]) * (design_[i][j] - avg[j]);
    }
  }
  for (int i = 0; i < k_var_; ++i) {
    for (int j = i + 1; j < k_var_; ++j) {
      double rho = 0;
      for (int o = 0; o < n_run_; ++o) {
        rho += (design_[o][i] - avg[i]) * (design_[o][j] - avg[j]);
      }
      rho /= std::sqrt(var[i] * var[j]);
      rho_avg += 2 * std::fabs(rho);
    }
  }
  rho_avg /= k_var_ * (k_var_ - 1);
  return rho_avg;
}

void Design::Display() {
  std::cout << n_run_ << " " << k_var_ << "\n";
  int w = std::floor(std::log10(n_run_)) + 1;
  for (int i = 0; i < n_run_; ++i) {
    for (int j = 0; j < k_var_; ++j) {
      std::cout << std::setw(w) << design_[i][j] << " \n"[j == k_var_ - 1];
    }
  }
  std::cout << "PhiP: " << CalcPhiP() << "\n";
  std::cout << "MaximinL1: " << CalcMinL1() << "\n";
  std::cout << "MaximinL2: " << CalcMinL2() << "\n";
  std::cout << "RhoMax: " << CalcRhoMax() << "\n";
  std::cout << "RhoAvg: " << CalcRhoAvg() << std::endl;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cout << "usage: ${exe} ${file}" << std::endl;
    return 1;
  }
  std::ifstream in(argv[1]);
  if (in.fail()) {
    std::cout << "open " << argv[1] << " failed." << std::endl;
    return 1;
  }
  int n, k;
  in >> n >> k;
  Design::VecInt2D a(n, std::vector<int>(k));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      in >> a[i][j];
    }
  }
  Design(a).Display();
  return 0;
}