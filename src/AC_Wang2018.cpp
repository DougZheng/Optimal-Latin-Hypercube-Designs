#include <cmath>
#include <cassert>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>

/*
Optimal maximin L1-distance Latin hypercube designs based on good lattice point designs
http://www.stat.ucla.edu/~hqxu/pub/WangXiaoXu2018.pdf

Time cost: O(kn^3logn)
Memory cost: O(kn)
*/

namespace ACWang2018 {
class Design {
 public:
  using VecInt2D = std::vector<std::vector<int>>;
  Design(const VecInt2D& a);
  Design(VecInt2D&& a);
  inline Design Copy() const { return *this; };
  inline VecInt2D GetA() const { return a_; }
  inline int GetN() const { return n_; }
  inline int GetK() const { return k_; }
  int GetMinL1();
  Design& WilliamTrans();
  Design& ModifiedWilliamTrans();
  Design& CyclicPlusB(int b);
  Design& Resize(int n, int k);
  Design& AppendZeros();
  void Display();
 private:
  VecInt2D a_;
  int n_;
  int k_;
};

Design::Design(const VecInt2D& a) {
  assert(a.size() > 0);
  n_ = a.size();
  k_ = a[0].size();
  a_ = a;
}

Design::Design(VecInt2D&& a) {
  assert(a.size() > 0);
  n_ = a.size();
  k_ = a[0].size();
  a_ = std::move(a);
}

int Design::GetMinL1() {
  if (n_ == 1) {
    return 0;
  }
  int min_l1 = __INT_MAX__;
  for (int i = 0; i < n_; ++i) {
    for (int j = i + 1; j < n_; ++j) {
      int l1 = 0;
      for (int o = 0; o < k_ && l1 < min_l1; ++o) {
        l1 += std::abs(a_[i][o] - a_[j][o]);
      }
      min_l1 = std::min(min_l1, l1);
    }
  }
  return min_l1;
}

Design& Design::WilliamTrans() {
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < k_; ++j) {
      a_[i][j] = a_[i][j] < (n_ + 1) / 2 ? 
        2 * a_[i][j] : 
        2 * (n_ - a_[i][j]) - 1;
    }
  }
  return *this;
}

Design& Design::ModifiedWilliamTrans() {
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < k_; ++j) {
      a_[i][j] = a_[i][j] < (n_ + 1) / 2 ? 
        2 * a_[i][j] : 2 * (n_ - a_[i][j]);
    }
  }
  return *this;
}

Design& Design::CyclicPlusB(int b) {
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < k_; ++j) {
      (a_[i][j] += b) %= n_;
    }
  }
  return *this;
}

Design& Design::Resize(int n, int k) {
  if (k_ != k) {
    for (int i = 0; i < n_; ++i) {
      while (a_[i].size() > k) {
        a_[i].pop_back();
      }
    }
    k_ = k;
  }
  if (n_ != n) {
    while (a_.size() > n) {
      a_.pop_back();
    }
    std::vector<int> cnt(n_, 0);
    for (int i = 0; i < k; ++i) {
      for (int j = 0; j < n; ++j) {
        ++cnt[a_[j][i]];
      }
      for (int j = 1; j < n_; ++j) {
        cnt[j] += cnt[j - 1];
      }
      for (int j = 0; j < n; ++j) {
        a_[j][i] = cnt[a_[j][i]] - 1;
      }
      std::fill(cnt.begin(), cnt.end(), 0);
    }
    n_ = n;
  }
  return *this;
}

Design& Design::AppendZeros() {
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < k_; ++j) {
      ++a_[i][j];
    }
  }
  a_.emplace_back(std::vector<int>(k_, 0));
  ++n_;
  return *this;
}

void Design::Display() {
  int w = std::floor(std::log10(n_)) + 1;
  for (int i = 0; i < n_; ++i) {
    for (int j = 0; j < k_; ++j) {
      std::cout << std::setw(w) << (a_[i][j] + 1) << " \n"[j == k_ - 1];
    }
  }
  std::cout << "maximin l1: " << GetMinL1() << std::endl;
}

class Solver {
 public:
  Design Algorithm1(int n, int k); // n x k design, n > 0 and k <= φ(n)
  Design Algorithm2(int n); // n x n design, 2n + 1 is a prime
  Design Algorithm3(int n); // (n + 1) x n design, 2n + 1 is a prime
  Design Algorithm4(int m, int n, int k); // n x k design from m x k design
  Design Solve(int n, int k); // For any n x k design
 private:
  inline bool IsPrime(int n) { return GetPhiN(n) == n - 1; }
  int GetPhiN(int n);
  std::vector<int> GetCoPrimeList(int n);
  Design GetInitDesign(int n, int k);
};

Design Solver::Algorithm1(int n, int k) {
  Design design = GetInitDesign(n, k);
  int maximin_l1 = -1;
  int opt_b = 0;
  bool opt_wt = false;
  for (int i = 0; i < n; ++i) {
    Design tmp_design = design.Copy().CyclicPlusB(i);
    int d1_l1 = tmp_design.GetMinL1();
    int d2_l1 = tmp_design.WilliamTrans().GetMinL1();
    if (d1_l1 > maximin_l1) {
      maximin_l1 = d1_l1;
      opt_b = i;
      opt_wt = false;
    }
    if (d2_l1 > maximin_l1) {
      maximin_l1 = d2_l1;
      opt_b = i;
      opt_wt = true;
    }
  }
  design.CyclicPlusB(opt_b);
  if (opt_wt) {
    design.WilliamTrans();
  }
  return design;
}

Design Solver::Algorithm2(int n) {
  assert(IsPrime(2 * n + 1));
  Design design = GetInitDesign(2 * n + 1, 2 * n);
  return design.ModifiedWilliamTrans().Resize(n, n);
}

Design Solver::Algorithm3(int n) {
  assert(IsPrime(2 * n + 1));
  Design design = GetInitDesign(2 * n + 1, 2 * n);
  return design.ModifiedWilliamTrans().Resize(n, n).AppendZeros();
}

Design Solver::Algorithm4(int m, int n, int k) {
  Design design = GetInitDesign(m, k);
  int maximin_l1 = -1;
  int opt_b = 0;
  bool opt_wt = false;
  for (int i = 0; i < m; ++i) {
    Design tmp_design = design.Copy().CyclicPlusB(i);
    int d1_l1 = tmp_design.Copy().Resize(n, k).GetMinL1();
    int d2_l1 = tmp_design.WilliamTrans().Resize(n, k).GetMinL1();
    if (d1_l1 > maximin_l1) {
      maximin_l1 = d1_l1;
      opt_b = i;
      opt_wt = false;
    }
    if (d2_l1 > maximin_l1) {
      maximin_l1 = d2_l1;
      opt_b = i;
      opt_wt = true;
    }
  }
  design.CyclicPlusB(opt_b);
  if (opt_wt) {
    design.WilliamTrans();
  }
  return design.Resize(m, k);
}

Design Solver::Solve(int n, int k) {
  std::vector<Design> designs;

  if (GetPhiN(n) >= k) {
    designs.emplace_back(Algorithm1(n, k));
  }

  int m = n;
  do {
    ++m;
    if (GetPhiN(m) >= k) {
      designs.emplace_back(Algorithm4(m, n, k));
    }
  } while (!IsPrime(m) || k > m - 1);

  m = n - 1;
  while (!IsPrime(2 * m + 1) || k > m) {
    ++m;
  }
  designs.emplace_back(Algorithm3(m).Resize(n, k));

  int maximin_l1 = -1;
  int design_idx = 0;
  for (int i = 0; i < designs.size(); ++i) {
    int l1 = designs[i].GetMinL1();
    if (l1 > maximin_l1) {
      maximin_l1 = l1;
      design_idx = i;
    }
  }
  return designs[design_idx];
}

int Solver::GetPhiN(int n) {
  int phi_n = n;
  for (int i = 2; i * i <= n; ++i) {
    if (n % i) continue;
    phi_n = phi_n / i * (i - 1);
    while (n % i == 0) {
      n /= i;
    }
  }
  if (n > 1) {
    phi_n = phi_n / n * (n - 1);
  }
  return phi_n;
}

std::vector<int> Solver::GetCoPrimeList(int n) {
  std::vector<int> ret_list;
  for (int i = 1; i <= n; ++i) {
    if (std::__gcd(n, i) == 1) {
      ret_list.push_back(i);
    }
  }
  return ret_list;
}

Design Solver::GetInitDesign(int n, int k) {
  std::vector<int> coprime_list = GetCoPrimeList(n);
  assert(k <= coprime_list.size());
  Design::VecInt2D a(n, std::vector<int>(k));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      a[i][j] = (i + 1) * coprime_list[j] % n;
    }
  }
  return Design(std::move(a));
}
} // namespace ACWang2018

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cout << "Usage: ${exe} ${n} ${k}" << std::endl;
    return 1;
  }
  int n = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  // ACWang2018::Solver().Algorithm1(n, k).Display();
  // ACWang2018::Solver().Algorithm2(n).Display();
  // ACWang2018::Solver().Algorithm3(n).Display();
  // ACWang2018::Solver().Algorithm4(n, k).Display();
  ACWang2018::Solver().Solve(n, k).Display();
  return 0;
}