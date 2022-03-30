#pragma once

#include "construction_algorithm.h"

/*
Optimal maximin L1-distance Latin hypercube designs based on good lattice point designs
http://www.stat.ucla.edu/~hqxu/pub/WangXiaoXu2018.pdf
*/

namespace LHD {
class Wang2018 : public ConstructionAlgorithm {
 private:
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
 public:
  Design::VecInt2D Solve(int n, int k) override; // For any n x k design
  void DisplayResult(const Design::VecInt2D& d) const override {
    Design(d).Display();
  }
 private:
  inline bool IsPrime(int n) { 
    return GetPhiN(n) == n - 1;
  }
  int GetPhiN(int n);
  std::vector<int> GetCoPrimeList(int n);
  Design GetInitDesign(int n, int k);
  Design Algorithm1(int n, int k);        // n x k design, n > 0 and k <= φ(n)
  Design Algorithm2(int n);               // n x n design, 2n + 1 is a prime
  Design Algorithm3(int n);               // (n + 1) x n design, 2n + 1 is a prime
  Design Algorithm4(int m, int n, int k); // n x k design from m x k design
};

} // namespace LHD