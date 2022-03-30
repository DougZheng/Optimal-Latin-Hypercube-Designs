#pragma once

#include "search_algorithm.h"

/*
An efficient algorithm for constructing optimal design of computer experiments
https://sci-hub.yncjkj.com/https://www.sciencedirect.com/science/article/pii/S0378375804001922
*/

namespace LHD {
class ESE : public SearchAlgorithm {
 public:
  ESE(int n, int k);
  inline void SetMCol(int m_col) { m_col_ = m_col; }
  inline void SetJColPair(int j_col_pair) { j_col_pair_ = j_col_pair; }
  inline void SetAlpha1(double alpha1) { alpha1_ = alpha1; }
  inline void SetAlpha2(double alpha2) { alpha2_ = alpha2; }
  inline void SetAlpha3(double alpha3) { alpha3_ = alpha3; }
  Design::VecInt2D Search() override;
  Design::VecInt2D IncrementalSearch(const Design::VecInt2D& ori_design);
 protected:
  void InitDefaultParam();
 private:
  int m_col_;
  int j_col_pair_;
  double alpha1_;
  double alpha2_;
  double alpha3_;
};
} // namespace LHD