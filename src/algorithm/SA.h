#pragma once

#include "search_algorithm.h"

/*
Exploratory designs for computational experiments
https://sci-hub.yncjkj.com/10.1016/0378-3758(94)00035-t
*/

namespace LHD {
class SA : public SearchAlgorithm {
 public:
  SA(int n, int k);
  inline void SetInitTemp(double init_temp) { init_temp_ = init_temp; }
  inline void SetRate(double rate) { rate_ = rate; }
  inline void SetIMax(int i_max) { i_max_ = i_max; }
  Design::VecInt2D Search() override;
  Design::VecInt2D IncrementalSearch(Design design);
 protected:
  void InitDefaultParam();
 private:
  double init_temp_;
  double rate_;
  int i_max_;
};
} // namespace LHD