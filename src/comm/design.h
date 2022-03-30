#pragma once

#include <string>
#include <vector>

#include "criterion/multicriteria.h"
#include "utils.h"

namespace LHD {
class Design {
 public:
  using VecInt2D = std::vector<std::vector<int>>;
  static VecInt2D ReadDesign(const std::string& file);
  Design(const VecInt2D& design);
  Design(int n, int k, int seed = 0);
  Design(const Design& design) {
    *this = design;
  }
  Design& operator=(const Design& design);
  inline int GetN() const { 
    return n_run_; 
  }
  inline int GetK() const { 
    return k_var_;
  }
  inline VecInt2D GetDesign() const { 
    return design_;
  }
  inline const VecInt2D& GetDesignRef() const {
    return design_;
  }
  inline double GetCriterion() const {
    return criteria_.GetCriterion();
  }
  inline void SwapInCol(int col, int r1, int r2) {
    criteria_.SwapInCol(col, r1, r2);
    std::swap(design_[r1][col], design_[r2][col]);
  }
  inline double PreSwapInCol(int col, int r1, int r2) const {
    return criteria_.PreSwapInCol(col, r1, r2);
  }
  void InitCriteria(const std::string& criteria_str);
  void Display() const;
  std::string GetCriteriaLogString() const;
 private:
  VecInt2D design_;
  int n_run_;
  int k_var_;
  MultiCriteria criteria_;
};
} // namespace LHD