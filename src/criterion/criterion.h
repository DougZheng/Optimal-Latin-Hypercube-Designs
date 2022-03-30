#pragma once

#include <cassert>
#include <string>
#include <vector>
#include <utility>
#include <memory>

namespace LHD {
class Criterion {
 public:
  explicit Criterion(const std::vector<std::vector<int>>* design) : design_(design) {
    assert(design->size() > 0);
    n_run_ = design->size();
    k_var_ = (*design)[0].size();
  }
  virtual ~Criterion() = default;
  virtual inline void SetDesign(const std::vector<std::vector<int>>* design) {
    design_ = design;
  }
  virtual inline std::shared_ptr<Criterion> Clone() const = 0;
  virtual inline double GetCriterion() const = 0;
  virtual inline std::string GetCriterionName() const = 0;
  virtual inline std::pair<double, double> GetCriterionBound() const = 0;
  virtual void SwapInCol(int col, int r1, int r2) = 0;
  virtual double PreSwapInCol(int col, int r1, int r2) const = 0;
 protected:
  const std::vector<std::vector<int>>* design_;
  int n_run_;
  int k_var_;
};
} // namespace LHD