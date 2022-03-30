#pragma once

#include "criterion.h"

namespace LHD {
class MaxAbsCorCriterion : public Criterion {
 public:
  MaxAbsCorCriterion(const std::vector<std::vector<int>>* design);
  virtual ~MaxAbsCorCriterion() = default;
  inline std::shared_ptr<Criterion> Clone() const override {
    return std::make_shared<MaxAbsCorCriterion>(*this);
  }
  inline double GetCriterion() const override {
    return rho_max_;
  }
  inline std::string GetCriterionName() const override {
    return "MaxAbsCor";
  }
  inline std::pair<double, double> GetCriterionBound() const override {
    return {0, 1};
  }
  void SwapInCol(int col, int r1, int r2) override;
  double PreSwapInCol(int col, int r1, int r2) const override;
 private:
  double kCorDenominator;
  std::vector<std::vector<double>> cor_;
  double rho_max_;
  int rho_max_i_;
  int rho_max_j_;
};
} // namespace LHD