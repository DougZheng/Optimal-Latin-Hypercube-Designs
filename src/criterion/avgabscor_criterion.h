#pragma once

#include "criterion.h"

namespace LHD {
class AvgAbsCorCriterion : public Criterion {
 public:
  AvgAbsCorCriterion(const std::vector<std::vector<int>>* design);
  virtual ~AvgAbsCorCriterion() = default;
  inline std::shared_ptr<Criterion> Clone() const override {
    return std::make_shared<AvgAbsCorCriterion>(*this);
  }
  inline double GetCriterion() const override {
    return rho_avg_;
  }
  inline std::string GetCriterionName() const override {
    return "AvgAbsCor";
  }
  inline std::pair<double, double> GetCriterionBound() const override {
    return {0, 1};
  }
  void SwapInCol(int col, int r1, int r2) override;
  double PreSwapInCol(int col, int r1, int r2) const override;
 private:
  double kCorDenominator;
  double kCorPairNum;
  std::vector<std::vector<double>> cor_;
  double rho_avg_;
};
} // namespace LHD