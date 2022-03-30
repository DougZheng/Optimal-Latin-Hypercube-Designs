#pragma once

#include "criterion.h"

namespace LHD {
class PhiPL1Criterion : public Criterion {
 public:
  PhiPL1Criterion(const std::vector<std::vector<int>>* design, int power = 15);
  virtual ~PhiPL1Criterion() = default;
  inline std::shared_ptr<Criterion> Clone() const override {
    return std::make_shared<PhiPL1Criterion>(*this);
  }
  inline double GetCriterion() const override {
    return phi_p_;
  }
  inline std::string GetCriterionName() const override {
    return "PhiP" + std::to_string(kPower) + "L1";
  }
  inline std::pair<double, double> GetCriterionBound() const override;
  void SwapInCol(int col, int r1, int r2) override;
  double PreSwapInCol(int col, int r1, int r2) const override;
 private:
  int kPower;
  std::vector<std::vector<int>> dis_;
  double phi_p_;
};
} // namespace LHD