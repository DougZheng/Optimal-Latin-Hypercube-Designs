#include "multicriteria.h"

namespace LHD {
MultiCriteria& MultiCriteria::operator=(const MultiCriteria& o) {
  weights_ = o.weights_;
  criteria_ = o.criteria_;
  for (auto& c : criteria_) {
    c = c->Clone();
  }
  return *this;
}

double MultiCriteria::GetCriterion() const {
  double val = 0;
  for (int i = 0; i < criteria_.size(); ++i) {
    double vi = criteria_[i]->GetCriterion();
    auto bound = criteria_[i]->GetCriterionBound();
    double norm_vi = (vi - bound.first) / (bound.second - bound.first);
    val += weights_[i] * norm_vi;
  }
  return val;
}

void MultiCriteria::SwapInCol(int col, int r1, int r2) {
  for (auto& c : criteria_) {
    c->SwapInCol(col, r1, r2);
  }
}

double MultiCriteria::PreSwapInCol(int col, int r1, int r2) const {
  double val = 0;
  for (int i = 0; i < criteria_.size(); ++i) {
    double vi = criteria_[i]->PreSwapInCol(col, r1, r2);
    auto bound = criteria_[i]->GetCriterionBound();
    double norm_vi = (vi - bound.first) / (bound.second - bound.first);
    val += weights_[i] * norm_vi;
  }
  return val;
}

void MultiCriteria::SetDesign(const std::vector<std::vector<int>>* design) {
  for (int i = 0; i < criteria_.size(); ++i) {
    criteria_[i]->SetDesign(design);
  }
}

void MultiCriteria::AddCritirion(std::shared_ptr<Criterion> criterion, double weight) {
  criteria_.push_back(criterion);
  weights_.push_back(weight);
}

std::string MultiCriteria::GetOneLineLogString() const {
  std::string s;
  for (int i = 0; i < criteria_.size(); ++i) {
    s += criteria_[i]->GetCriterionName() +
      ": " + std::to_string(criteria_[i]->GetCriterion()) + ", ";
  }
  s += "Val: " + std::to_string(GetCriterion());
  return s;
}

std::string MultiCriteria::GetDescriptionString() const {
  std::string s;
  for (int i = 0; i < criteria_.size(); ++i) {
    s += criteria_[i]->GetCriterionName() + 
      ": " + std::to_string(criteria_[i]->GetCriterion()) +
      " (weight = " + std::to_string(weights_[i]) + ")" + "\n";
  }
  return s;
}

void MultiCriteria::Clear() {
  weights_.clear();
  criteria_.clear();
}
} // namespace LHD