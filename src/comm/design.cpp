#include "design.h"

#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <random>
#include <memory>

#include "criterion/phipl1_criterion.h"
#include "criterion/phipl2_criterion.h"
#include "criterion/maxabscor_criterion.h"
#include "criterion/avgabscor_criterion.h"

namespace LHD {
std::vector<std::vector<int>> Design::ReadDesign(const std::string& file) {
  std::ifstream in(file);
  if (in.fail()) {
    std::cerr << "open " << file << " failed." << std::endl;
    return VecInt2D();
  }
  int n, k;
  in >> n >> k;
  VecInt2D a(n, std::vector<int>(k));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < k; ++j) {
      in >> a[i][j];
    }
  }
  return a;
}

Design::Design(const VecInt2D& design) : design_(design) {
  assert(design.size() > 0);
  n_run_ = design.size();
  k_var_ = design[0].size();
}

Design::Design(int n, int k, int seed) : n_run_(n), k_var_(k) {
  std::mt19937 rng_(seed);
  design_.resize(n, std::vector<int>(k));
  std::vector<int> rnd_list(n);
  std::iota(rnd_list.begin(), rnd_list.end(), 1);
  for (int i = 0; i < k; ++i) {
    std::shuffle(rnd_list.begin(), rnd_list.end(), rng_);
    for (int j = 0; j < n; ++j) {
      design_[j][i] = rnd_list[j];
    }
  }
}

Design& Design::operator=(const Design& design) {
  design_ = design.design_;
  n_run_ = design.n_run_;
  k_var_ = design.k_var_;
  criteria_ = design.criteria_;
  criteria_.SetDesign(&design_);
  return *this;
}

void Design::InitCriteria(const std::string& criteria_str) {
  criteria_.Clear();
  const auto& tokens = Utils::Split(criteria_str, "+");
  for (const auto& token : tokens) {
    std::size_t idx = 0;
    double weight = std::stof(token, &idx);
    std::shared_ptr<Criterion> criterion;
    switch (token[idx]) {
      case 'P' : {
        int power = std::stoi(token.substr(idx + 4));
        int dis_level = token.back() == '1' ? 1 : 2;
        if (dis_level == 1) {
          criterion = std::make_shared<PhiPL1Criterion>(&design_, power);
        } else {
          criterion = std::make_shared<PhiPL2Criterion>(&design_, power);
        }
        break;
      }
      case 'M' : {
        criterion = std::make_shared<MaxAbsCorCriterion>(&design_);
        break;
      }
      case 'A' : {
        criterion = std::make_shared<AvgAbsCorCriterion>(&design_);
        break;
      }
      default : break;
    }
    criteria_.AddCritirion(criterion, weight);
  }
}

void Design::Display() const {
  std::cout << n_run_ << " " << k_var_ << "\n";
  int w = std::floor(std::log10(n_run_)) + 1;
  for (int i = 0; i < n_run_; ++i) {
    for (int j = 0; j < k_var_; ++j) {
      std::cout << std::setw(w) << design_[i][j] << " \n"[j == k_var_ - 1];
    }
  }
  std::cout << criteria_.GetDescriptionString();
  std::cout.flush();
}

std::string Design::GetCriteriaLogString() const {
  return criteria_.GetOneLineLogString();
}
}