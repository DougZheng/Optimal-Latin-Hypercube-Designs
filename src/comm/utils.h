#pragma once

#include <cassert>
#include <string>
#include <vector>
#include <random>

namespace Utils {
// "1.0PhiP15L1";
static std::string SPhiPL1(double weight, int power = 15) { 
  return std::to_string(weight) + "PhiP" + std::to_string(power) + "L1";
}
// "1.0PhiP15L2";
static std::string SPhiPL2(double weight, int power = 15) { 
  return std::to_string(weight) + "PhiP" + std::to_string(power) + "L2";
}
// "1.0MaxAbsCor";
static std::string SMaxAbsCor(double weight) { 
  return std::to_string(weight) + "MaxAbsCor";
}
// "1.0AvgAbsCor";
static std::string SAvgAbsCor(double weight) { 
  return std::to_string(weight) + "AvgAbsCor";
}
// "0.5PhiP15L1+0.5MaxAbsCor";
static std::string SMultiCriteria(const std::string& s1, const std::string& s2) {
  return s1 + "+" + s2;
}

static std::vector<std::string> Split(const std::string& s, const std::string& delimiters = ",") {
  std::vector<std::string> tokens;
  auto last_pos = s.find_first_not_of(delimiters, 0);
  auto pos = s.find_first_of(delimiters, last_pos);
  while (pos != std::string::npos || last_pos != std::string::npos) {
    tokens.emplace_back(s.substr(last_pos, pos - last_pos));
    last_pos = s.find_first_not_of(delimiters, pos);
    pos = s.find_first_of(delimiters, last_pos);
  }
  return tokens;
}

static inline double QuickPow15(double x) {
  double y = x;
  y *= y, y *= y, y *= y, y *= y;
  return y / x;
}

static inline double QuickPow(double x, int p) {
  assert(p >= 0);
  if (p == 15) return QuickPow15(x);
  double y = 1;
  while (p > 0) {
    if (p & 1) y *= x;
    x *= x;
    p >>= 1;
  }
  return y;
}

template<typename T>
static void ShuffleM(std::vector<T>& list, int m, std::mt19937& rng) {
  int n = list.size();
  assert(m <= n);
  for (int i = 0; i < m; ++i) {
    int j = rng() % (n - i);
    std::swap(list[i], list[j]);
  }
}
} // namespace Utils