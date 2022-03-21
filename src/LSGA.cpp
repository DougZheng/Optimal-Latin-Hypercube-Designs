#include "LSGA.h"

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cout << "Usage: ${exe} ${n} ${k} ${w} ${cnt}" << std::endl;
    return 1;
  }
  // ./LSGA 20 8 0.5 1000 >LHD
  int n = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  double w = std::stod(argv[3]);
  int cnt = std::stoi(argv[4]);
  auto searcher = SearchingAlgorithm::LSGA(n, k);
  searcher.SetW(w);
  searcher.SetIterateCnt(cnt);
  searcher.Search().Display();
  return 0;
}