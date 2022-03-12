#include "SA.h"

int main(int argc, char** argv) {
  if (argc < 5) {
    std::cout << "Usage: ${exe} ${n} ${k} ${w} ${cnt} [${file}]" << std::endl;
    return 1;
  }
  // ./SA 20 8 0.5 10 >LHD
  // ./SA 20 8 0.5 10 LHD
  int n = std::stoi(argv[1]);
  int k = std::stoi(argv[2]);
  double w = std::stod(argv[3]);
  int cnt = std::stoi(argv[4]);
  auto searcher = SearchingAlgorithm::SA(n, k);
  searcher.SetW(w);
  searcher.SetRepeatedCnt(cnt);
  if (argc > 5) {
    auto design = SearchingAlgorithm::Design::ReadDesign(argv[5]);
    searcher.IncrementalSearch(design).Display();
  } else {
    searcher.Search().Display();
  }
  return 0;
}