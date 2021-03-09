#include <vector>
#include <iostream>

int main() {
  constexpr size_t N = 1000;
  std::vector<double> data(N);

  for(size_t i = 0; i < N; ++i) {
     data[i] = i;
  }

  double sum = 0.;
  for(size_t i = 0; i <= N; ++i) {
    sum += data[i];
  }

  std::cout << (N * (N-1) / 2.) << " == " << sum << std::endl;
}
