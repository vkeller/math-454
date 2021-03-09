#ifndef GRID_HH
#define GRID_HH

/* -------------------------------------------------------------------------- */
#include <vector>
/* -------------------------------------------------------------------------- */

class Grid {
public:
  Grid(int m, int n);

  /// access the value [i][j] of the grid
  inline float & operator()(int i, int j) { return m_storage[i * m_n + j]; }
  inline const float & operator()(int i, int j) const {
    return m_storage[i * m_n + j];
  }

  /// set the grid to 0
  void clear();

  int m() const;
  int n() const;

private:
  int m_m, m_n;
  std::vector<float> m_storage;
};

#endif /* GRID_HH */
