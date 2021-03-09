/* -------------------------------------------------------------------------- */
#include "simulation.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
Simulation::Simulation(int m, int n)
    : m_global_m(m), m_global_n(n), m_epsilon(1e-7), m_h_m(1. / m),
      m_h_n(1. / n), m_grids(m, n), m_f(m, n),
      m_dumper(new DumperASCII(m_grids.old())) {}

/* -------------------------------------------------------------------------- */
void Simulation::set_initial_conditions() {
  for (int i = 0; i < m_global_m; i++) {
    for (int j = 0; j < m_global_n; j++) {
      m_f(i, j) = -2. * 100. * M_PI * M_PI * std::sin(10. * M_PI * i * m_h_m) *
                  std::sin(10. * M_PI * j * m_h_n);
    }
  }
}

/* -------------------------------------------------------------------------- */
std::tuple<float, int> Simulation::compute() {
  int s = 0;
  float l2 = 0;
  do {
    l2 = compute_step();

    m_grids.swap();

    m_dumper->dump(s);
    ++s;
  } while (l2 > m_epsilon);

  return std::make_tuple(l2, s);
}

/* -------------------------------------------------------------------------- */
void Simulation::set_epsilon(float epsilon) { m_epsilon = epsilon; }

/* -------------------------------------------------------------------------- */
float Simulation::epsilon() const { return m_epsilon; }

/* -------------------------------------------------------------------------- */
float Simulation::compute_step() {
  float l2 = 0.;

  Grid & u = m_grids.current();
  Grid & uo = m_grids.old();

  for (int i = 1; i < m_global_m - 1; i++) {
    for (int j = 1; j < m_global_n - 1; j++) {
      // computation of the new step
      u(i, j) = 0.25 * (uo(i - 1, j) + uo(i + 1, j) + uo(i, j - 1) +
                        uo(i, j + 1) - m_f(i, j) * m_h_m * m_h_n);

      // L2 norm
      l2 += (uo(i, j) - u(i, j)) * (uo(i, j) - u(i, j));
    }
  }

  return l2;
}
