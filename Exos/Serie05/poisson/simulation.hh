#ifndef SIMULATION_HH
#define SIMULATION_HH

/* -------------------------------------------------------------------------- */
#include "double_buffer.hh"
#include "dumpers.hh"
#include "grid.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
class Simulation {
public:
  Simulation(int m, int n);

  /// set the initial conditions, Dirichlet and source term
  virtual void set_initial_conditions();

  /// perform the simulation
  std::tuple<float, int> compute();

  /// access the precision
  void set_epsilon(float epsilon);
  float epsilon() const;

protected:
  /// compute one step and return an error
  virtual float compute_step();

private:
  /// Global problem size
  int m_global_m, m_global_n;

  /// Precision to achieve
  float m_epsilon;

  /// grid spacing
  float m_h_m;
  float m_h_n;

  /// Grids storage
  DoubleBuffer m_grids;

  /// source term
  Grid m_f;

  /// Dumper to use for outputs
  std::unique_ptr<Dumper> m_dumper;
};

#endif /* SIMULATION_HH */
