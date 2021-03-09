#ifndef SIMULATION_HH
#define SIMULATION_HH

/* -------------------------------------------------------------------------- */
#include "double_buffer.hh"
#include "dumpers.hh"
#include "grid.hh"
/* -------------------------------------------------------------------------- */
#include <tuple>
#include <memory>
#include <array>
#include <mpi.h>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
class Simulation {
public:
  Simulation(int m, int n, MPI_Comm communicator);

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

  /// compute one row
  float compute_row(int i);

private:
  /// Global problem size
  int m_global_m, m_global_n;

  /// Local problem size
  int m_local_m, m_local_n;

  /// local offset
  int m_offset_m, m_offset_n;

  /// local neighbors
  int m_north_prank, m_south_prank;

  /// proc rank
  int m_prank;

  /// communicator size
  int m_psize;

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

  /// communicator
  MPI_Comm m_communicator;
};

#endif /* SIMULATION_HH */
