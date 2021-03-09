/* -------------------------------------------------------------------------- */
#include "simulation.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
Simulation::Simulation(int m, int n, MPI_Comm communicator)
    : m_global_m(m), m_global_n(n), m_epsilon(1e-7), m_h_m(1. / m),
      m_h_n(1. / n), m_dumper(new DumperBinary(m_grids.old(), communicator)),
      m_communicator(communicator) {

  // retrieving the number of proc and the rank in the proc pool
  MPI_Comm_rank(m_communicator, &m_prank);
  MPI_Comm_size(m_communicator, &m_psize);

  // computation of the local size of the grid the remainder is spread equally
  // on the first processors
  m_local_m = m / m_psize + (m_prank < m % m_psize ? 1 : 0);
  m_local_n = n;

  // adding the ghosts lines if needed
  if (m_psize > 1)
    m_local_m += (m_prank == 0 || m_prank == m_psize - 1) ? 1 : 2;

  // computing the offsets of the local grid in the global one
  m_offset_m =
      (m / m_psize) * m_prank + (m_prank < m % m_psize ? m_prank : m % m_psize);
  m_offset_n = 0;

  // resizing the different grids
  m_grids.resize(m_local_m, m_local_n);
  m_f.resize(m_local_m, m_local_n);

  // determining the rank of the neighbors
  m_north_prank = (m_prank == 0 ? MPI_PROC_NULL : m_prank - 1);
  m_south_prank = (m_prank == (m_psize - 1) ? MPI_PROC_NULL : m_prank + 1);

  // Some info if needed to debug
  // std::cout << m_prank << " " << m_global_m << " " << m_global_n << " "
  //           << m_local_m << " " << m_local_n << " " << m_offset_m << " "
  //           << m_offset_n << " " << m_north_prank << " " << m_south_prank
  //           << std::endl;
}

/* -------------------------------------------------------------------------- */
void Simulation::set_initial_conditions() {
  int i_start = (m_prank == 0 ? 0 : 1);
  int i_end   = (m_prank == m_psize - 1 ? m_local_m : m_local_m - 1);

  for (int i = i_start; i < i_end; i++) {
    for (int j = 0; j < m_local_n; j++) {
      m_f(i, j) = -2. * 100. * M_PI * M_PI *
                  std::sin(10. * M_PI * (m_offset_m + i - i_start) * m_h_m) *
                  std::sin(10. * M_PI * (m_offset_n + j) * m_h_n);
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

    // m_dumper->dump(s);
    ++s;
  } while (l2 > m_epsilon);

  return std::make_tuple(l2, s);
}

/* -------------------------------------------------------------------------- */
void Simulation::set_epsilon(float epsilon) { m_epsilon = epsilon; }

/* -------------------------------------------------------------------------- */
float Simulation::epsilon() const { return m_epsilon; }

/* -------------------------------------------------------------------------- */
inline float Simulation::compute_row(int i) {
  float l2 = 0;

  Grid & u = m_grids.current();
  Grid & uo = m_grids.old();

  for (int j = 1; j < m_local_n - 1; j++) {
    // computation of the new step
    u(i, j) = 0.25 * (uo(i - 1, j) + uo(i + 1, j) + uo(i, j - 1) +
                      uo(i, j + 1) - m_f(i, j) * m_h_m * m_h_n);

    // L2 norm
    l2 += (uo(i, j) - u(i, j)) * (uo(i, j) - u(i, j));
  }

  return l2;
}

/* -------------------------------------------------------------------------- */
float Simulation::compute_step() {
  float l2 = 0.;

  Grid & uo = m_grids.old();

  // Taking care of communications going up (so receiving from bottom)
  MPI_Sendrecv(&uo(1, 0), m_local_n, MPI_FLOAT, m_north_prank, 0,
               &uo(m_local_m - 1, 0), m_local_n, MPI_FLOAT, m_south_prank, 0,
               m_communicator, MPI_STATUS_IGNORE);

    // Taking care of communications going down (so receiving from top)
  MPI_Sendrecv(&uo(m_local_m - 2, 0), m_local_n, MPI_FLOAT, m_south_prank, 0,
               &uo(0, 0), m_local_n, MPI_FLOAT, m_north_prank, 0,
               m_communicator, MPI_STATUS_IGNORE);


  // computing all the rows
  for (int i = 1; i < m_local_m - 1; i++) {
    l2 += compute_row(i);
  }

  // Summing the value of all the processors together
  MPI_Allreduce(MPI_IN_PLACE, &l2, 1, MPI_FLOAT, MPI_SUM, m_communicator);

  return l2;
}
