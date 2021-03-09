#include <mpi.h>

#ifndef DUMPERS_HH
#define DUMPERS_HH

/* -------------------------------------------------------------------------- */
class Grid;

/* -------------------------------------------------------------------------- */
class Dumper {
public:
  explicit Dumper(const Grid & grid, MPI_Comm communicator)
      : m_grid(grid), m_min(-1.), m_max(1.), m_communicator(communicator) {}

  virtual void dump(int step) = 0;

  void set_min(float min);
  void set_max(float max);

  float min() const;
  float max() const;

protected:
  const Grid & m_grid;
  float m_min, m_max;

  MPI_Comm m_communicator;
};

/* -------------------------------------------------------------------------- */
class DumperASCII : public Dumper {
public:
  explicit DumperASCII(const Grid & grid, MPI_Comm communicator)
      : Dumper(grid, communicator) {}

  virtual void dump(int step);
};

/* -------------------------------------------------------------------------- */
class DumperBinary : public Dumper {
public:
  explicit DumperBinary(const Grid & grid, MPI_Comm communicator)
      : Dumper(grid, communicator) {}

  virtual void dump(int step);
};

#endif /* DUMPERS_HH */
