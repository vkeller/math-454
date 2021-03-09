
/* -------------------------------------------------------------------------- */
#include "dumpers.hh"
#include "grid.hh"
/* -------------------------------------------------------------------------- */
#include <array>
#include <fstream>
#include <iomanip>
#include <sstream>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void Dumper::set_min(float min) { m_min = min; }

void Dumper::set_max(float max) { m_max = max; }

float Dumper::min() const { return m_min; }

float Dumper::max() const { return m_max; }

/* -------------------------------------------------------------------------- */
void DumperASCII::dump(int step) {
  std::ofstream fout;
  std::stringstream sfilename;

  sfilename << "out_" << std::setfill('0') << std::setw(5) << step << ".pgm";
  fout.open(sfilename.str());

  int m = m_grid.m();
  int n = m_grid.n();

  fout << "P2" << std::endl << "# CREATOR: Poisson program" << std::endl;
  fout << m << " " << n << std::endl;
  fout << 255 << std::endl;

  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      int v = 255. * (m_grid(i, j) - m_min) / (m_max - m_min);
      v = std::min(v, 255);
      fout << v << std::endl;
    }
  }
}

/* -------------------------------------------------------------------------- */
void DumperBinary::dump(int step) {
  std::stringstream sfilename;

  int prank, psize;
  MPI_Comm_rank(m_communicator, &prank);
  MPI_Comm_size(m_communicator, &psize);

  sfilename << "out_" << std::setfill('0') << std::setw(5) << step << ".bmp";

  int h = m_grid.m();
  int w = m_grid.n();

  int start = 0;

  // removing the ghosts from the height
  if (psize > 1) {
    h = (prank == 0 || prank == psize - 1 ? h - 1 : h - 2);
    start = (prank == 0 ? 0 : 1);
  }

  // Gathering the size of every processors, this could be done as in the
  // constructor of the Simulation instead
  std::vector<int> size_per_proc(psize);
  MPI_Allgather(&h, 1, MPI_INT, size_per_proc.data(), 1, MPI_INT,
                m_communicator);

  // determining the local offset
  int offset_h = 0;
  for (int i = 0; i < prank; ++i) {
    offset_h += size_per_proc[i];
  }

  int total_h = offset_h;
  for (int i = prank; i < psize; ++i) {
    total_h += size_per_proc[i];
  }

  int row_size = 3 * w;
  // if the file width (3*w) is not a multiple of 4 adds enough bytes to make it
  // a multiple of 4
  int padding = (4 - (row_size) % 4) % 4;
  row_size += padding;

  int filesize = 54 + row_size * total_h;

  std::vector<char> img(row_size * h);
  std::fill(img.begin(), img.end(), 0);

  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      float v = ((m_grid(i + start, j) - m_min) / (m_max - m_min));

      float r = v * 255; // Red channel
      float g = v * 255; // Green channel
      float b = v * 255; // Red channel

      r = std::min(r, 255.f);
      g = std::min(g, 255.f);
      b = std::min(b, 255.f);

      img[row_size * i + 3 * j + 2] = r;
      img[row_size * i + 3 * j + 1] = g;
      img[row_size * i + 3 * j + 0] = b;
    }
  }

  std::array<char, 14> bmpfileheader = {'B', 'M', 0, 0,  0, 0, 0,
                                        0,   0,   0, 54, 0, 0, 0};
  std::array<char, 40> bmpinfoheader = {40, 0, 0, 0, 0, 0, 0,  0,
                                        0,  0, 0, 0, 1, 0, 24, 0};

  bmpfileheader[2] = filesize;
  bmpfileheader[3] = filesize >> 8;
  bmpfileheader[4] = filesize >> 16;
  bmpfileheader[5] = filesize >> 24;

  bmpinfoheader[4] = w;
  bmpinfoheader[5] = w >> 8;
  bmpinfoheader[6] = w >> 16;
  bmpinfoheader[7] = w >> 24;
  bmpinfoheader[8] = total_h;
  bmpinfoheader[9] = total_h >> 8;
  bmpinfoheader[10] = total_h >> 16;
  bmpinfoheader[11] = total_h >> 24;
  bmpinfoheader[20] = (filesize - 54);
  bmpinfoheader[21] = (filesize - 54) >> 8;
  bmpinfoheader[22] = (filesize - 54) >> 16;
  bmpinfoheader[23] = (filesize - 54) >> 24;

  MPI_File fh;
  MPI_Status status;

  // opening a file in write and create mode
  MPI_File_open(MPI_COMM_WORLD, sfilename.str().c_str(),
                MPI_MODE_WRONLY | MPI_MODE_CREATE,
                MPI_INFO_NULL, &fh);
  // defining the size of the file
  MPI_File_set_size(fh, filesize);

  // rank 0 writes the header
  if (prank == 0) {
    MPI_File_write_at(fh, 0, bmpfileheader.data(), 14, MPI_CHAR, &status);
    MPI_File_write_at(fh, 14, bmpinfoheader.data(), 40, MPI_CHAR, &status);
  }

  int offset = 54 + row_size * offset_h;

  // We also could write that data with a write_at, the view is just to show
  // different possibilities
  MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, "native",
                      MPI_INFO_NULL);

  MPI_File_write(fh, img.data(), img.size(), MPI_CHAR, &status);
  MPI_File_close(&fh);
}
