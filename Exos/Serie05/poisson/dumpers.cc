/* -------------------------------------------------------------------------- */
#include "dumpers.hh"
#include "grid.hh"
/* -------------------------------------------------------------------------- */
#include <iomanip>
#include <sstream>
#include <fstream>
#include <array>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
void Dumper::set_min(float min) {
  m_min = min;
}

void Dumper::set_max(float max) {
  m_max = max;
}

float Dumper::min() const {
  return m_min;
}

float Dumper::max() const {
  return m_max;
}

/* -------------------------------------------------------------------------- */
void DumperASCII::dump(int step) {
  std::ofstream fout;
  std::stringstream sfilename;

  sfilename << "output/out_" << std::setfill('0') << std::setw(5) << step << ".pgm";
  fout.open(sfilename.str());

  int m = m_grid.m();
  int n = m_grid.n();

  fout <<  "P2" << std::endl << "# CREATOR: Poisson program" << std::endl;
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
  std::ofstream fout;
  std::stringstream sfilename;

  sfilename << "out_" << std::setfill('0') << std::setw(5) << step << ".bmp";
  fout.open(sfilename.str(), std::ios_base::binary);

  int h = m_grid.m();
  int w = m_grid.n();

  int row_size = 3 * w;
  // if the file width (3*w) is not a multiple of 4 adds enough bytes to make it
  // a multiple of 4
  int padding = (4 - (row_size) % 4) % 4;
  row_size += padding;

  int filesize = 54 + (row_size)*h;

  std::vector<char> img(row_size*h);
  std::fill(img.begin(), img.end(), 0);


  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      float v = ((m_grid(h - 1 - i, j) - m_min) / (m_max - m_min));

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

  bmpinfoheader[4]  = w;
  bmpinfoheader[5]  = w >> 8;
  bmpinfoheader[6]  = w >> 16;
  bmpinfoheader[7]  = w >> 24;
  bmpinfoheader[8]  = h;
  bmpinfoheader[9]  = h >> 8;
  bmpinfoheader[10] = h >> 16;
  bmpinfoheader[11] = h >> 24;
  bmpinfoheader[20] = (filesize - 54);
  bmpinfoheader[21] = (filesize - 54) >> 8;
  bmpinfoheader[22] = (filesize - 54) >> 16;
  bmpinfoheader[23] = (filesize - 54) >> 24;

  fout.write(bmpfileheader.data(), 14);
  fout.write(bmpinfoheader.data(), 40);

  fout.write(img.data(), h * row_size);
}
