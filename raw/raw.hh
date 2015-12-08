#ifndef RAW_HH
#define RAW_HH

// Author V. Balagura, balagura@cern.ch (19.11.2012)

#include <vector>
#include <map>
#include <iostream>

using namespace std;

struct ChipADC {
  bool gain, hit;
  unsigned short int adc;
  ChipADC() : gain(false), hit(false), adc(0) {}
  ChipADC(unsigned short int d) { gain = (d >> 13) & 1; hit = (d >> 12) & 1; adc = (d & 0xfff); }
};
struct ChipSCA {
  vector<ChipADC> high, low;
  unsigned short int chip_id;
  unsigned short int isca;
  ChipSCA() : chip_id(0), isca(0) { high.resize(64); low.resize(64); }
};
struct ReadSpill {
  void open(istream& is) { m_input = &is; m_data.reserve(1024); }

  ReadSpill() {}
  ReadSpill(istream& is) { open(is); }
  ReadSpill(const ReadSpill&) {} // assumes no copies are used; needed for vector<ReadSpill>

  bool next(); // The main function:
  // finds next start/end spill marks in the file,
  // if successful, stores acquisition_number and binary data between the spill marks in the file (see data()),
  // splits this data into chip blocks and stores the positions of boundaries (see chip()),
  // then converts every block into ChipSCA structure above, grouping them according to BX ID,
  // ie. creates associative array of (key-value) pairs (see sca()), with key == BX ID and value = vector<ChipSCA>;
  // returns true.
  // Otherwise (eg. at EOF) returns false
  // If user is interested only in final data, they may be obtained with sca() function, the rest is not needed.

  const map<unsigned short int, vector<ChipSCA> >& sca() const { return m_sca; } // associative array of (key-value) pairs with key == BX ID and
  // value = vector<ChipSCA>
  unsigned int acquisition_number() const { return m_acquisition_number; }
  typedef vector<unsigned short int>::const_iterator data_iter;
  struct ChipBounds { data_iter begin, end; };
  const vector<ChipBounds>& chip() const { return m_chip; } // structure keeping chip block boundaries in data() below
  const vector<unsigned short int>& data() const { return m_data; } // contains all spill data
  ~ReadSpill() {}

protected:
  istream* m_input;
  unsigned int m_acquisition_number;
  vector<unsigned short int> m_data;
  vector<ChipBounds> m_chip;
  map<unsigned short int, vector<ChipSCA> > m_sca; // m_sca[bxid] = vector<ChipSCA> structures
  bool find_spill_start();
  bool find_spill_end();
  void find_chip_bounds();
  void fill_sca();
};

#endif
