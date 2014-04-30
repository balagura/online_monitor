// Author V. Balagura, balagura@cern.ch (19.11.2012)

#include <raw.hh>
#include <deque>
#include <cstring> // for strlen()
#include <algorithm> // for reverse_copy

bool ReadSpill::find_spill_start() {
  // find the pattern (any 4 bytes + "SPIL  ") marking the spill start in the file, stop just after it, return true,
  // 4 bytes contain the acquisition number, store it in the "m_acquisition_number" variable;
  // if the pattern is not found: return false;
  deque<unsigned short int> pattern;
  unsigned short int i;
  m_data.clear(); // empty previous data
  m_sca.clear();
  static int pattern_length = strlen("SPIL  ")/2 + 2;
  while (m_file.read((char*)&i, 2)) {
    // printf("%x ",i);
    pattern.push_back(i); // pattern stores last bytes from the file
    if (pattern.size()  > pattern_length) pattern.pop_front();
    if (pattern.size() == pattern_length &&
	pattern[2] == ('S' | (((unsigned short int) 'P') << 8)) && // note swapped bytes
	pattern[3] == ('I' | (((unsigned short int) 'L') << 8)) &&
	pattern[4] == (' ' | (((unsigned short int) ' ') << 8))) {
      m_acquisition_number = ((((unsigned int) pattern[0])<<16) | pattern[1]);
      // printf("\n%x %x %x %x %x\n", pattern[0],pattern[1],pattern[2],pattern[3],pattern[4]);
      // printf("Acquisition number: %i\n", m_acquisition_number);
      return true;
    }
  }
  return false;
}
bool ReadSpill::find_spill_end() {
  // find the pattern ("    " + 0xffff) marking the spill end in the file, stop just after it, return true,
  // if the pattern is not found: return false;
  static int end_pattern_length = strlen("    ")/2 + 1;
  static int start_pattern_length = strlen("SPIL  ")/2 + 2; // note, there may be spill start without spill end
  unsigned short int i;
  // printf("search spill end");
  while (m_file.read((char*)&i, 2)) {
    // printf("%x ",i);
    m_data.push_back(i); // store everything in m_data
    unsigned int len = m_data.size();
    // first, check there is no new spill start
    if (len >= start_pattern_length &&
	m_data[len-3] == ('S' | (((unsigned short int) 'P') << 8)) && // note swapped bytes
	m_data[len-2] == ('I' | (((unsigned short int) 'L') << 8)) &&
	m_data[len-1] == (' ' | (((unsigned short int) ' ') << 8))) {
      m_acquisition_number = ((((unsigned int) m_data[len-5])<<16) | m_data[len-4]);
      // printf("\nNew spill start found %x %x %x %x %x\n", m_data[len-5],m_data[len-4],m_data[len-3],m_data[len-2],m_data[len-1]);
      // printf("Acquisition number: %i\n", m_acquisition_number);
      m_data.clear(); len = 0; // reset everything, assign new acqusition number
    }
    if (len >= end_pattern_length &&
	m_data[len-3] == 0x2020 && m_data[len-2] == 0x2020 && m_data[len-1] == 0xffff) { // 0x20 == ' ' (space), note also,
      // bytes are swapped in two-bytes unsigned short int, but here bytes are the same inside 3 pairs: ' ' or 0xff
      m_data.resize(len-3); // remove match pattern in the end
      return true;
    }
  }
  return false;
}
void ReadSpill::find_chip_bounds() {
  static int start_pattern_length = 2 + strlen("CHIP  ")/2; // pattern 0xfd 0xff [0x01-\0x04] \0xff "CHIP  " marks start of CHIP data block
  static int end_pattern_length = 3;    // pattern (any 2 bytes) 0xfe 0xff [[\x01-\x04]\xff marks end of CHIP data block
  m_chip.clear();
  data_iter start = m_data.begin(), end = m_data.end();
  bool found;
  do {
    found = false;
    for ( ; start + start_pattern_length <= end; ++start) {
      if (*start == 0xfffd &&
	  0xff01 <= *(start+1) && *(start+1) <= 0xff04 &&
	  *(start+2) == ('C' | (((unsigned short int) 'H') << 8)) && // note swapped bytes
	  *(start+3) == ('I' | (((unsigned short int) 'P') << 8)) &&
	  *(start+4) == (' ' | (((unsigned short int) ' ') << 8))) { start += start_pattern_length; found = true; break; }
    }
    if (!found) break;
    found = false;
    for (data_iter finish = start; finish + end_pattern_length <= end; ++finish) {
      if (*(finish+1) == 0xfffe &&
	  0xff01 <= *(finish+2) && *(finish+2) <= 0xff04 &&
	  start <= finish-1 && *(finish-1) < 256) { // in addition to pattern match, require chip id == last entry in data block < 256
	ChipBounds cb; cb.begin = start; cb.end = finish;
	m_chip.push_back(cb);
	start =  finish + end_pattern_length;
	found = true;
	break;
      }
    }
  } while (found);
}
void ReadSpill::fill_sca() {
  ChipSCA s;
  for (int i=0; i<int(m_chip.size()); ++i) {
    data_iter start = m_chip[i].begin, end = m_chip[i].end;
    // printf("end-start-1 = %f, chip = %x\n",double(end-start-1)/129., *(end-1));
    s.chip_id = *(end - 1);
    int nSCA = (end - start - 1) / 129; // chip id takes one word in the end, every SCA takes 64*2 (for ADC) + 1 (for BXID)
    for (int i=0; i<nSCA; ++i) {
      s.isca = nSCA - 1 - i; // last SCA is stored first
      reverse_copy(start + 128*i,      start + 128*i +  64, s.high.begin()); // copy uses ChipADC(unsigned short int d)
      reverse_copy(start + 128*i + 64, start + 128*i + 128, s.low.begin());  // reverse_copy since last channel comes first in data
      unsigned short int bxid = *(start + 128*nSCA + i); // last bxid is also first
      m_sca[bxid].push_back(s);
    }
  }
}
bool ReadSpill::next() {
  if (find_spill_start() & find_spill_end()) {
    find_chip_bounds();
    fill_sca();
    return true;
  } else return false;
}
