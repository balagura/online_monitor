// Author V. Balagura, balagura@cern.ch (07.12.2015)

#include <raw.hh>
#include <fstream>
#include <cstring> // for strlen
#include <cstdlib> // for atoi
#include <sys/stat.h>

class SpillStart : public ReadSpill {
public:
  SpillStart() {}
  SpillStart(istream& is) : ReadSpill(is) {}

  streampos find_spill_start() {
    static int pattern_length = strlen("SPIL  ")/2 + 3;
    bool ok = ReadSpill::find_spill_start();
    if (!ok) return -1;
    return m_input->tellg() - streampos(2*pattern_length);
  }
};
  
int main(int argc, char** argv) {
  if (argc != 3) {
    cout << "Usage: " << argv[0] << "  <file_name>  <maximal_slice_size_in_bytes>\n"
      "The program splits the raw file into N = ceiling(file_size / <maximal_slice_size_in_bytes>) slices and returns two lists:\n"
      "acq_1 position_1\n"
      "acq_2 position_2\n"
      "...\n"
      "acq_N position_N\n\n"
      "where acq_i is the minimal spill (==acquision) number in slice i and\n"
      "position_i is the offset of the acq_i spill block in bytes from the beginning of the file\n"
      "(use C++ function \'seekg(position_i)\' to go there)\n";
    return 0;
  }

  char* filename = argv[1];
  long long int max_size = (atoll(argv[2])+1) / 2 * 2; // in bytes, ensure it is multiple of 2 as RAW file is written in 2-byte words,
  // shift by one byte completely spoils reading with ReadSpill.
  // in case of rounding, if argv[2]!=2*n, do not make more slices, ie. choose larger slice size 2n+1, 


  struct stat st;
  stat(filename, &st);
  int size = st.st_size;

  int N = 1 + int((size-1) / max_size); // == ceiling(size / max_size)
  //  cerr << "N=" << N << ", filesize=" << size << ", max_size=" << max_size << endl;
  
  ifstream file;
  file.open(filename, ios::in | ios::binary);

  for (long long int i=0; i<N; ++i) {
    file.seekg(i*max_size);
    SpillStart raw(file);
    streampos pos = raw.find_spill_start();
    cout << raw.acquisition_number() << " " << pos << endl;
  }
  return 0;
}
