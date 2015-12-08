// Author V. Balagura, balagura@cern.ch (19.11.2012)

#include <raw.hh>
#include <fstream>
#include <cstdlib>
#include <string>
#include <climits> // for INT_MAX

// // taken from http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
// // to get command line options:

// #include <algorithm>
// char* getCmdOption(char ** begin, char ** end, const std::string & option) {
//     char ** itr = std::find(begin, end, option);
//     if (itr != end && ++itr != end) return *itr;
//     return 0;
// }
// bool cmdOptionExists(char** begin, char** end, const std::string& option) {
//     return std::find(begin, end, option) != end;
// }
// // ------------------------------------------------------------

template<class T> struct DeleteMe { // helper class to ensure that new()'ly created object is deleted in the end
  DeleteMe() : p(NULL) {}
  void assign(T* fP) { p = fP; }
  ~DeleteMe() { if (p != NULL) delete(p); }
protected:
  T* p;
};

int main(int argc, char** argv) {
  char* filename;
  int nMaxLines = INT_MAX;
  double pedestal_suppression = 0.05;
  string format = "";
  streampos skipNBytes = 0;
  long int acqEnd = 0;
  if (argc < 2 || argc > 6) {
    cout << "Usage: " << argv[0] << " <file_name> [max_N_lines]\n"
      "or " << argv[0] << " <file_name> max_N_lines high_gain_triggers [pedestal_suppression]\n"
      "or " << argv[0] << " <file_name> high_gain_triggers pedestal_suppression skipNBytes acqEnd\n\n"
      "<file_name> can be \'-\' for stdin (note: then skipNBytes can not be used and should be zero)\n"
      "If max_N_lines is absent or 0, all file is processed\n"
      "If pedestal_suppression is given, all not triggered data\n"
      "is postscaled by this factor (between 0 and 1, default 0.05)\n"
      "skipNBytes defines how many first input bytes to skip\n"
      "acqEnd selects only spill (==acquisition) numbers less than acqEnd (0 to process everything)\n";
    return 0;
  }
  else {
    filename = argv[1];
    if (argc < 6) {
      if (argc >= 3) {
	int n = atoi(argv[2]);
	if (n != 0) nMaxLines = n;
      }
      if (argc >= 4) format = argv[3];
      if (argc >= 5) pedestal_suppression = atof(argv[4]);
    }
    else {
      format               = argv[2];
      pedestal_suppression = atof(argv[3]);
      skipNBytes           = atoi(argv[4]);
      acqEnd               = atoi(argv[5]);
    }
  }
  ReadSpill spill;
  DeleteMe<ifstream> closeFileInTheEnd;
  istream* is;
  if (string(filename) == string("-"))
    is = &cin;
  else {
    ifstream* file = new(ifstream);
    closeFileInTheEnd.assign(file);
    file->open(filename, ios::in | ios::binary);
    is = file;
  }
  if (skipNBytes != 0) is->seekg(skipNBytes);
  spill.open(*is);
  int iLine=0;
  
  while (spill.next()) {
    if ( acqEnd != 0 & spill.acquisition_number() >= acqEnd ) return 0;
    const map<unsigned short int, vector<ChipSCA> >& sca = spill.sca();
    for (map<unsigned short int, vector<ChipSCA> >::const_iterator isca=sca.begin(); isca!=sca.end(); ++isca) {
      int bx  =  isca->first;
      for (int i=0; i<int(isca->second.size()); ++i) {
	const ChipSCA& s = isca->second[i];
	if (format == "") {
	  cout << spill.acquisition_number() << " " << bx << " " << s.chip_id << " " << s.isca+1;
	  for (int ich=0; ich<int(s.high.size()); ++ich) {
	    ChipADC ah =s.high[ich], al = s.low[ich];
	    /*          cout << " " << (ah.gain ? 'T' : 'F') << " " << (ah.hit ? 'T' : 'F') << " " << ah.adc
			<< " " << (al.gain ? 'T' : 'F') << " " << (al.hit ? 'T' : 'F') << " " << al.adc; */
	    cout << " " << (ah.hit ? 'T' : 'F') << " " << ah.adc
		 << " " << (al.hit ? 'T' : 'F') << " " << al.adc;
	  }
	  cout << '\n';
	  ++iLine; if (iLine >= nMaxLines) return 0;
	}
	else if (format == "high_gain_triggers") {
	  for (int ich=0; ich<int(s.high.size()); ++ich) {
	    ChipADC ah =s.high[ich], al = s.low[ich];
	    if (ah.hit || drand48() <= pedestal_suppression) { // eg. for pedestal_suppression=0.05: suppress non triggers by 20
	      cout << spill.acquisition_number() << " " << bx << " " << s.isca+1 << " "
		   << s.chip_id << " " << ich << " " << ah.adc << " " << (ah.hit ? 'T' : 'F') << '\n';
	      ++iLine; if (iLine >= nMaxLines) return 0;
	    }
	  }
	}
      }
    }
  }
  return 0;
}
