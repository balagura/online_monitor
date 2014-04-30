// Author V. Balagura, balagura@cern.ch (19.11.2012)

#include <raw.hh>
#include <iostream>
#include <cstdlib>
#include <string>
#include <climits> // for INT_MAX

int main(int argc, char** argv) {
  int nMaxLines = INT_MAX;
  double pedestal_suppression = 0.05;
  string format = "";
  if (argc < 2 || argc > 5) {
    cout << "Usage: " << argv[0] << " <file_name> [max_N_lines]\n"
      "or " << argv[0] << " <file_name> max_N_lines high_gain_triggers [pedestal_suppression]\n"
      "If max_N_lines is absent or 0, all file is processed\n"
      "If pedestal_suppression is given, all not triggered data\n"
      "is postscaled by this factor (between 0 and 1, default 0.05)\n";
    return 0;
  }
  else {
    if (argc >= 3) {
      int n = atoi(argv[2]);
      if (n != 0) nMaxLines = n;
    }
    if (argc >= 4) format = argv[3];
    if (argc >= 5) pedestal_suppression = atof(argv[4]);
  }
  ReadSpill spill(argv[1]);
  int iLine=0;

  int prev_acq = -1; // acquisition_number has 16 bits; when it goes down (eg. from 65535 to 0): add 65536
  int acq;

  while (spill.next()) {
    if      (prev_acq == -1)                         acq  =  spill.acquisition_number();
    else if (spill.acquisition_number() >= prev_acq) acq += (spill.acquisition_number() - prev_acq);
    else                                             acq += (spill.acquisition_number() - prev_acq + 65536);
    prev_acq = spill.acquisition_number();

    const map<unsigned short int, vector<ChipSCA> >& sca = spill.sca();
    for (map<unsigned short int, vector<ChipSCA> >::const_iterator isca=sca.begin(); isca!=sca.end(); ++isca) {
      int bx  =  isca->first;
      for (int i=0; i<int(isca->second.size()); ++i) {
	const ChipSCA& s = isca->second[i];
	if (format == "") {
	  cout << acq << " " << bx << " " << s.chip_id << " " << s.isca+1;
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
	      cout << acq << " " << bx << " " << s.isca+1 << " "
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
