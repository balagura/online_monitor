all : raw

raw : raw_main.C raw.C raw.hh
	g++ -I. -O3 raw_main.C raw.C -o raw
print_spill : print_spill.C raw_new.C raw_new.hh
	g++ -I. -g print_spill.C raw_new.C -o print_spill
raw_event : raw_event.C raw.C raw.hh
	g++ -I. -O3 raw_event.C raw.C -o raw_event
