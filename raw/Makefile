all : raw split_raw

raw : raw_main.C raw.C raw.hh
	g++ -I. -O3 raw_main.C raw.C -o raw
split_raw : split_raw.C raw.C raw.hh
	g++ -I. -O3 split_raw.C raw.C -o split_raw
