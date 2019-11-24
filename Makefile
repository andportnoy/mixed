check: check_design

# we can get -0 vs +0, so binary comparison won't work here
check_design: X-py.bin Z-py.bin L-py.bin X-r.bin Z-r.bin L-r.bin
	python3 cmp.py -a X-py.bin -b X-r.bin
	python3 cmp.py -a Z-py.bin -b Z-r.bin
	python3 cmp.py -a L-py.bin -b L-r.bin

X-r.bin Z-r.bin L-r.bin: data.feather formula.txt
	Rscript build_matrices.r --data data.feather --formula formula.txt \
		--X X-r.bin --Z Z-r.bin --L L-r.bin

X-py.bin Z-py.bin L-py.bin: data.feather formula.txt
	python3 build_matrices.py --data data.feather --formula formula.txt \
		--X X-py.bin --Z Z-py.bin --L L-py.bin

data.feather: columns.txt
	python3 generate_data.py --columns $^ --data $@

clean:
	rm -rf $(wildcard *.feather *.bin)
