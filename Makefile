X-py.bin Z-py.bin L-py.bin: data.feather formula.txt
	python3 build_matrices.py --data data.feather --formula formula.txt \
		-X X-py.bin -Z Z-py.bin -L L-py.bin

data.feather: columns.txt
	python3 generate_data.py --columns $^ --data $@

clean:
	rm -rf $(wildcard *.feather *.bin)
