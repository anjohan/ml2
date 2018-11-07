sources = $(shell find -name "*.f90")
SHELL := /usr/bin/bash
deps = sources2.bib data/J_Ridge_1.png data/J_ols_1.png data/J_LASSO_1000_1.png data/J_LASSO_400_1.png data/states.bin data/reg_nn_test_couplings.dat data/J_nn_1.png

all:
	mkdir -p data
	make build
	make -j $(deps)
	make report.pdf

build: $(sources)
	mkdir -p build && cd build && FC=caf cmake -DCMAKE_BUILD_TYPE=Release .. && make

%.pdf: %.tex $(deps)
	latexmk -pdflua -time -g -shell-escape $*

%.pdf: %.asy
	asy -maxtile "(400,400)" -o $@ $<

sources2.bib: sources.bib
	betterbib -l $< $@

data/J_%_1.png: ./programs/matrix_to_png.py data/J_%.dat
	python $< $*

data/J_ols.dat: build/linreg
	./$<

data/J_Ridge.dat: data/J_ols.dat

data/J_LASSO_%.dat: programs/lasso.py
	python $< $*

build/%: build programs/%.f90

debug: $(shell find . -name "*.f90")
	mkdir -p debug && cd debug && FC=caf cmake .. -DCMAKE_BUILD_TYPE=Debug && make -j

data/states.pkl:
	wget https://physics.bu.edu/~pankajm/ML-Review-Datasets/isingMC/Ising2DFM_reSample_L40_T=All.pkl -O $@

data/labels.pkl:
	wget https://physics.bu.edu/~pankajm/ML-Review-Datasets/isingMC/Ising2DFM_reSample_L40_T=All_labels.pkl -O $@

data/states.bin: programs/pkl2bin.py data/labels.pkl data/states.pkl
	python $<

data/logreg_table.dat: build/logreg programs/logreg.f90 build
	mpirun ./$<

data/%.dat: build/% programs/%.f90 build
	./$<

data/J_nn.dat: data/reg_nn_test_couplings.dat

data/reg_nn_test_couplings.dat: build/reg_nn_test_couplings programs programs/reg_nn_test_couplings.f90 build
	mpirun ./$<

clean:
	latexmk -c
	rm -rf *.run.xml *.bbl build
	find -name "*.mod" -delete
