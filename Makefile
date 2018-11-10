sources = $(shell find -name "*.f90")
deps = sources2.bib data/J_ols_1000_1.png data/J_ols_1600_1.png data/J_Ridge_1000_1.png data/J_Ridge_1600_1.png data/J_LASSO_1000_1.png data/J_LASSO_400_1.png data/states.bin data/reg_nn_test_couplings.dat data/J_nn_1.png data/reg_nn_test_spins.dat

.PRECIOUS: *.dat

all:
	mkdir -p data
	make build
	make $(deps)
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

data/J_Ridge_%.dat data/J_ols_%.dat: build/linreg
	./$< <<< $*

data/J_LASSO_%.dat: programs/lasso.py
	python $< $*

build/%: build programs/%.f90
	cd build && make $*

debug: $(shell find . -name "*.f90")
	mkdir -p debug && cd debug && FC=caf cmake .. -DCMAKE_BUILD_TYPE=Debug && make

data/states.pkl:
	wget https://physics.bu.edu/~pankajm/ML-Review-Datasets/isingMC/Ising2DFM_reSample_L40_T=All.pkl -O $@

data/labels.pkl:
	wget https://physics.bu.edu/~pankajm/ML-Review-Datasets/isingMC/Ising2DFM_reSample_L40_T=All_labels.pkl -O $@

data/states_%.bin: programs/pkl2bin.py data/labels.pkl data/states.pkl
	python $< $*

data/logreg_table.dat: build/logreg programs/logreg.f90 build data/states_1.bin
	mpirun ./$<

data/reg_nn_convergence.dat: build/reg_nn_convergence
	mpirun ./$<

data/%.dat: build/% programs/%.f90 build
	./$<

data/J_nn.dat: data/reg_nn_test_couplings.dat

data/reg_nn_test_%.dat: build/reg_nn_test_%
	mpirun ./$<

clean:
	latexmk -c
	rm -rf *.run.xml *.bbl build
	find -name "*.mod" -delete
