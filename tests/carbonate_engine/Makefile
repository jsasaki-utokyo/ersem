FC = gfortran
FFLAGS ?= -O2 -fcheck=all
SRC_DIR = ../../src

test_driver: fabm_types_mock.F90 $(SRC_DIR)/carbonate_engine.F90 test_driver.f90
	$(FC) $(FFLAGS) -o $@ $^

test: test_driver
	./test_driver reference_cases.csv

clean:
	rm -f test_driver *.o *.mod

.PHONY: test clean
