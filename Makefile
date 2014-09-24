PYTHON := /usr/bin/python3
CYTHON := /usr/bin/cython3


INSTALL_OPT = '--prefix=.'


.PHONY: install build cython_files clean

install:
	${PYTHON} setup.py install ${INSTALL_OPT}

build:
	${PYTHON} setup.py build

cython_files: 
	rm -f planetc/keporb.c
	rm -f planetc/ma02.c
	${CYTHON} planetc/keporb.pyx
	${CYTHON} planetc/ma02.pyx

clean:
	rm -rf build/
