PYTHON := /usr/bin/python
CYTHON := /usr/bin/cython


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
