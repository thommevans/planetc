

INSTALL_OPT = '--prefix=.'


.PHONY: install build cython_files clean

install:
	python setup.py install ${INSTALL_OPT}

build:
	python setup.py build

cython_files: 
	rm -f planetc/keporb.c
	rm -f planetc/ma02.c
	cython planetc/keporb.pyx
	cython planetc/ma02.pyx

clean:
	rm -rf build/
