CXXFLAGS=-Wall -mtune=native -O3 -g -fPIC -lboost_python3 -fopenmp -I. -I/usr/include/python3.7m -I/usr/lib/python3.7/site-packages/numpy/core/include/ -shared -Wl,-soname,upw_loops.so 

ALL: upw_loops.so

upw_loops.so: upw_loops.cc demangle.cc
	g++ upw_loops.cc demangle.cc -o upw_loops.so ${CXXFLAGS}
