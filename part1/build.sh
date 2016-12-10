rm rwnettest.so
f2py --f90flags='-fopenmp' -lgomp -c part1_dev.f90 -m rwnettest
