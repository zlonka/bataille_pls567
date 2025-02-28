bn:	bataille.o.cpp
	gcc -DRANGEMENT=NAT -O3 -march=native -o bn bataille.o.cpp && echo "bn ready"
bo:	bataille.o.cpp
	gcc -DRANGEMENT=OPT -O3 -march=native -o bo bataille.o.cpp && echo "bo ready"
all:	bn bo
