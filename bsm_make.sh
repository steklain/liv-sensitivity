#!/bin/bash
myname=$@
gcc -Wall -c bsm.c -std=c99
g++ -Wall -c "$myname".cc
g++ -Wall "$myname".o bsm.o -lglobes -lgsl -lgslcblas -O3 -o "$myname"
./"$myname"
rm *.o