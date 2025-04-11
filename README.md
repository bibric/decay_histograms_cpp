#It's part of my homework. Skript makes diagrams of particle decay.
## Instructions
###My version of C++ compilator:
 * Apple clang version 15.0.0 (clang-1500.3.9.4)
 * Target: arm64-apple-darwin23.6.0
 * Thread model: posix
###My ROOT:
 * 6.32.08  https://root.cern
###How to compile:
 * g++ -O2 -o vector `root-config --cflags --libs` -lGenVector vector.cpp
###Output files will be 1 histogram in two files:
 * invmass.pdf
 * invmass.png 
