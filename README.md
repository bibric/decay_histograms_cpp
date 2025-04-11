# It's part of my homework. Skript makes diagrams of particle decay.

## Instructions

### My version of C++ compiler:
* Apple clang version 15.0.0 (clang-1500.3.9.4)  
* Target: arm64-apple-darwin23.6.0  
* Thread model: posix  

### My ROOT:
* 6.32.08 [https://root.cern](https://root.cern)

### How to compile:
```bash
g++ -O2 -o vector `root-config --cflags --libs` -lGenVector vector.cpp

