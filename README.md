# CKKS-Reference

Codes for the empirical evidence that composite scaling factors work. 

This script originated from the fork of ISO standardization codes (2023.04.19).



## How to compile?

You can build and compile the script using Microsoft Visual Studio. You can open the solution file (CKKS_RNS.sln) with Microsoft Visual Studio and then build and compile the projects. 

### Projects

There are four projects in the solution file CKKS_RNS.sln. 

- "CKKS_all combined" project
  
  This project is the original ISO standardization project (forked on 2023.04.19). 

- "64-bit Prime" project
  
  This project is for the error benchmark of CKKS with scaling factor = $2^{58}$  and 64-bit RNS primes. See 'main-58bit.cpp' for the details.

- "32-bit Prime" project
  
  This project is for the error benchmark of CKKS with scaling factor = $2^{58}$ and 32-bit RNS primes. See 'main-32bit.cpp' for the details. For each rescaling and multiplication, it consumes two RNS primes. 

- "32-bit Precision" project
  
  This project is for the error benchmark of CKKS with scaling factor = $2^{28}, 2^{29}, 2^{30}, 2^{31}$, and 32-bit RNS primes. See 'main-32pre.cpp' for the details. 

### Develop Environment

- Visual Studio @17.5.3

- Windows 10


