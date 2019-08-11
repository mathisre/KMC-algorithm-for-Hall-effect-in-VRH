# KMC-algorithm-for-Hall-effect-in-VRH

This repository contains a C++ program used to simulate VRH dynamics using KMC algorithms. The code is mainly written by Andreas Glatz and Martin Kirkengen. 

run_T_simul.py can be used to start the C++ programs to simulate x and y conductivity for a range of temperature for given input parameters.

analyze_directory.py is a python script that analyzes data produced by the C++ program. It assumes the files that belong together are in the same directory, using the naming convention defined in the run_T_simul.py script and in the C++ program.

