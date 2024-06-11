# EC-Solver:

## Overview:

This code is built upon the openly available code [FBS-Solver](https://github.com/DMGW-Goethe/FBS-Solver). See there for more details about the code itself.

The purpose of EC-Solver is to compute stable configurations of neutron stars in Einstein-Cartan theory.
For more details about the theory, we refer to our publication (see below).

## Publication

This code has been used in our [publication [arXiv:2406.05851]](https://arxiv.org/abs/2406.05851). All raw data used in this publication and all figures can be generated using the code provided in this repository.

## Compiling

To compile, call

    make fbs

This will compile the c++ code without the Cython components that are included in FBS-Solver.

In case dependencies are needed, run:

    sudo apt install libboost-dev
    pip3 install numpy matplotlib cython pandas

## Running The Code

Before running, you must create a folder named `output` in the top level folder of this repository. To generate all raw data, run `./main.out`. The code execution will take 15-25 minutes on a laptop, depending on your machine.
The data can then be plotted using the python plotting scripts in the folder `EC_plotting_scripts`. To produce the figures from our publication, simply run the individual plotting scripts using `python3 figure1.py` (and analogously for figure2-4.py).

If you want to run cases with your own parameters, see the `main.cpp` for some examples of how to integrate a single star, or create an MR-curve.
