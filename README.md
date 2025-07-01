# TRACE - Topological Ring and Additive-Coordinated Cage Explorer
## Introduction

TRACE (Topological Ring and Additive-Coordinated Cage Explorer) is a high-performance algorithm developed to identify and classify hydrate cage structures from molecular dynamics (MD) simulations of clathrate hydrates. Designed to analyze GROMACS-format trajectories (`.gro`), TRACE provides comprehensive detection of structural motifs including rings, cups, incomplete cages, and complete cages—particularly those formed through additive–water interactions.

TRACE achieves both high accuracy and computational efficiency across a wide range of topologies. Its key features include:

- Detection of complete cages and near-complete (incomplete) cages.
- Ring enumeration up to 12-membered rings (user-adjustable).
- Custom hydrogen bond definitions to identify additive-coordinated cages (currently supports one additive species, user-extensible).
- Automated identification of cage clusters.
- Full tracking of cage-forming water and additive molecules, enabling post-analysis of cage lifetime and additive residence time.
- Crystallinity estimation of the entire system.
- Guest molecule occupancy analysis for each cage.

This tool was developed as part of an academic study on hydrate formation dynamics. For more information, please refer to the upcoming publication [reference pending].

## Usage

### 1. Clone the repository

Clone the repository to your desired directory (replace `/path/to/TRACE` with your chosen path):

`git clone https://github.com/qwe88866/TRACE.git /path/to/TRACE`

### 2. Compile the program

Compile TRACE using `g++` with OpenMP support:
`g++ -fopenmp /path/to/TRACE/TRACE.cpp -o /path/to/TRACE/TRACE.exe`
- A C++ compiler supporting at least C++11 standard (e.g., g++ 4.8 or later).
- OpenMP support enabled (g++ 4.8 or later recommended).

### 3. Run TRACE

`cd /path/to/TRACE/example`

TRACE automatically adjusts computation based on the input files provided. Some example usages:

- Compute only water cages:

`/path/to/TRACE/TRACE.exe -w example_H2O.gro`

- Compute cages for water and guest molecules:

`/path/to/TRACE/TRACE.exe -w example_H2O.gro -g example_guest.gro`

- Compute cages for water, urea additives, and guests:

`/path/to/TRACE/TRACE.exe -w example_H2O.gro -g example_guest.gro -a example_urea.gro -h hbond_urea.txt`

### Visualization with VMD

To visualize the results, you need to install VMD (Visual Molecular Dynamics).

After installing and launching VMD, load the visualization script by running the following command in the VMD console (Tcl prompt):

'source /path/to/TRACE/visualize.tcl'

Note: When calculating additive-coordinated cages, a valid hydrogen bond definition file conforming to our specifications must be provided. The format and details of this file will be discussed in the "Command Line Options" section.

Note: Water, guest molecules, and additives must be provided as separate `.gro` files. You can use GROMACS' `gmx_mpi make_ndx` tool to create these separate index groups from your full system trajectory.

Note: We generate a total of 9 output files for analysis.
For detailed information, see Section 5: Output Files and Analysis.
