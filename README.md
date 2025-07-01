# TRACE - Topological Ring and Additive-Coordinated Cage Explorer
## 1. Introduction

TRACE (Topological Ring and Additive-Coordinated Cage Explorer) is a high-performance algorithm written in `C++`, developed to identify and classify hydrate cage structures from molecular dynamics (MD) simulations of clathrate hydrates. Designed to analyze GROMACS-format trajectories (`.gro`), TRACE provides comprehensive detection of structural motifs including rings, cups, incomplete cages, and complete cages—particularly those formed through additive–water interactions.

TRACE achieves both high accuracy and computational efficiency across a wide range of system sizes and topologies, from hundreds to millions of particles. Its key features include:

- Detection of complete cages and near-complete (incomplete) cages.
- Ring enumeration up to 12-membered rings (user-adjustable).
- Custom hydrogen bond definitions to identify additive-coordinated cages (currently supports one additive species, user-extensible).
- Automated identification of cage clusters.
- Full tracking of cage-forming water and additive molecules, enabling post-analysis of cage lifetime and additive residence time.
- Crystallinity estimation of the entire system.
- Guest molecule occupancy analysis for each cage.

This tool was developed as part of an academic study on hydrate formation dynamics. For more information, please refer to the upcoming publication [reference pending].

## 2. Usage

### 2.1 Clone the repository

Clone the repository to your desired directory (replace `/path/to/TRACE` with your chosen path):

`git clone https://github.com/qwe88866/TRACE.git /path/to/TRACE`

### 2.2 Compile the program

Compile TRACE using `g++` with OpenMP support:
`g++ -fopenmp /path/to/TRACE/TRACE.cpp -o /path/to/TRACE/TRACE.exe`
- A C++ compiler supporting at least C++11 standard (e.g., g++ 4.8 or later).
- OpenMP support enabled (g++ 4.8 or later recommended).

### 2.3 Run TRACE

`cd /path/to/TRACE/example`

TRACE automatically adjusts computation based on the input files provided. Some example usages:

- Compute only water cages:

`/path/to/TRACE/TRACE.exe -w example_H2O.gro`

- Compute cages for water and guest molecules:

`/path/to/TRACE/TRACE.exe -w example_H2O.gro -g example_guest.gro`

- Compute cages for water, urea additives, and guests:

`/path/to/TRACE/TRACE.exe -w example_H2O.gro -g example_guest.gro -a example_urea.gro -h hbond_urea.txt`

**Note:** See [Section 3: File Format](#3-file-format) for details on file structure.

**Note:** We generate a total of 9 output files for analysis.  
All output files are saved by default in the current working directory.  
The output file names are fixed (non-randomized), so please be careful to avoid overwriting existing files.

For detailed information, see Section 4. Output Files and Analysis.

### 2.4 Visualization with VMD

To visualize the results, you need to install VMD (Visual Molecular Dynamics).

Once installed, launch VMD from the terminal:

`vmd`

In the VMD Tk Console, run the following command to load the visualization script:

`source /path/to/TRACE/visualize.tcl`

## 3. File Format

Water, guest molecules, and additives must be provided as separate `.gro` files.  
You can use GROMACS' `trjconv` and `make_ndx` tools to extract and create these separate index groups from your full system trajectory.

For details on the `.gro` file format used by TRACE, please refer to the official GROMACS documentation:  
[GROMACS 5.0.4 `.gro` file format](https://manual.gromacs.org/archive/5.0.4/online/gro.html)

Because the program does **not** require users to input the number of atoms per molecule (e.g., for H2O, guest, additive),  
it relies on **sequential and consistent (`resid`)** IDs  to identify molecule boundaries. This is the default behavior for GROMACS `.gro` output.

For example:
| Correct format        | Wrong format          |
|-----------------------|-----------------------|
|    1H2O   OICE  ...   |    1H2O   OICE  ...   |
|    1H2O   HICE  ...   |    1H2O   HICE  ...   |
|    1H2O   HICE  ...   |    1H2O   HICE  ...   |
|    2H2O   OICE  ...   |    1H2O   OICE  ...   |
|    2H2O   HICE  ...   |    1H2O   HICE  ...   |
|    2H2O   HICE  ...   |    1H2O   HICE  ...   |
|    3H2O   OICE  ...   |    1H2O   OICE  ...   |
|    ...                |    ...                |

### 3.1 H₂O file 

- The first atom must be **O** (oxygen).  
- If present, **H1** and **H2** atoms must immediately follow as the second and third atoms, respectively.  
- Supported water models:  
  - **O**       : 1-point water  
  - **OHH**     : 3-point water  
  - **OHHM**    : 4-point water  
  - **OHHLL**   : 5-point water
    
For example:
| Correct format       | Wrong format         | Wrong format         |
|----------------------|----------------------|----------------------|
|    1H2O   OICE ...   |    1H2O   HICE ...   |   1H2O   HICE ...    |
|    1H2O   HICE ...   |    1H2O   OICE ...   |   1H2O   HICE ...    |
|    1H2O   HICE ...   |    1H2O   HICE ...   |   1H2O   OICE ...    |

> *Note: Dummy atoms (if any) will be ignored.*
> *If hydrogen atoms (H1, H2) are present, hydrogen bonds and bond angles will be considered; otherwise, these features are disabled.*

### 3.2 Guest file 

- The guest molecule uses the **first atom as its center point**.  
- To include multiple guest molecule types (e.g., CO₂, CH₄), extract a single representative atom per molecule (such as the carbon atom for CO₂ and CH₄) to serve as the molecule’s center point.
- Use `make_ndx` to create index groups for each guest type, then use `trjconv` to extract these single-point representations and merge them into a single `.gro` file for processing by TRACE.

### 3.3 Additive Hydrogen Bond Definition File

TRACE allows users to define custom hydrogen bonds between additive and water molecules by providing an external file via the `-h` option.

Each hydrogen bond entry must start with `!` and follow the format:
`! donor_type acceptor_type donor_atom_index acceptor_atom_index donor_H_atom_index theta_cutoff`
Where:
- `donor_type`, `acceptor_type`: either `w` (water) or `a` (additive)  
- `donor_atom_index`, `acceptor_atom_index`, and `donor_H_atom_index`: indices according to the `.gro` file  
- `theta_cutoff`: angular cutoff in degrees  
- The order of entries does not affect the outcome

Atom indices **must match the order in the `.gro` file**. For example, water must be ordered as `O H1 H2 ...`, and the additive atoms should follow their appearance in the file.  

For H₂O.gro example:
|---------------------------|
|    1H2O  OICE1    1 ...   | 
|    1H2O  HICE2    2 ...   | 
|    1H2O  HICE3    3 ...   |
|    2H2O  OICE1    4 ...   |
|    ...                    |

For urea.gro example:
|---------------------------|
|    1Urea    C1    1 ...   |  
|    1Urea    O2    2 ...   | 
|    1Urea    N3    3 ...   |
|    1Urea    N4    4 ...   |  
|    1Urea    H5    5 ...   | 
|    1Urea    H6    6 ...   |
|    1Urea    H7    7 ...   |  
|    1Urea    H8    8 ...   | 
|    2Urea    C1    9 ...   |
|    ...                    |

#### Sample content of `hbond_urea.txt`
<pre>
# water-additive interaction
! w a 0 1 1 35 #OICE1(0)-O2(1)-HICE2(1)
! w a 0 1 2 35 #OICE1(0)-O2(1)-HICE3(2)
! a w 2 0 4 35 #N3(2)-OICE1(0)-H5(4)
! a w 2 0 5 35 #N3(2)-OICE1(0)-H6(5)
! a w 3 0 6 35 #N4(3)-OICE1(0)-H7(6)
! a w 3 0 7 35 #N4(3)-OICE1(0)-H8(7)

# additive-additive interaction
! a a 2 1 4 35 #N3(2)-O2(1)-H5(4)
! a a 2 1 5 35 #N3(2)-O2(1)-H6(5)
! a a 3 1 6 35 #N4(3)-O2(1)-H7(6)
! a a 3 1 7 35 #N4(3)-O2(1)-H8(7)
</pre>

> 💡 **Note:** water–water hydrogen bonds are built-in and do not need to be listed.
You can refer to the provided [`hbond_urea.txt`](./example/hbond_urea.txt) file for a complete working example.
