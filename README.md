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
