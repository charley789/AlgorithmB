# Algorithm B for Correlation-Enhanced BP-OSD Decoding

This repository contains the implementation of **Algorithm B** for correlation-enhanced **BP-OSD decoding** of **CSS quantum codes**.

## Overview

The project focuses on a decoding framework that combines:

- **Belief Propagation (BP)** for iterative soft-information decoding
- **Order Statistic Decoding (OSD)** as a post-processing step
- **Correlation enhancement** between the two binary components in CSS quantum codes

The main implementation is written in C.

## Project Structure

The main source files are organized as follows:

- `ldpcQaun/BPOSD.c`  
  Main program and core decoding flow for Algorithm B.

- `ldpcQaun/ldpc_parm.h`  
  Header file containing the main parameter settings and compile-time configurations.

- `ldpcQaun/bp_dec/bpdec.c`  
  Implementation of the **belief propagation (BP)** decoder.

- `ldpcQaun/OSD/OSD.c`  
  Implementation of the **order statistic decoder (OSD)**.

## Flowchart

The following figure shows the overall procedure of Algorithm B.

![Algorithm B flowchart](block_diagram2.jpg)

## Notes

- The decoding process is mainly controlled from `BPOSD.c`.
- Most experiment settings and decoding parameters can be adjusted in `ldpc_parm.h`.
- The BP and OSD modules are separated into independent source files for easier modification and testing.

## Purpose

This codebase was developed for research on **correlation-enhanced BP-OSD decoding for CSS quantum codes**, with a particular focus on **Algorithm B**.
