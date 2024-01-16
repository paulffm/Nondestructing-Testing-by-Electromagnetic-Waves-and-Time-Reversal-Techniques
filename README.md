# Nondestructing-Testing-by-Electromagnetic-Waves-and-Time-Reversal-Techniques

The goal of this project is to detect objects in a cube using Time Reversal Techniques inspired by the principles described in the paper "Time-reversal acoustics" by Mathias Fink [1]. Time-Reversal Acoustics is a method based on the principles of wave propagation reciprocity, which has wide-ranging applications in the field of acoustics, including medical imaging, underwater communication, and non-destructive testing.

## Overview

This repository contains the implementation of Time-Reversal Acoustics (TRA) techniques for detecting objects within a resonant box. The project is divided into multiple modules, each focusing on specific aspects of electromagnetic wave simulations and TRA-based object detection.

## Modules

### Module 1: Eigenvalue FIT Solver

In this module, we provide analytical and numerical solutions to the eigenvalue problem within a box resonator. Additionally, we conduct a convergence study to assess the accuracy and efficiency of the solver.

### Module 2: Frequency Domain FIT Solver

Module 2 extends the capabilities of the FIT solver to the frequency domain, considering spatial distributions. A convergence study is also included to evaluate the solver's performance.

### Module 3: Time Domain FIT Solver with Chirp

This module introduces a time domain FIT solver with electrical stimulation using chirp signals.

### Module 4: Time Reversal Techniques

Module 4 explores Time Reversal Techniques to restore and focus signals, a fundamental concept in Time-Reversal Acoustics.

### Module 5: Cube Detection Inside a Box Resonator

In this final module, we apply Time-Reversal Techniques to detect the presence and position of cubes within the box resonator. The objective is to locate cubes using TRA-based methods.

[1]: https://iopscience.iop.org/article/10.1088/1742-6596/118/1/012001/pdf
