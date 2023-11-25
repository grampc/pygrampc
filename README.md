# GRAMPC Python Interface

## Description

A Python interface for the GRAMPC MPC solver.

## Installation
Clone this repository and install the package via pip. You need a suitable C++ compiler installed listed at the [pybind11](https://github.com/pybind/pybind11) project.
```
pip install "path to interface"
```

## Problem Description
The problem description can either be stated in C++ or Python.

### Python Problem Description
In Python you can inherit from the C++ `ProblemBinding` class. 
For an example consider the Python template or the implemented problems inside the tests folder. 
You only have to define the functions you need, the rest don't need to be stated explicitly.

The input variables are all numpy arrays, so they natively can be used inside numpy functions. 

To set the output data of your function do write to the preallocated memory inside the out array with
```
out[0] = ...
out[:] = ...
```

statements like
```
out = out * 2
```
do not work because Python creates a new variable `out` which doesn't points to the allocated memory from GRAMPC

### C++ Problem Description
If extra speed is desired, the Python interface can be used with a compiled problem description. 
Inside the template folder there is a C++ template with a CMakeLists.txt file for compilation. 
As in the Python interface, you only have to override the functions you need.

Using the C++ problem description can result in a 100 times speedup per time step.
