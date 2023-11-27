# PyGRAMPC

PyGRAMPC is a Python interface for the GRAMPC solver build by pybind11. It features at least the 
same functionality as the Matlab interface. For questions regarding GRAMPC refer to the [grampc](https://github.com/grampc/grampc)
repository.

## Features 

 - Independent installation of GRAMPC via pip.
 - Writing the problem description in Python or C++.
 - Debug your Python problem description directly in Python
 - C allocated arrays are exposed via numpy arrays.

## Installation
Clone this repository and install the package via pip. 
You need a suitable C++ compiler installed listed at the [pybind11](https://github.com/pybind/pybind11) project.
```
pip install "path to interface"
```

## Usage
Please refer to the two examples inside the `examples` folder.

## Problem Description
The problem description can either be stated in C++ or Python.

### Python Problem Description
In Python you can inherit from the C++ `ProblemBase` class supplied by the pygrampc package. 
This class redirects the C function calls to Python, so native Python code can be executed. 
Please note that the Python problem description is intended for rapid prototyping. 
If the full speed of GRAMPC is desired, one must write the problem description in C++

#### Initialization
The attributes Nx, Nu, Np, Ng, Nh, NgT and NhT must be initalized during `__init__()` to circumvent undefined behaviour. 
Setting these values during `__init__()` is equivalent to the `void ocp_dim(...)` function from the C interface.

An example would be
```python
class MyProblem(ProblemBase):
    def __init__(self):
        ProblemBase.__init__(self, *args)
        self.Nx = 2
        self.Nu = 1
        self.Np = 0
        self.Ng = 0
        self.Nh = 0
        self.NgT = 2
        self.NhT = 0
        ...
```
When writing your problem description in Python, only `__init__()`, `ffct()` and `dfdx()` are mandatory to implement. 
If you are unsure about the function arguments, please refer to `MyProblem.py`.

#### Usage
All function arguments are numpy arrays. The output of the functions need to be written at the allocated memory of `out`.
To set the output data of your function do write to the preallocated memory inside the out array with
```python
out[0] = ...
out[:] = ...
```

statements like
```python
out = out * 2
```
do not work because Python creates a new variable `out` which doesn't points to the allocated memory from GRAMPC

Because the arguments are numpy arrays, one should use numpy methods for calculation. 
Here is an example for computing the integral cost term:
```python
def lfct(self, out, t, x, u, p, xdes, udes):
    out[0] = np.dot(self.Q, np.power(x - xdes, 2)) + np.dot(self.R, np.power(u - udes, 2))
```

### C++ Problem Description
If extra speed is desired, the Python interface can be used with a compiled problem description. 
Inside the template folder there is a C++ template with a CMakeLists.txt file for compilation. 
As in the Python interface, you only have to override the functions you need.

Using the C++ problem description can result in a 100 times speedup per time step.

For an example, please refer to the `Crane2D` example. 
There exists a Python and C++ problem description, which showcases the achieved speedup.


