from pygrampc import ProblemBase, Grampc
import pytest
import numpy as np


class Problem(ProblemBase):
    def __init__(self):
        ProblemBase.__init__(self)
        self.Ng = 0
        self.NgT = 0
        self.Nh = 0
        self.NhT = 0
        self.Np = 1
        self.Nu = 1
        self.Nx = 1

    def ffct(self, out, t, x, u, p):
        out[0] = -x[0] + u[0]

    def dfdx_vec(self, out, t, x, vec, u, p):
        out[0] = -vec[0]
    
    def dfdu_vec(self, out, t, x, vec, u, p):
        out[0] = vec[0]

    def lfct(self, out, t, x, u, p, xdes, udes):
        out[0] = x[0]**2 + u[0]**2

    def dldx(self, out, t, x, u, p, xdes, udes):
        out[0] = 2*x[0]

    def dldu(self, out, t, x, u, p, xdes, udes):
        out[0] = 2*u[0]


def test1_run():
    problem = Problem()
    grampc = Grampc(problem)

    with pytest.raises(RuntimeError):
        grampc.run()

    grampc.set_param({"dt": 1.0})
    with pytest.raises(RuntimeError):
        grampc.run()

    grampc.set_param({"Thor": 2.0})
    grampc.run()


def test2_estimate_penalty():
    problem = Problem()
    grampc = Grampc(problem)
    
    grampc.estim_penmin(False)
    with pytest.raises(RuntimeError):
        grampc.estim_penmin(True)

    grampc.set_param({"dt": 1.0})
    with pytest.raises(RuntimeError):
        grampc.estim_penmin(True)

    grampc.set_param({"Thor": 2.0})
    grampc.estim_penmin(True)


def test3_set_get_parameters():
    problem = Problem()
    grampc = Grampc(problem)

    parameters = {
        "x0": [0.1],
        "xdes": [1.0],
        "u0": [0.1],
        "udes": [0.0],
        "umax": [1.0],
        "umin": [-1.0],
        "p0": [3.0],
        "pmax": [4.0],
        "pmin": [2.0],
        "Thor": 3600.0,
        "Tmax": 100.0,
        "Tmin": 1.0,
        "dt": 30.0,
        "t0": 2.0
    }

    grampc.set_param(parameters)
    for key, value in parameters.items():
        assert getattr(grampc.param, key) == parameters[key]

    with pytest.raises(KeyError):
        grampc.set_param({"MyKey": 10})
    
    with pytest.raises(ValueError, match="x0: dimension does not match problem description. Expected 1, got 2"):
        grampc.set_param({"x0": [1.0, 2.0]})


def test4_set_rws():
    problem = Problem()
    grampc = Grampc(problem)

    parameters = {
        "x0": [0.1],
        "xdes": [1.0],
        "u0": [0.1],
        "udes": [0.0],
        "umax": [1.0],
        "umin": [-1.0],
        "p0": [3.0],
        "pmax": [4.0],
        "pmin": [2.0],
        "Thor": 3600.0,
        "Tmax": 100.0,
        "Tmin": 1.0,
        "dt": 30.0,
        "t0": 2.0
    }

    grampc.set_param(parameters)
    with pytest.raises(ValueError, match="Wrong dimensions detected. Expected \\(1, 30\\), got \\(1, 4\\)"):
        grampc.set_rws_u(np.zeros((1, 4)))

    grampc.set_rws_u(np.zeros((1, 30)))


if __name__ == "__main__":
    test4_set_rws()
