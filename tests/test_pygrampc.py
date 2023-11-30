import pygrampc as mpc
import pytest

class Problem(mpc.ProblemBase):
    def __init__(self):
        mpc.ProblemBase.__init__(self)
        self.Ng = 0
        self.NgT = 0
        self.Nh = 0
        self.NhT = 0
        self.Np = 0
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


def test_construction():
    problem = Problem()
    grampc = mpc.Grampc(problem)

def test_run():
    problem = Problem()
    grampc = mpc.Grampc(problem)

    with pytest.raises(RuntimeError):
        grampc.run()

    grampc.set_param({"dt": 1.0})
    with pytest.raises(RuntimeError):
        grampc.run()

    grampc.set_param({"Thor": 2.0})
    grampc.run()

def test_estimate_penalty():
    problem = Problem()
    grampc = mpc.Grampc(problem)
    
    grampc.estim_penmin(False)
    with pytest.raises(RuntimeError):
        grampc.estim_penmin(True)

    grampc.set_param({"dt": 1.0})
    with pytest.raises(RuntimeError):
        grampc.estim_penmin(True)

    grampc.set_param({"Thor": 2.0})
    grampc.estim_penmin(True)

if __name__ == "__main__":
    test_estimate_penalty()
