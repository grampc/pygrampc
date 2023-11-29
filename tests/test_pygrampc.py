import pygrampc as mpc


def test_construction():
    class TestProblem(mpc.ProblemBase):
        def __init__(self):
            mpc.ProblemBase.__init__(self)
            self.Ng = 0
            self.NgT = 0
            self.Nh = 0
            self.NhT = 0
            self.Np = 0
            self.Nu = 0
            self.Nx = 1

        def ffct(self, out, t, x, u, p):
            out[0] = -x[0]

        def dfdx_vec(self, out, t, x, vec, u, p):
            out[0] = -vec[0]

    problem = TestProblem()
    grampc = mpc.Grampc(problem)

    grampc.set_param({
        "dt": 1.0,
        "Thor": 2.0
        })

    grampc.run()
