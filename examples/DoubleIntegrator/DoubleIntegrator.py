from pygrampc import ProblemDescription, Grampc, GrampcResults
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


class DoubleIntegrator(ProblemDescription):
    def __init__(self):
        ProblemDescription.__init__(self)
        self.Nx = 2
        self.Nu = 1
        self.Np = 0
        self.Ng = 0
        self.Nh = 0
        self.NgT = 2
        self.NhT = 0

        self.CostIntegral = 0.1
        self.CostTerminal = 1.0

    def ffct(self, out, t, x, u, p):
        out[0] = x[1]
        out[1] = u[0]

    def dfdx_vec(self, out, t, x, vec, u, p):
        out[0] = 0
        out[1] = vec[0]

    def dfdu_vec(self, out, t, x, vec, u, p):
        out[0] = vec[1]

    def lfct(self, out, t, x, u, p, xdes, udes):
        out[0] = self.CostIntegral * (u[0] - udes[0]) ** 2

    def dldu(self, out, t, x, u, p, xdes, udes):
        out[0] = 2 * self.CostIntegral * (u[0] - udes[0])

    def Vfct(self, out, T, x, p, xdes):
        out[0] = self.CostTerminal * T

    def dVdT(self, out, T, x, p, xdes):
        out[0] = self.CostTerminal

    def gTfct(self, out, T, x, p):
        out[0] = x[0]
        out[1] = x[1]

    def dgTdx_vec(self, out, T, x, p, vec):
        out[0] = vec[0]
        out[1] = vec[1]


if __name__ == "__main__":
    Tsim = 8
    plotSteps = 150
    options = "DoubleIntegrator.json"

    # initialize problem and GRAMPC
    Problem = DoubleIntegrator()
    grampc = Grampc(Problem, options, plot_prediction=False)

    # construct solution structure
    vec = GrampcResults(grampc, Tsim, plot_results=True, plot_statistics=True)

    dt = grampc.param.dt

    for i, t in enumerate(vec.t):
        vec.CPUtime[i] = grampc.run()
        vec.update(grampc, i)

        if i + 1 > len(vec.t) or vec.t[i + 1] > Tsim:
            break

        # simulate system
        sol = solve_ivp(grampc.ffct, [t, t + dt], grampc.param.x0,
                        args=(grampc.sol.unext, grampc.sol.pnext))

        # set current time and state
        grampc.set_param({"x0": sol.y[:, -1],
                         "t0": t + dt})

        if grampc.sol.Tnext <= grampc.param.Tmin + grampc.param.dt and bool(grampc.opt.OptimTime):
            grampc.set_opt({"OptimTime": "off"})
            Tsim = vec.t[i + 1]

        # plots of the grampc predictions
        if i % plotSteps == 0:
            grampc.plot()
            vec.plot()

    vec.plot()
    plt.show()
