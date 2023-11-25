import numpy as np
import matplotlib.pyplot as plt
import pygrampc as mpc
from scipy.integrate import solve_ivp


class Crane2D(mpc.ProblemBase):
    def __init__(self, Q: np.array, R: np.array, ScaleConstraint, MaxConstraintHeight, MaxAngularDeflection):
        mpc.ProblemBase.__init__(self)
        self.Nx = 6
        self.Nu = 2
        self.Np = 0
        self.Ng = 0
        self.Nh = 3
        self.NgT = 0
        self.NhT = 0

        self.Q = Q
        self.R = R
        self.ScaleConstraint = ScaleConstraint
        self.MaxConstraintHeight = MaxConstraintHeight
        self.MaxAngularDeflection = MaxAngularDeflection

    def ffct(self, out, t, x, u, p):
        out[0] = x[1]
        out[1] = u[0]
        out[2] = x[3]
        out[3] = u[1]
        out[4] = x[5]
        out[5] = -((9.81 * np.sin(x[4]) + np.cos(x[4]) * u[0] + 2 * x[3] * x[5]) / x[2])

    def dfdx_vec(self, out, t, x, vec, u, p):
        sinX = np.sin(x[4])
        cosX = np.cos(x[4])
        g = 9.81

        out[0] = 0
        out[1] = vec[0]
        out[2] = (g * sinX + cosX * u[0] + 2 * x[3] * x[5]) * vec[5] / x[2]**2
        out[3] = vec[2] - (2 * x[5] * vec[5]) / x[2]
        out[4] = -((g * cosX - sinX * u[0]) * vec[5] / x[2])
        out[5] = vec[4] - (2 * x[3] * vec[5]) / x[2]

    def dfdu_vec(self, out, t, x, vec, u, p):
        out[0] = vec[1] - (np.cos(x[4]) * vec[5]) / x[2]
        out[1] = vec[3]

    def lfct(self, out, t, x, u, p, xdes, udes):
        out[0] = np.dot(self.Q, np.power(x - xdes, 2)) + np.dot(self.R, np.power(u - udes, 2))

    def dldx(self, out, t, x, u, p, xdes, udes):
        out[:] = 2 * self.Q * (x - xdes)

    def dldu(self, out, t, x, u, p, xdes, udes):
        out[:] = 2 * self.R * (u - udes)

    def hfct(self, out, t, x, u, p):
        Position = x[0] + np.sin(x[4]) * x[2]

        out[0] = np.cos(x[4]) * x[2] - self.ScaleConstraint * np.power(Position, 2) - self.MaxConstraintHeight
        out[1] = x[5] - self.MaxAngularDeflection
        out[2] = -x[5] - self.MaxAngularDeflection

    def dhdx_vec(self, out, t, x, u, p, vec):
        tmp = self.ScaleConstraint * (x[0] + np.sin(x[4]) * x[2])

        out[0] = tmp * vec[0]
        out[1] = 0
        out[2] = (np.sin(x[4]) * tmp + np.cos(x[4])) * vec[0]
        out[3] = 0
        out[4] = (np.cos(x[4]) * x[2] * tmp - np.sin(x[4]) * x[2]) * vec[0]
        out[5] = 0 + vec[1] - vec[2]

    def dhdu_vec(self, out, t, x, u, p, vec):
        out[0] = 0


if __name__ == "__main__":
    Tsim = 12.5
    plotSteps = 150
    path = "Crane2D.json"

    Q = np.array([1.0, 2.0, 2.0, 1.0, 1.0, 4.0])
    R = np.array([0.05, 0.05])
    Param = [0.2, 1.25, 0.3]
    problem = Crane2D(Q, R, Param)

    # initialize grampc object
    grampc = mpc.Grampc(problem, path)

    # estimate penaltyMin and set option
    grampc.estim_penmin(True)
    grampc.print_opts()
    grampc.print_params()

    # construct solution structure
    vec = mpc.GrampcResults(grampc, Tsim, plot_results=True, plot_statistics=True)

    dt = grampc.param.dt

    for i, t in enumerate(vec.t):
        vec.CPUtime[i] = grampc.run()
        vec.update(grampc, i)

        if i + 1 > len(vec.t):
            break

        # simulate the system
        sol = solve_ivp(grampc.ffct, [t, t+dt], grampc.param.x0,
                        args=(grampc.sol.unext, grampc.sol.pnext))

        # set current time and current state
        grampc.set_param({"x0": sol.y[:, -1],
                          "t0": t + dt})

        # evaluate inequality constraints
        vec.constr[i, :] = grampc.hfct(vec.t[i], vec.x[i, :], vec.u[i, :], vec.p[i, :])

        # plots of the grampc predictions
        if i % plotSteps == 0:
            grampc.plot()
            vec.plot()

    vec.plot()
    plt.show()
