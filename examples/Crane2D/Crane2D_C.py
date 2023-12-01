import numpy as np
import matplotlib.pyplot as plt
import pygrampc as mpc
from scipy.integrate import solve_ivp
# import the C++ problem description. Might be a different path depending on the used compiler. Look for .pyd or .so files
from build.Release.crane_problem import Crane2D


if __name__ == "__main__":
    Tsim = 12.5
    plotSteps = 150
    path = "Crane2D.json"

    Q = np.array([1.0, 2.0, 2.0, 1.0, 1.0, 4.0])
    R = np.array([0.05, 0.05])
    problem = Crane2D(Q, R, 0.2, 1.25, 0.3)

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

        # plots of the grampc predictions
        if i % plotSteps == 0:
            grampc.plot()
            vec.plot()

    vec.plot()
    plt.show()
