import numpy as np
import matplotlib.pyplot as plt
from pygrampc import ProblemBase, Grampc, GrampcResults
from scipy.integrate import solve_ivp

# import the C++ problem description. Might be a different path depending on the used compiler. Look for .pyd or .so files
from crane_problem import Crane2D


if __name__ == "__main__":
    Tsim = 12.5
    plotSteps = 150
    options = "Crane2D.json"

    Q = np.array([1.0, 2.0, 2.0, 1.0, 1.0, 4.0])
    R = np.array([0.05, 0.05])

    # initialize problem and GRAMPC
    problem = Crane2D(Q, R, 0.2, 1.25, 0.3)
    grampc = Grampc(problem, options)

    # estimate penaltyMin and set option
    grampc.estim_penmin(True)
    grampc.print_opts()
    grampc.print_params()

    # construct solution structure
    vec = GrampcResults(grampc, Tsim, plot_results=True, plot_statistics=True)

    dt = grampc.param.dt

    for i, t in enumerate(vec.t):
        vec.CPUtime[i] = grampc.run()
        vec.update(grampc, i)

        if i + 1 > len(vec.t):
            break

        # simulate system
        sol = solve_ivp(grampc.ffct, [t, t+dt], grampc.param.x0,
                        args=(grampc.sol.unext, grampc.sol.pnext))

        # set current time and state
        grampc.set_param({"x0": sol.y[:, -1],
                          "t0": t + dt})

        # plots of the grampc predictions
        if i % plotSteps == 0:
            grampc.plot()
            vec.plot()

    vec.plot()
    plt.show()
