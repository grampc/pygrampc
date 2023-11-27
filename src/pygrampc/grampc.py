import json
from math import nan, inf
from typing import NamedTuple
import matplotlib.pyplot as plt
import numpy as np

from ._core import GrampcBinding, ProblemBase


class ValidatorOptions(NamedTuple):
    allowedValues: tuple
    Class: str
    dimension: int = 0


class Grampc(GrampcBinding):
    """
    An interface for the GRAMPC solver. Handles the initialization and destruction of the GRAMPC solver.

    Attributes:
        problem: The problem description for this instance.
        fig: handle for the prediction plot figure.
        param: The GRAMPC parameters struct.
        opt: The GRAMPC options struct.
        rws: The GRAMPC real workspace struct.
        sol: The GRAMPC solution struct.
    """
    _parameters_real = {
        "dt": ValidatorOptions((0.0, inf), "OpenInterval"),
        "Thor": ValidatorOptions((0.0, inf), "OpenInterval"),
        "Tmax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "Tmin": ValidatorOptions((0.0, inf), "OpenInterval"),
        "t0": ValidatorOptions((-inf, inf), "OpenInterval")}

    _options_real = {
        "IntegratorRelTol": ValidatorOptions((0.0, inf), "OpenInterval"),
        "IntegratorAbsTol": ValidatorOptions((0.0, inf), "OpenInterval"),
        "IntegratorMinStepSize": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchMax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchMin": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchInit": ValidatorOptions((0.0, inf), "OpenInterval"),
        "LineSearchAdaptAbsTol": ValidatorOptions((0.0, inf), "LeftClosedInterval"),
        "LineSearchAdaptFactor": ValidatorOptions((1.0, inf), "OpenInterval"),
        "LineSearchIntervalTol": ValidatorOptions((0.0, 0.5), "OpenInterval"),
        "LineSearchIntervalFactor": ValidatorOptions((0.0, 1.0), "OpenInterval"),
        "OptimParamLineSearchFactor": ValidatorOptions((0.0, inf), "OpenInterval"),
        "OptimTimeLineSearchFactor": ValidatorOptions((0.0, inf), "OpenInterval"),
        "TScale": ValidatorOptions((0.0, inf), "OpenInterval"),
        "TOffset": ValidatorOptions((-inf, inf), "OpenInterval"),
        "JScale": ValidatorOptions((0.0, inf), "OpenInterval"),
        "MultiplierMax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "MultiplierDampingFactor": ((0.0, 1.0), "ClosedInterval"),
        "PenaltyMax": ValidatorOptions((0.0, inf), "OpenInterval"),
        "PenaltyMin": ValidatorOptions((0.0, inf), "OpenInterval"),
        "PenaltyIncreaseFactor": ValidatorOptions((1.0, inf), "LeftClosedInterval"),
        "PenaltyDecreaseFactor": ValidatorOptions((0.0, 1.0), "ClosedInterval"),
        "PenaltyIncreaseThreshold": ValidatorOptions((0.0, inf), "LeftClosedInterval"),
        "AugLagUpdateGradientRelTol": ValidatorOptions((0.0, 1.0), "ClosedInterval"),
        "ConvergenceGradientRelTol": ValidatorOptions((0.0, 1.0), "ClosedInterval")}

    _options_int = {
        "Nhor": ValidatorOptions((2, inf), "LeftClosedInterval"),
        "MaxGradIter": ValidatorOptions((1, inf), "LeftClosedInterval"),
        "MaxMultIter": ValidatorOptions((1, inf), "LeftClosedInterval"),
        "IntegratorMaxSteps": ValidatorOptions((1, inf), "LeftClosedInterval")}

    _options_str = {
        "ShiftControl": ValidatorOptions(("on", "off"), "Enum"),
        "IntegralCost": ValidatorOptions(("on", "off"), "Enum"),
        "TerminalCost": ValidatorOptions(("on", "off"), "Enum"),
        "IntegratorCost": ValidatorOptions(("trapezodial", "simpson"), "Enum"),
        "Integrator": ValidatorOptions(("euler", "modeuler", "heun", "ruku45", "rodas"), "Enum"),
        "LineSearchType": ValidatorOptions(("adaptive", "explicit1", "explicit2"), "Enum"),
        "LineSearchExpAutoFallback": ValidatorOptions(("on", "off"), "Enum"),
        "OptimControl": ValidatorOptions(("on", "off"), "Enum"),
        "OptimParam": ValidatorOptions(("on", "off"), "Enum"),
        "OptimTime": ValidatorOptions(("on", "off"), "Enum"),
        "ScaleProblem": ValidatorOptions(("on", "off"), "Enum"),
        "EqualityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "InequalityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "TerminalEqualityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "TerminalInequalityConstraints": ValidatorOptions(("on", "off"), "Enum"),
        "ConstraintsHandling": ValidatorOptions(("auglag", "extpen"), "Enum"),
        "ConvergenceCheck": ValidatorOptions(("on", "off"), "Enum")}

    def __init__(self, problem: ProblemBase, options_path: str = None, plot_prediction=False):
        """
        Initilizes the interface and the underlying GRAMPC object.

        Args:
            problem: The concrete problem description object.
            options_path (str): The path to the options.json file.
            plot_prediction (bool): If True enables the prediction plot utility. 
        """
        super().__init__(problem)
        self._parameters_real_vec = {
            "x0": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "xdes": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "u0": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "udes": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "umax": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "umin": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "p0": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "pmax": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "pmin": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np)}
        self._options_real_vec = {
            "xScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "xOffset": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nx),
            "uScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "uOffset": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nu),
            "pScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "pOffset": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Np),
            "cScale": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nc),
            "ConstraintsAbsTol": ValidatorOptions((-inf, inf), "OpenInterval", self.param.Nc)}

        self._options_int_vec = {"FlagsRodas": ValidatorOptions((0, self.param.Nx), "Rodas", 8)}

        if options_path is not None:
            with open(options_path) as f:
                grampc_options = json.load(f)
            self.set_param(grampc_options["Parameters"])
            self.set_opt(grampc_options["Options"])

        if plot_prediction:
            self.fig = plt.figure(figsize=(10, 5))
        else:
            self.fig = None

    @staticmethod
    def _len(value):
        """Returns the length of a list or numpy array. For scalar inputs the function returns 0"""
        try:
            return len(value)
        except TypeError:
            return 0

    def _validator(self, key: str, value, options: ValidatorOptions) -> bool:
        """
        Helper function to validate the inputs for the GRAMPC parameters and options.

        Args:
            key (str): The option/parameter key to check.
            value: The value of the option/parameter.
            options (ValidatorOptions): NamedTuple with the option class and allowed values.

        Returns:
            bool: True if value is valid.
        
        Raises:
            ValueError: Raised if value is not valid.
            ValueError: Raised if dimension does not match.
        """
        valid_input = False
        if options.Class == "Enum":
            valid_input = value in options.allowedValues
        elif options.Class == "OpenInterval":
            if self._len(value) == options.dimension:
                if isinstance(value, (list, np.ndarray)):
                    valid_input = all(
                        [options.allowedValues[0] < item < options.allowedValues[1] for item in value])
                else:
                    valid_input = options.allowedValues[0] < value < options.allowedValues[1]
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected {options.dimension}, "
                    f"got {self._len(value)}")

        elif options.Class == "LeftClosedInterval":
            if self._len(value) == options.dimension:
                if isinstance(value, (list, np.ndarray)):
                    valid_input = all(
                        [options.allowedValues[0] <= item < options.allowedValues[1] for item in value])
                else:
                    valid_input = options.allowedValues[0] <= value < options.allowedValues[1]
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected {options.dimension}, "
                    f"got {self._len(value)}")

        elif options.Class == "ClosedInterval":
            if self._len(value) == options.dimension:
                if isinstance(value, (list, np.ndarray)):
                    valid_input = all(
                        [options.allowedValues[0] <= item <= options.allowedValues[1] for item in value])
                else:
                    valid_input = options.allowedValues[0] <= value <= options.allowedValues[1]
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected {options.dimension}, "
                    f"got {self._len(value)}")

        elif options.Class == "Rodas":
            if self._len(value) == 8:
                # Check if IFCN, IDFX, IJAC and IMAS are either 1 or 0
                valid_input = all([value[i] == 1 or value[i] == 0 for i in range(4)])
                # check if MLJAC, MUJAX, MLMAS and MUMAS are between 0 and Nx
                valid_input = valid_input and all(
                    [options.allowedValues[0] <= value[i + 4] <= options.allowedValues[1] for i in range(4)])
            else:
                raise ValueError(
                    f"{key}: dimension does not match problem description. Expected 8, got {self._len(value)}")

        if not valid_input:
            raise ValueError(f"{key}: unallowed value detected. Expected {options.allowedValues}, got {value}")
        return valid_input

    def set_param(self, parameters: dict):
        """
        Sets the GRAMPC parameters with value validation.

        Args:
            parameters (dict): Mapping with parameter name and values.
        
        Raises: 
            ValueError: if the dimensions of vector values like the states or controls don't match the problem description
            KeyError: if a parameter key is not a valid GRAMPC parameter
        """
        for key, value in parameters.items():
            if key in self._parameters_real:
                if self._validator(key, value, self._parameters_real[key]):
                    self._set_param_real(key, value)
            elif key in self._parameters_real_vec:
                if self._validator(key, value, self._parameters_real_vec[key]):
                    self._set_param_real_vec(key, value)
            else:
                raise KeyError(f"{key} not found as GRAMPC parameter")

    def set_opt(self, options: dict):
        """
        Sets the GRAMPC options with value validation.

        Args:
            options (dict): Mapping with parameter name and values.
        
        Raises: 
            ValueError: if the dimensions of vector values like the states or controls don't match the problem description
            KeyError: if a parameter key is not a valid GRAMPC parameter
        """
        for key, value in options.items():
            if key in self._options_str:
                if self._validator(key, value, self._options_str[key]):
                    self._set_opt_str(key, value)
            elif key in self._options_int:
                if self._validator(key, value, self._options_int[key]):
                    self._set_opt_int(key, value)
            elif key in self._options_int_vec:
                self._set_opt_int_vec(key, value)
            elif key in self._options_real:
                if self._validator(key, value, self._options_real[key]):
                    self._set_opt_real(key, value)
            elif key in self._options_real_vec:
                if self._validator(key, value, self._options_real_vec[key]):
                    self._set_opt_real_vec(key, value)
            else:
                raise KeyError(f"{key} not found as GRAMPC option")

    def plot(self):
        """
        Draws the current prediction plot with the values of the rws struct.
        """
        if self.fig is not None:
            plt.figure(self.fig.number)
            self.fig.clf()
            self.fig.suptitle("Prediction")

            self.fig.add_subplot(2, 3, 1)
            for val in self.rws.x:
                plt.plot(self.rws.t, val)
            plt.xlim(left=self.rws.t[0])
            plt.grid()
            plt.title("Predicted states")

            self.fig.add_subplot(2, 3, 2)
            for val in self.rws.adj:
                plt.plot(self.rws.t, val)
            plt.xlim(left=self.rws.t[0])
            plt.grid()
            plt.title("Predicted adjoint states")

            self.fig.add_subplot(2, 3, 3)
            for val in self.rws.u:
                plt.plot(self.rws.t, val)
            plt.xlim(left=self.rws.t[0])
            plt.grid()
            plt.title("Predicted controls")

            self.fig.add_subplot(2, 3, 4)
            if self.param.Nc == self.param.NgT + self.param.NhT:
                for val in self.rws.cfct:
                    plt.scatter(self.rws.t[-1], val[-1])
            else:
                for val in self.rws.cfct:
                    plt.plot(self.rws.t, val)
            plt.xlim(left=self.rws.t[0])
            plt.grid()
            plt.title("Predicted constraints")

            self.fig.add_subplot(2, 3, 5)
            if self.param.Nc == self.param.NgT + self.param.NhT:
                for val in self.rws.mult:
                    plt.scatter(self.rws.t[-1], val[-1])
            else:
                for val in self.rws.mult:
                    plt.plot(self.rws.t, val)
            plt.xlim(left=self.rws.t[0])
            plt.grid()
            plt.title("Predicted Lagrange multipliers")

            self.fig.add_subplot(2, 3, 6)
            if self.param.Nc == self.param.NgT + self.param.NhT:
                for val in self.rws.pen:
                    plt.scatter(self.rws.t[-1], val[-1])
            else:
                for val in self.rws.pen:
                    plt.plot(self.rws.t, val)
            plt.xlim(left=self.rws.t[0])
            plt.grid()
            plt.title("Predicted penalty parameters")

            plt.draw()
            plt.pause(0.0001)


class GrampcResults:
    """
    Dataclass for storing the results of the GRAMPC solver.

    Attributes:
        t (ndarray): time array.
        x (ndarray): state array.
        u (ndarray): control array.
        p (ndarray): parameter array.
        adj (ndarray): adjoint states array.
        J (ndarray): cost array.
        CPUtime (ndarray): CPU wall clock time for every time step
    """

    class index:
        def __init__(self, grampc):
            self.g = grampc.param.Ng
            self.h = grampc.param.Nh + self.g
            self.gT = grampc.param.NgT + self.h
            self.hT = grampc.param.NhT + self.gT

    def __init__(self, grampc: Grampc, Tsim: float, plot_statistics: bool = False, plot_results: bool = False):
        num_data_points = int(Tsim / grampc.param.dt) + 1
        self._id = self.index(grampc)
        self.t = np.linspace(grampc.param.t0, Tsim, num_data_points)
        self.x = np.full((num_data_points, grampc.param.Nx), nan)
        self.u = np.full((num_data_points, grampc.param.Nu), nan)
        self.p = np.full((num_data_points, grampc.param.Nu), nan)
        self.adj = np.full((num_data_points, grampc.param.Nx), nan)
        self.J = np.full((num_data_points, 2), nan)

        self.CPUtime = np.full((num_data_points, 1), nan)
        self.lsExplicit = np.full((num_data_points, 1), nan)
        self.continuous_constraints = grampc.param.Ng + grampc.param.Nh > 0

        if grampc.param.Ng + grampc.param.Nh > 0:
            self.constr = np.full((num_data_points, grampc.param.Ng + grampc.param.Nh), nan)
            self.mult = np.full((num_data_points, grampc.param.Ng + grampc.param.Nh), nan)
            self.pen = np.full((num_data_points, grampc.param.Ng + grampc.param.Nh), nan)

        if grampc.param.NgT + grampc.param.NhT > 0:
            self.constrT = np.full((num_data_points, grampc.param.NgT + grampc.param.NhT), nan)
            self.multT = np.full((num_data_points, grampc.param.NgT + grampc.param.NhT), nan)
            self.penT = np.full((num_data_points, grampc.param.NgT + grampc.param.NhT), nan)

        if plot_statistics:
            self.fig_statistics = plt.figure()
        else:
            self.fig_statistics = None

        if plot_results:
            self.fig_results = plt.figure()
        else:
            self.fig_results = None

    def update(self, grampc: Grampc, index: int):
        self.x[index, :] = grampc.rws.x[:, 0] * grampc.opt.xScale + grampc.opt.xOffset
        self.u[index, :] = grampc.rws.u[:, 0] * grampc.opt.uScale + grampc.opt.uOffset
        self.adj[index, :] = grampc.rws.adj[:, 0] / grampc.opt.xScale
        self.J[index, :] = grampc.sol.J
        if grampc.opt.LineSearchType != 0:
            self.lsExplicit[index] = grampc.rws.lsExplicit[2]

        # retrieve the equality and inequality constraints
        if grampc.param.Ng + grampc.param.Nh > 0:
            if grampc.param.Ng:
                self.constr[index, 0:self._id.g] = grampc.gfct(self.t[index], self.x[index, :], self.u[index, :], self.p[index, :])
            if grampc.param.Nh:
                self.constr[index, self._id.g:self._id.h] = grampc.hfct(self.t[index], self.x[index, :], self.u[index, :], self.p[index, :])

            self.mult[index, :] = grampc.rws.mult[0:self._id.h, 0]
            self.pen[index, :] = grampc.rws.pen[0:self._id.h, 0]

        if grampc.param.NgT + grampc.param.NhT > 0:
            if grampc.param.NgT:
                self.constr[index, self._id.h:self._id.gT] = grampc.gTfct(self.t[index], self.x[index, :], self.p[index, :])
            if grampc.param.NhT:
                self.constr[index, self._id.gT:self._id.hT] = grampc.hTfct(self.t[index], self.x[index, :], self.p[index, :])

            self.multT[index, :] = grampc.rws.mult[self._id.h:self._id.hT, -1]
            self.penT[index, :] = grampc.rws.pen[self._id.h:self._id.gT, -1]

    def plot(self):
        if self.fig_statistics is not None:
            plt.figure(self.fig_statistics.number)
            self.fig_statistics.clf()
            self.fig_statistics.suptitle("Statistics")

            self.fig_statistics.add_subplot(1, 3, 1)
            plt.plot(self.t, self.J)
            plt.xlim(self.t[0], self.t[-1])
            plt.title("Costs")
            plt.grid()

            self.fig_statistics.add_subplot(1, 3, 2)
            plt.plot(self.t, self.CPUtime)
            plt.xlim(self.t[0], self.t[-1])
            plt.title("Computation time in ms")
            plt.grid()

            self.fig_statistics.add_subplot(1, 3, 3)
            plt.plot(self.t, self.lsExplicit)
            plt.xlim(self.t[0], self.t[-1])
            plt.title("Line search step size")
            plt.yscale("log")
            plt.grid()

        if self.fig_results is not None:
            plt.figure(self.fig_results.number)
            self.fig_results.clf()
            self.fig_results.suptitle("Results")
            numRows = 1

            if self.continuous_constraints:
                numRows = 2

            self.fig_results.add_subplot(numRows, 3, 1)
            plt.plot(self.t, self.x)
            plt.xlim(self.t[0], self.t[-1])
            plt.title("States")
            plt.grid()

            self.fig_results.add_subplot(numRows, 3, 2)
            plt.plot(self.t, self.adj)
            plt.xlim(self.t[0], self.t[-1])
            plt.title("Adjoint states")
            plt.grid()

            self.fig_results.add_subplot(numRows, 3, 3)
            plt.plot(self.t, self.u)
            plt.xlim(self.t[0], self.t[-1])
            plt.title("Controls")
            plt.grid()

            if self.continuous_constraints:
                self.fig_results.add_subplot(numRows, 3, 4)
                plt.plot(self.t, self.constr)
                plt.xlim(self.t[0], self.t[-1])
                plt.title("Constraints")
                plt.grid()

                self.fig_results.add_subplot(numRows, 3, 5)
                plt.plot(self.t, self.mult)
                plt.xlim(self.t[0], self.t[-1])
                plt.title("Lagrange multipliers")
                plt.grid()

                self.fig_results.add_subplot(numRows, 3, 6)
                plt.plot(self.t, self.pen)
                plt.xlim(self.t[0], self.t[-1])
                plt.title("Penalty parameters")
                plt.grid()

        plt.draw()
        plt.pause(0.0001)
