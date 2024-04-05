import numpy as np

class GrampcBinding:
    class _grampc_opt:
        def __init__(self) -> None: ...
        AugLagUpdateGradientRelTol: float
        ConstraintsAbsTol: np.ndarray
        ConstraintsHandling: int
        ConvergenceCheck: int
        ConvergenceGradientRelTol: float
        EqualityConstraints: int
        FlagsRodas: list[int]
        InequalityConstraints: int
        IntegralCost: int
        Integrator: int
        IntegratorAbsTol: float
        IntegratorCost: int
        IntegratorMaxSteps: int
        IntegratorMinStepSize: float
        IntegratorRelTol: float
        JScale: float
        LineSearchAdaptAbsTol: float
        LineSearchAdaptFactor: float
        LineSearchExpAutoFallback: int
        LineSearchInit: float
        LineSearchIntervalFactor: float
        LineSearchIntervalTol: float
        LineSearchMax: float
        LineSearchMin: float
        LineSearchType: int
        MaxGradIter: int
        MaxMultIter: int
        MultiplierDampingFactor: float
        MultiplierMax: float
        Nhor: int
        OptimControl: int
        OptimParam: int
        OptimParamLineSearchFactor: float
        OptimTime: int
        OptimTimeLineSearchFactor: float
        PenaltyDecreaseFactor: float
        PenaltyIncreaseFactor: float
        PenaltyIncreaseThreshold: float
        PenaltyMax: float
        PenaltyMin: float
        ScaleProblem: int
        ShiftControl: int
        TOffset: float
        TScale: float
        TerminalCost: int
        TerminalEqualityConstraints: int
        TerminalInequalityConstraints: int
        TimeDiscretization: int
        cScale: np.ndarray
        pOffset: np.ndarray
        pScale: np.ndarray
        uOffset: np.ndarray
        uScale: np.ndarray
        xOffset: np.ndarray
        xScale: np.ndarray

    class _grampc_param:
        def __init__(self) -> None: ...
        Nc: int
        Ng: int
        NgT: int
        Nh: int
        NhT: int
        Np: int
        Nu: int
        Nx: int
        Thor: float
        Tmax: float
        Tmin: float
        dt: float
        p0: np.ndarray
        pmax: np.ndarray
        pmin: np.ndarray
        t0: float
        u0: np.ndarray
        udes: np.ndarray
        umax: np.ndarray
        umin: np.ndarray
        x0: np.ndarray
        xdes: np.ndarray

    class _grampc_rws:
        def __init__(self) -> None: ...
        T: float
        Tprev: float
        adj: np.ndarray
        cfct: np.ndarray
        cfctAbsTol: np.ndarray
        cfctprev: np.ndarray
        dcdp: np.ndarray
        dcdt: float
        dcdu: np.ndarray
        dcdx: np.ndarray
        gradT: float
        gradp: np.ndarray
        gradpprev: np.ndarray
        gradu: np.ndarray
        graduprev: np.ndarray
        iparRodas: list[int]
        iworkRodas: list[int]
        liworkRodas: int
        lrwsGeneral: int
        lsAdapt: np.ndarray
        lsExplicit: np.ndarray
        lworkRodas: int
        mult: np.ndarray
        p: np.ndarray
        pen: np.ndarray
        pls: np.ndarray
        pprev: np.ndarray
        rparRodas: np.ndarray
        rwsGeneral: np.ndarray
        rwsScale: np.ndarray
        t: np.ndarray
        tls: np.ndarray
        u: np.ndarray
        uls: np.ndarray
        uprev: np.ndarray
        workRodas: np.ndarray
        x: np.ndarray

    class _grampc_sol:
        def __init__(self) -> None: ...
        J: np.ndarray
        Tnext: float
        cfct: float
        iter: list[int]
        pen: float
        pnext: np.ndarray
        status: int
        unext: np.ndarray
        xnext: np.ndarray

    problem: ProblemBase
    opt: _grampc_opt
    param: _grampc_param
    rws: _grampc_rws
    sol: _grampc_sol

    def __init__(self, problem: ProblemBase) -> None: ...

    def run(self) -> float:
        """
        Calls grampc_run which executes one MPC step. Updates the values inside rws and sol.

        Returns:
            float: CPU wall clock time of one function call in milliseconds.

        Raises:
            RuntimeError: if dt or Thor don't have valid values
        """
        
    def estim_penmin(self, run_grampc: bool) -> None:
        """
        Estimates the minimal penalty parameter value. 
        
        Args: 
            run_grampc (bool): Specifies if grampc_run() shall be called.
        
        Raises:
            RuntimeError: if dt or Thor don't have valid values
        """

    def set_rws_u(self, u_new: np.ndarray) -> None:
        """
        Sets the data in rws.u with u_new

        Args:
            u_new (np.ndarray): New data.

        Raises:
            ValueError: if the dimensions of u_new don't match rws.u
        """

    def set_rws_multiplier(self, multiplier_new: np.ndarray) -> None:
        """
        Sets the data in rws.mult with multiplier_new

        Args:
            multiplier_new (np.ndarray): New data.

        Raises:
            ValueError: if the dimensions of multiplier_new don't match rws.mult
        """

    def set_rws_penalty(self, penalty_new: np.ndarray) -> None:
        """
        Sets the data in rws.pen with penalty_new

        Args:
            penalty_new (np.ndarray): New data.

        Raises:
            ValueError: if the dimensions of penalty_new don't match rws.pen
        """

    def ffct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> np.ndarray: ...
    def lfct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, xdes: np.ndarray, udes: np.ndarray) -> float: ...
    def Vfct(self, T: float, x: np.ndarray, p: np.ndarray, xdes: np.ndarray) -> float: ...

    def gTfct(self, T: float, x: np.ndarray, p: np.ndarray) -> np.ndarray: ...
    def gfct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> np.ndarray: ...
    def hTfct(self, T: float, x: np.ndarray, p: np.ndarray) -> np.ndarray: ...
    def hfct(self, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> np.ndarray: ...

    def print_opts(self) -> None:
        """
        Prints the options of the underlying grampc_opt struct.
        """

    def print_params(self) -> None:
        """
        Prints the parameters of the underlying grampc_param struct.
        """

    def print_status(self) -> None:
        """
        Prints the current status of GRAMPC.
        """

class ProblemBase:
    Ng: int
    NgT: int
    Nh: int
    NhT: int
    Np: int
    Nu: int
    Nx: int

    def __init__(self) -> None: ...
    def ffct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dfdx_vec(self, out: np.ndarray, t: float, x: np.ndarray, vec: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dfdu_vec(self, out: np.ndarray, t: float, x: np.ndarray, vec: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dfdp_vec(self, out: np.ndarray, t: float, x: np.ndarray, vec: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...

    def lfct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, xdes: np.ndarray, udes: np.ndarray) -> None: ...
    def dldp(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, xdes: np.ndarray, udes: np.ndarray) -> None: ...
    def dldu(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, xdes: np.ndarray, udes: np.ndarray) -> None: ...
    def dldx(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, xdes: np.ndarray, udes: np.ndarray) -> None: ...

    def Vfct(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, xdes: np.ndarray) -> None: ...
    def dVdT(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, xdes: np.ndarray) -> None: ...
    def dVdp(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, xdes: np.ndarray) -> None: ...
    def dVdx(self, out: np.ndarray, T: float, x: np.ndarray, p: np.ndarray, xdes: np.ndarray) -> None: ...

    def gfct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dgdp_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray) -> None: ...
    def dgdu_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray) -> None: ...
    def dgdx_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray) -> None: ...

    def hfct(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dhdp_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray) -> None: ...
    def dhdu_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray) -> None: ...
    def dhdx_vec(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray, vec: np.ndarray) -> None: ...

    def gTfct(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray) -> None: ...
    def dgTdT_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dgTdp_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dgTdx_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...

    def hTfct(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray) -> None: ...
    def dhTdT_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dhTdp_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dhTdx_vec(self, out: np.ndarray, T: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...

    def dfdx(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dfdxtrans(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dfdt(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, p: np.ndarray) -> None: ...
    def dHdxdt(self, out: np.ndarray, t: float, x: np.ndarray, u: np.ndarray, vec: np.ndarray, p: np.ndarray) -> None: ...
    def Mfct(self, out: np.ndarray) -> None: ...
    def Mtrans(self, out: np.ndarray) -> None: ...
