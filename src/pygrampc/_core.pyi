import numpy.typing as npt

class GrampcBinding:
    class _grampc_opt:
        def __init__(self) -> None: ...
        AugLagUpdateGradientRelTol: float
        ConstraintsAbsTol: npt.NDArray
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
        cScale: npt.NDArray
        pOffset: npt.NDArray
        pScale: npt.NDArray
        uOffset: npt.NDArray
        uScale: npt.NDArray
        xOffset: npt.NDArray
        xScale: npt.NDArray

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
        p0: npt.NDArray
        pmax: npt.NDArray
        pmin: npt.NDArray
        t0: float
        u0: npt.NDArray
        udes: npt.NDArray
        umax: npt.NDArray
        umin: npt.NDArray
        x0: npt.NDArray
        xdes: npt.NDArray

    class _grampc_rws:
        def __init__(self) -> None: ...
        T: float
        Tprev: float
        adj: npt.NDArray
        cfct: npt.NDArray
        cfctAbsTol: npt.NDArray
        cfctprev: npt.NDArray
        dcdp: npt.NDArray
        dcdt: float
        dcdu: npt.NDArray
        dcdx: npt.NDArray
        gradT: float
        gradp: npt.NDArray
        gradpprev: npt.NDArray
        gradu: npt.NDArray
        graduprev: npt.NDArray
        iparRodas: list[int]
        iworkRodas: list[int]
        liworkRodas: int
        lrwsGeneral: int
        lsAdapt: npt.NDArray
        lsExplicit: npt.NDArray
        lworkRodas: int
        mult: npt.NDArray
        p: npt.NDArray
        pen: npt.NDArray
        pls: npt.NDArray
        pprev: npt.NDArray
        rparRodas: npt.NDArray
        rwsGeneral: npt.NDArray
        rwsScale: npt.NDArray
        t: npt.NDArray
        tls: npt.NDArray
        u: npt.NDArray
        uls: npt.NDArray
        uprev: npt.NDArray
        workRodas: npt.NDArray
        x: npt.NDArray

    class _grampc_sol:
        def __init__(self) -> None: ...
        J: npt.NDArray
        Tnext: float
        cfct: float
        iter: list[int]
        pen: float
        pnext: npt.NDArray
        status: int
        unext: npt.NDArray
        xnext: npt.NDArray

    problem: ProblemBase
    opt: _grampc_opt
    param: _grampc_param
    rws: _grampc_rws
    sol: _grampc_sol

    def __init__(self, problem: ProblemBase) -> None: ...

    def run(self) -> float: ...
    def estim_penmin(self, run_grampc: bool) -> None: ...

    def ffct(self, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> npt.NDArray: ...
    def lfct(self, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray, udes: npt.NDArray) -> float: ...
    def Vfct(self, T: float, x: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray) -> float: ...

    def gTfct(self, t: float, x: npt.NDArray, p: npt.NDArray) -> npt.NDArray: ...
    def gfct(self, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> npt.NDArray: ...
    def hTfct(self, t: float, x: npt.NDArray, p: npt.NDArray) -> npt.NDArray: ...
    def hfct(self, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> npt.NDArray: ...

    def print_opts(self) -> None: ...
    def print_params(self) -> None: ...
    def print_status(self) -> None: ...

class ProblemBase:
    Ng: int
    NgT: int
    Nh: int
    NhT: int
    Np: int
    Nu: int
    Nx: int

    def __init__(self) -> None: ...
    def ffct(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dfdx_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, vec: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dfdu_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, vec: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dfdp_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, vec: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...

    def lfct(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray, udes: npt.NDArray) -> None: ...
    def dldp(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray, udes: npt.NDArray) -> None: ...
    def dldu(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray, udes: npt.NDArray) -> None: ...
    def dldx(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray, udes: npt.NDArray) -> None: ...

    def Vfct(self, out: npt.NDArray, T: float, x: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray) -> None: ...
    def dVdT(self, out: npt.NDArray, T: float, x: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray) -> None: ...
    def dVdp(self, out: npt.NDArray, T: float, x: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray) -> None: ...
    def dVdx(self, out: npt.NDArray, T: float, x: npt.NDArray, p: npt.NDArray, xdes: npt.NDArray) -> None: ...

    def gfct(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dgdp_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, vec: npt.NDArray) -> None: ...
    def dgdu_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, vec: npt.NDArray) -> None: ...
    def dgdx_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, vec: npt.NDArray) -> None: ...

    def hfct(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dhdp_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, vec: npt.NDArray) -> None: ...
    def dhdu_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, vec: npt.NDArray) -> None: ...
    def dhdx_vec(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray, vec: npt.NDArray) -> None: ...

    def gTfct(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray) -> None: ...
    def dgTdT_vec(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dgTdp_vec(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dgTdx_vec(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...

    def hTfct(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray) -> None: ...
    def dhTdT_vec(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dhTdp_vec(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dhTdx_vec(self, out: npt.NDArray, T: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...

    def dfdx(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dfdxtrans(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dfdt(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, p: npt.NDArray) -> None: ...
    def dHdxdt(self, out: npt.NDArray, t: float, x: npt.NDArray, u: npt.NDArray, vec: npt.NDArray, p: npt.NDArray) -> None: ...
    def Mfct(self, out: npt.NDArray) -> None: ...
    def Mtrans(self, out: npt.NDArray) -> None: ...
