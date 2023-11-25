from pygrampc import ProblemBase


class MyProblem(ProblemBase):
    def __init__(self):
        ProblemBase.__init__(self)
        self.Nx = 0
        self.Nu = 0
        self.Np = 0
        self.Ng = 0
        self.Nh = 0
        self.NgT = 0
        self.NhT = 0

    def ffct(self, out, t, x, u, p):
        """
        :param out: only modify with element access like out[:] or out[0]!
        :param t, x, u, p: immutable arrays on C++ side, do not modify!
        """
        pass

    def dfdx_vec(self, out, t, x, vec, u, p):
        pass

    def dfdu_vec(self, out, t, x, vec, u, p):
        pass

    def lfct(self, out, t, x, u, p, xdes, udes):
        pass

    def dldu(self, out, t, x, u, p, xdes, udes):
        pass

    def Vfct(self, out, T, x, p, xdes):
        pass

    def dVdT(self, out, T, x, p, xdes):
        pass

    def gfct(self, out, t, x, u, p):
        pass

    def dgdx_vec(self, out, t, x, u, p, vec):
        pass

    def dgdu_vec(self, out, t, x, u, p, vec):
        pass

    def dgdp_vec(self, out, t, x, u, p, vec):
        pass

    def hfct(self, out, t, x, u, p):
        pass

    def dhdx_vec(self, out, t, x, u, p, vec):
        pass

    def dhdu_vec(self, out, t, x, u, p, vec):
        pass

    def dhdp_vec(self, out, t, x, u, p, vec):
        pass

    def gTfct(self, out, T, x, p):
        pass

    def dgTdx_vec(self, out, T, x, p, vec):
        pass

    def dgTdp_vec(self, out, T, x, p, vec):
        pass

    def dgTdT_vec(self, out, T, x, p, vec):
        pass

    def hTfct(self, out, T, x, p):
        pass

    def dhTdx_vec(self, out, T, x, p, vec):
        pass

    def dhTdp_vec(self, out, T, x, p, vec):
        pass

    def dhTdT_vec(self, out, T, x, p, vec):
        pass

    def dfdx(self, out, t, x, u, p):
        pass

    def dfdxtrans(self, out, t, x, u, p):
        pass

    def dfdt(self, out, t, x, u, p):
        pass

    def dHdxdt(self, out, t, x, u, vec, p):
        pass

    def Mfct(self, out):
        pass

    def Mtrans(self, out):
        pass
