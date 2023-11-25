#include "Crane2D.hpp"
#include <cmath>

// insert your problem description here

Crane2D::Crane2D(Vector Q, Vector R, typeRNum ScaleConstraint, typeRNum MaxConstraintHeight, typeRNum MaxAngularDeflection)
 : ProblemBase()
{
    Nx_ = 6;
    Nu_ = 2;
    Np_ = 0;
    Ng_ = 0;
    Nh_ = 3;
    NgT_ = 0;
    NhT_ = 0;
    this->Q = Q;
    this->R = R;
    this->ScaleConstraint = ScaleConstraint;
    this->MaxConstraintHeight = MaxConstraintHeight;
    this->MaxAngularDeflection = MaxAngularDeflection;
}

/** System function f(t,x,u,p)
------------------------------------ **/
void Crane2D::ffct(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p) 
{
    out[0] = x[1];
    out[1] = u[0];
    out[2] = x[3];
    out[3] = u[1];
    out[4] = x[5];
    out[5] = -((9.81 * sin(x[4]) + cos(x[4]) * u[0] + 2 * x[3] * x[5]) / x[2]);
};

/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void Crane2D::dfdx_vec(VectorRef out, const double t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) 
{
    typeRNum sinX = sin(x[4]);
    typeRNum cosX = cos(x[4]);
    typeRNum g = 9.81;

    out[0] = 0;
    out[1] = vec[0];
    out[2] = (g * sinX + cosX * u[0] + 2 * x[3] * x[5]) * vec[5] / (x[2]*x[2]);
    out[3] = vec[2] - (2 * x[5] * vec[5]) / x[2];
    out[4] = -((g * cosX - sinX * u[0]) * vec[5] / x[2]);
    out[5] = vec[4] - (2 * x[3] * vec[5]) / x[2];
};

/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void Crane2D::dfdu_vec(VectorRef out, const double t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) 
{
    out[0] = vec[1] - (cos(x[4]) * vec[5]) / x[2];
    out[1] = vec[3];
};


/** Integral cost l(t,x(t),u(t),p,xdes,udes)
-------------------------------------------------- **/
void Crane2D::lfct(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) 
{
    out[0] = 0.0;
    for (typeInt i = 0; i < Q.size(); i++)
    {
        out[0] += Q[i] * pow(x[i] - xdes[i], 2.0);
    }
    for (typeInt i = 0; i < R.size(); i++)
    {
        out[0] += R[i] * pow(u[i] - udes[i], 2.0);
    }
};
/** Gradient dl/dx **/
void Crane2D::dldx(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) 
{
    for (typeInt i = 0; i < Q.size(); i++)
    {
        out[i] = 2 * Q[i] * (x[i] - xdes[i]);
    }
};
/** Gradient dl/du **/
void Crane2D::dldu(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) 
{
    for (typeInt i = 0; i < R.size(); i++)
    {
        out[i] = 2 * R[i] * (u[i] - udes[i]);
    }
};


/** Inequality constraints h(t,x,u,p) < 0 */
void Crane2D::hfct(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p) 
{
    typeRNum Position = x[0] + sin(x[4]) * x[2];

    out[0] = cos(x[4]) * x[2] - ScaleConstraint * pow(Position, 2.0) - MaxConstraintHeight;
    out[1] = x[5] - MaxAngularDeflection;
    out[2] = -x[5] - MaxAngularDeflection;
};

/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void Crane2D::dhdx_vec(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) 
{
    typeRNum tmp = ScaleConstraint * (x[0] + sin(x[4]) * x[2]);

    out[0] = tmp * vec[0];
    out[1] = 0;
    out[2] = (sin(x[4]) * tmp + cos(x[4])) * vec[0];
    out[3] = 0;
    out[4] = (cos(x[4]) * x[2] * tmp - sin(x[4]) * x[2]) * vec[0];
    out[5] = 0 + vec[1] - vec[2];
};

/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void Crane2D::dhdu_vec(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec)
{
    out[0] = 0;
};

PYBIND11_MODULE(crane_problem, m)
{
    pybind11::class_<Crane2D, ProblemBase>(m, "Crane2D")
        .def(pybind11::init<Vector, Vector, typeRNum, typeRNum, typeRNum>())
        .def_readonly("Nx", &Crane2D::Nx_)
        .def_readonly("Nu", &Crane2D::Nu_)
        .def_readonly("Np", &Crane2D::Np_)
        .def_readonly("Ng", &Crane2D::Ng_)
        .def_readonly("Nh", &Crane2D::Nh_)
        .def_readonly("NgT", &Crane2D::NgT_)
        .def_readonly("NhT", &Crane2D::NhT_)
        
        // make your custom fields available from python
        .def_readwrite("Q", &Crane2D::Q)
        .def_readwrite("R", &Crane2D::R)
        .def_readwrite("MaxAngularDeflection", &Crane2D::MaxAngularDeflection)
        .def_readwrite("ScaleConstraint", &Crane2D::ScaleConstraint)
        .def_readwrite("MaxConstraintHeight", &Crane2D::MaxConstraintHeight)

    // these functions are not needed for the Grampc interface, but provide an interface for python code
        .def("ffct", &Crane2D::ffct)
        .def("dfdx_vec", &Crane2D::dfdx_vec)
        .def("dfdu_vec", &Crane2D::dfdu_vec)
        .def("dfdp_vec", &Crane2D::dfdp_vec)

        .def("lfct", &Crane2D::lfct)
        .def("dldx", &Crane2D::dldx)
        .def("dldu", &Crane2D::dldu)
        .def("dldp", &Crane2D::dldp)

        .def("Vfct", &Crane2D::Vfct)
        .def("dVdx", &Crane2D::dVdx)
        .def("dVdp", &Crane2D::dVdp)
        .def("dVdT", &Crane2D::dVdT)

        .def("gfct", &Crane2D::gfct)
        .def("dgdx_vec", &Crane2D::dgdx_vec)
        .def("dgdu_vec", &Crane2D::dgdu_vec)
        .def("dgdp_vec", &Crane2D::dgdp_vec)

        .def("hfct", &Crane2D::hfct)
        .def("dhdx_vec", &Crane2D::dhdx_vec)
        .def("dhdu_vec", &Crane2D::dhdu_vec)
        .def("dhdp_vec", &Crane2D::dhdp_vec)

        .def("gTfct", &Crane2D::gTfct)
        .def("dgTdx_vec", &Crane2D::dgTdx_vec)
        .def("dgTdp_vec", &Crane2D::dgTdp_vec)
        .def("dgTdT_vec", &Crane2D::dgTdT_vec)

        .def("hTfct", &Crane2D::hTfct)
        .def("dhTdx_vec", &Crane2D::dhTdx_vec)
        .def("dhTdp_vec", &Crane2D::dhTdp_vec)
        .def("dhTdT_vec", &Crane2D::dhTdT_vec)

        .def("dfdx", &Crane2D::dfdx)
        .def("dfdxtrans", &Crane2D::dfdxtrans)
        .def("dfdt", &Crane2D::dfdt)
        .def("dHdxdt", &Crane2D::dHdxdt)
        .def("Mfct", &Crane2D::Mfct)
        .def("Mtrans", &Crane2D::Mtrans);
}