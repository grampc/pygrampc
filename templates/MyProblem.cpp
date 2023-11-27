/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#include "MyProblem.hpp"
#include <cmath>

// insert your problem description here

MyProblem::MyProblem()
 : ProblemBase()
{
    Nx_ = 2;
    Nu_ = 1;
    Np_ = 0;
    Ng_ = 0;
    Nh_ = 0;
    NgT_ = 2;
    NhT_ = 0;
}

/** System function f(t,x,u,p)
------------------------------------ **/
void MyProblem::ffct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {};
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void MyProblem::dfdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) {};
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void MyProblem::dfdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) {};
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void MyProblem::dfdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) {};


/** Integral cost l(t,x(t),u(t),p,xdes,udes)
-------------------------------------------------- **/
void MyProblem::lfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {};
/** Gradient dl/dx **/
void MyProblem::dldx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {};
/** Gradient dl/du **/
void MyProblem::dldu(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {};
/** Gradient dl/dp **/
void MyProblem::dldp(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {};


/** Terminal cost V(T,x,p) */
void MyProblem::Vfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {};
/** Gradient dV/dx **/
void MyProblem::dVdx(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {};
/** Gradient dV/dp **/
void MyProblem::dVdp(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {};
/** Gradient dV/dT **/
void MyProblem::dVdT(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {};


/** Equality constraints g(t,x,u,p) = 0 */
void MyProblem::gfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {};
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void MyProblem::dgdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {};
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void MyProblem::dgdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {};
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void MyProblem::dgdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {};


/** Inequality constraints h(t,x,u,p) < 0 */
void MyProblem::hfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {};
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void MyProblem::dhdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {};
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void MyProblem::dhdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {};
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void MyProblem::dhdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {};


/** Terminal equality constraints gT(T,x,p) = 0 */
void MyProblem::gTfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p) {};
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void MyProblem::dgTdx_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {};
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void MyProblem::dgTdp_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {};
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void MyProblem::dgTdT_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {};


/** Terminal inequality constraints hT(T,x,p) < 0 */
void MyProblem::hTfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p) {};
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void MyProblem::dhTdx_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {};
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void MyProblem::dhTdp_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {};
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void MyProblem::dhTdT_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {};


/** Additional functions required for semi-implicit systems
M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
------------------------------------------------------- **/
/** Jacobian df/dx in cVectorRef form (column-wise) **/
void MyProblem::dfdx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {};
/** Jacobian df/dx in vector form (column-wise) **/
void MyProblem::dfdxtrans(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {};
/** Jacobian df/dt **/
void MyProblem::dfdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {};
/** Jacobian d(dH/dx)/dt  **/
void MyProblem::dHdxdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef vec, cVectorRef p) {};
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
void MyProblem::Mfct(VectorRef out) {};
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void MyProblem::Mtrans(VectorRef out) {};

PYBIND11_MODULE(my_problem, m)
{
    pybind11::class_<MyProblem, ProblemBase>(m, "MyProblem")
        .def(pybind11::init<>())
        .def_readonly("Nx", &MyProblem::Nx_)
        .def_readonly("Nu", &MyProblem::Nu_)
        .def_readonly("Np", &MyProblem::Np_)
        .def_readonly("Ng", &MyProblem::Ng_)
        .def_readonly("Nh", &MyProblem::Nh_)
        .def_readonly("NgT", &MyProblem::NgT_)
        .def_readonly("NhT", &MyProblem::NhT_)

    // these functions are not needed for the Grampc interface, but provide an interface for python code
        .def("ffct", &MyProblem::ffct)
        .def("dfdx_vec", &MyProblem::dfdx_vec)
        .def("dfdu_vec", &MyProblem::dfdu_vec)
        .def("dfdp_vec", &MyProblem::dfdp_vec)

        .def("lfct", &MyProblem::lfct)
        .def("dldx", &MyProblem::dldx)
        .def("dldu", &MyProblem::dldu)
        .def("dldp", &MyProblem::dldp)

        .def("Vfct", &MyProblem::Vfct)
        .def("dVdx", &MyProblem::dVdx)
        .def("dVdp", &MyProblem::dVdp)
        .def("dVdT", &MyProblem::dVdT)

        .def("gfct", &MyProblem::gfct)
        .def("dgdx_vec", &MyProblem::dgdx_vec)
        .def("dgdu_vec", &MyProblem::dgdu_vec)
        .def("dgdp_vec", &MyProblem::dgdp_vec)

        .def("hfct", &MyProblem::hfct)
        .def("dhdx_vec", &MyProblem::dhdx_vec)
        .def("dhdu_vec", &MyProblem::dhdu_vec)
        .def("dhdp_vec", &MyProblem::dhdp_vec)

        .def("gTfct", &MyProblem::gTfct)
        .def("dgTdx_vec", &MyProblem::dgTdx_vec)
        .def("dgTdp_vec", &MyProblem::dgTdp_vec)
        .def("dgTdT_vec", &MyProblem::dgTdT_vec)

        .def("hTfct", &MyProblem::hTfct)
        .def("dhTdx_vec", &MyProblem::dhTdx_vec)
        .def("dhTdp_vec", &MyProblem::dhTdp_vec)
        .def("dhTdT_vec", &MyProblem::dhTdT_vec)

        .def("dfdx", &MyProblem::dfdx)
        .def("dfdxtrans", &MyProblem::dfdxtrans)
        .def("dfdt", &MyProblem::dfdt)
        .def("dHdxdt", &MyProblem::dHdxdt)
        .def("Mfct", &MyProblem::Mfct)
        .def("Mtrans", &MyProblem::Mtrans);
}