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

// insert your problem description here

MyProblem::MyProblem()
 : ProblemDescription()
{
    // Overwrite these parameters defined by ProblemBase
    Nx = 2;
    Nu = 1;
    Np = 0;
    Ng = 0;
    Nh = 0;
    NgT = 2;
    NhT = 0;
}

/** System function f(t,x,u,p)
------------------------------------ **/
virtual void MyProblem::ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {};
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
virtual void MyProblem::dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) {};
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
virtual void MyProblem::dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) {};
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
virtual void MyProblem::dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) {};


/** Integral cost l(t,x(t),u(t),p,xdes,udes)
-------------------------------------------------- **/
virtual void MyProblem::lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {};
/** Gradient dl/dx **/
virtual void MyProblem::dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {};
/** Gradient dl/du **/
virtual void MyProblem::dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {};
/** Gradient dl/dp **/
virtual void MyProblem::dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {};


/** Terminal cost V(T,x,p) */
virtual void MyProblem::Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {};
/** Gradient dV/dx **/
virtual void MyProblem::dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {};
/** Gradient dV/dp **/
virtual void MyProblem::dVdp(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {};
/** Gradient dV/dT **/
virtual void MyProblem::dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {};


/** Equality constraints g(t,x,u,p) = 0 */
virtual void MyProblem::gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {};
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
virtual void MyProblem::dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
virtual void MyProblem::dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
virtual void MyProblem::dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {};


/** Inequality constraints h(t,x,u,p) < 0 */
virtual void MyProblem::hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {};
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
virtual void MyProblem::dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
virtual void MyProblem::dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
virtual void MyProblem::dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {};


/** Terminal equality constraints gT(T,x,p) = 0 */
virtual void MyProblem::gTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p) {};
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
virtual void MyProblem::dgTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
virtual void MyProblem::dgTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
virtual void MyProblem::dgTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {};


/** Terminal inequality constraints hT(T,x,p) < 0 */
virtual void MyProblem::hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p) {};
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
virtual void MyProblem::dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
virtual void MyProblem::dhTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {};
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
virtual void MyProblem::dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {};


/** Additional functions required for semi-implicit systems
M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
------------------------------------------------------- **/
/** Jacobian df/dx in cVectorRef form (column-wise) **/
virtual void MyProblem::dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {};
/** Jacobian df/dx in vector form (column-wise) **/
virtual void MyProblem::dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {};
/** Jacobian df/dt **/
virtual void MyProblem::dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {};
/** Jacobian d(dH/dx)/dt  **/
virtual void MyProblem::dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef vec, VectorConstRef p) {};
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
virtual void MyProblem::Mfct(VectorRef out) {};
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
virtual void MyProblem::Mtrans(VectorRef out) {};

PYBIND11_MODULE(my_problem, m)
{
    pybind11::class_<MyProblem, ProblemDescription>(m, "MyProblem")
        .def(pybind11::init<>())
        .def_readonly("Nx", &MyProblem::Nx)
        .def_readonly("Nu", &MyProblem::Nu)
        .def_readonly("Np", &MyProblem::Np)
        .def_readonly("Ng", &MyProblem::Ng)
        .def_readonly("Nh", &MyProblem::Nh)
        .def_readonly("NgT", &MyProblem::NgT)
        .def_readonly("NhT", &MyProblem::NhT)

    // these functions are not needed for the GRAMPC interface, but provide an interface for python code
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