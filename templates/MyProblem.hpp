/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#include "problem_description.hpp"

using namespace grampc;

class MyProblem : public ProblemDescription
{
    public:
        // Define own variables here
    public:

        MyProblem(); // if needed, define your custom constructor

        ~MyProblem() {}

		/** System function f(t,x,u,p)
		------------------------------------ **/
		virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
		virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
		virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
		virtual void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override;


		/** Integral cost l(t,x(t),u(t),p,xdes,udes)
		-------------------------------------------------- **/
		virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
		/** Gradient dl/dx **/
		virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
		/** Gradient dl/du **/
		virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;
		/** Gradient dl/dp **/
		virtual void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override;


		/** Terminal cost V(T,x,p) */
		virtual void Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
		/** Gradient dV/dx **/
		virtual void dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
		/** Gradient dV/dp **/
		virtual void dVdp(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;
		/** Gradient dV/dT **/
		virtual void dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override;


		/** Equality constraints g(t,x,u,p) = 0 */
		virtual void gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
		virtual void dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;


		/** Inequality constraints h(t,x,u,p) < 0 */
		virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
		virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override;


		/** Terminal equality constraints gT(T,x,p) = 0 */
		virtual void gTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p) override;
		/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
		virtual void dgTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
		virtual void dgTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
		virtual void dgTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;


		/** Terminal inequality constraints hT(T,x,p) < 0 */
		virtual void hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p) override;
		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		virtual void dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
		virtual void dhTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		virtual void dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override;


		/** Additional functions required for semi-implicit systems
		M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
		------------------------------------------------------- **/
		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian df/dt **/
		virtual void dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override;
		/** Jacobian d(dH/dx)/dt  **/
		virtual void dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef vec, VectorConstRef p) override;
		/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mfct(VectorRef out) override;
		/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mtrans(VectorRef out) override;
};