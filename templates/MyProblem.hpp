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

//PYBIND11_EXPORT only needed if class shall be extended in Python. Then also a trampoline class is needed.
class PYBIND11_EXPORT MyProblem : public ProblemDescription
{
    public:
        // Define own variables here
    public:

        MyProblem(); // if needed, define your custom constructor

        ~MyProblem() {}

		/** System function f(t,x,u,p)
		------------------------------------ **/
		void ffct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
		void dfdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override;
		/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
		void dfdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override;
		/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
		void dfdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override;


		/** Integral cost l(t,x(t),u(t),p,xdes,udes)
		-------------------------------------------------- **/
		void lfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override;
		/** Gradient dl/dx **/
		void dldx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override;
		/** Gradient dl/du **/
		void dldu(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override;
		/** Gradient dl/dp **/
		void dldp(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override;


		/** Terminal cost V(T,x,p) */
		void Vfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) override;
		/** Gradient dV/dx **/
		void dVdx(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) override;
		/** Gradient dV/dp **/
		void dVdp(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) override;
		/** Gradient dV/dT **/
		void dVdT(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) override;


		/** Equality constraints g(t,x,u,p) = 0 */
		void gfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
		void dgdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
		void dgdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
		void dgdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;


		/** Inequality constraints h(t,x,u,p) < 0 */
		void hfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
		void dhdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
		void dhdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
		void dhdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;


		/** Terminal equality constraints gT(T,x,p) = 0 */
		void gTfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p) override;
		/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
		void dgTdx_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
		void dgTdp_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
		void dgTdT_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) override;


		/** Terminal inequality constraints hT(T,x,p) < 0 */
		void hTfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p) override;
		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		void dhTdx_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
		void dhTdp_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		void dhTdT_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) override;


		/** Additional functions required for semi-implicit systems
		M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
		------------------------------------------------------- **/
		/** Jacobian df/dx in vector form (column-wise) **/
		void dfdx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian df/dx in vector form (column-wise) **/
		void dfdxtrans(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian df/dt **/
		void dfdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian d(dH/dx)/dt  **/
		void dHdxdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef vec, cVectorRef p) override;
		/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
		void Mfct(VectorRef out) override;
		/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
		void Mtrans(VectorRef out) override;
};