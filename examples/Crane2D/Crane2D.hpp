/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#include "problem_base.hpp"

class PYBIND11_EXPORT Crane2D : public ProblemBase
{
    public:
        Vector Q;
        Vector R;
        typeRNum ScaleConstraint;
        typeRNum MaxConstraintHeight;
        typeRNum MaxAngularDeflection;
    public:

        Crane2D(Vector Q, Vector R, typeRNum ScaleConstraint, typeRNum MaxConstraintHeight, typeRNum MaxAngularDeflection);

        ~Crane2D() {}

		/** System function f(t,x,u,p)
		------------------------------------ **/
		void ffct(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
		void dfdx_vec(VectorRef out, const double t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override;
		/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
		void dfdu_vec(VectorRef out, const double t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override;


		/** Integral cost l(t,x(t),u(t),p,xdes,udes)
		-------------------------------------------------- **/
		void lfct(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override;
		/** Gradient dl/dx **/
		void dldx(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override;
		/** Gradient dl/du **/
		void dldu(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override;


		/** Inequality constraints h(t,x,u,p) < 0 */
		void hfct(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p) override;
		/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
		void dhdx_vec(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;
		/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
		void dhdu_vec(VectorRef out, const double t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override;
};