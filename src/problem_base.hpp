/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#ifndef PROBLEM_BINDINGS_HPP
#define PROBLEM_BINDINGS_HPP

extern "C"
{
    #include "grampc.h"
}

#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, 1> Vector;
typedef Eigen::Ref<Vector> VectorRef;
typedef Eigen::Ref<const Vector> cVectorRef;

class PYBIND11_EXPORT ProblemBase
{
    public:
        typeInt Nx_;
        typeInt Nu_;
        typeInt Np_;
        typeInt Ng_;
        typeInt Nh_;
        typeInt NgT_;
        typeInt NhT_;
        typeInt Rodas_Jac = 0;
        typeInt Rodas_M = 0;

    public:

        virtual ~ProblemBase() {}

		/** System function f(t,x,u,p,userparam)
		------------------------------------ **/
		virtual void ffct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) = 0;
		/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
		virtual void dfdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) = 0;
		/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
		virtual void dfdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) {}
		/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
		virtual void dfdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) {}


		/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
		-------------------------------------------------- **/
		virtual void lfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {}
		/** Gradient dl/dx **/
		virtual void dldx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {}
		/** Gradient dl/du **/
		virtual void dldu(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {}
		/** Gradient dl/dp **/
		virtual void dldp(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) {}


		/** Terminal cost V(T,x,p) */
		virtual void Vfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {}
		/** Gradient dV/dx **/
		virtual void dVdx(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {}
		/** Gradient dV/dp **/
		virtual void dVdp(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {}
		/** Gradient dV/dT **/
		virtual void dVdT(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef xdes) {}


		/** Equality constraints g(t,x,u,p) = 0 */
		virtual void gfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {}
		/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dgdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
		virtual void dgdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dgdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {}


		/** Inequality constraints h(t,x,u,p) < 0 */
		virtual void hfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {}
		/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
		virtual void dhdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
		virtual void dhdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
		virtual void dhdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) {}


		/** Terminal equality constraints gT(T,x,p) = 0 */
		virtual void gTfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p) {}
		/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
		virtual void dgTdx_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
		virtual void dgTdp_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
		virtual void dgTdT_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {}


		/** Terminal inequality constraints hT(T,x,p) < 0 */
		virtual void hTfct(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p) {}
		/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
		virtual void dhTdx_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
		virtual void dhTdp_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {}
		/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
		virtual void dhTdT_vec(VectorRef out, ctypeRNum T, cVectorRef x, cVectorRef p, cVectorRef vec) {}


		/** Additional functions required for semi-implicit systems
		M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
		------------------------------------------------------- **/
		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {}
		/** Jacobian df/dx in vector form (column-wise) **/
		virtual void dfdxtrans(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {}
		/** Jacobian df/dt **/
		virtual void dfdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) {}
		/** Jacobian d(dH/dx)/dt  **/
		virtual void dHdxdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef vec, cVectorRef p) {}
		/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mfct(VectorRef out) {}
		/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
		virtual void Mtrans(VectorRef out) {}
};

/*Template class needed for overriting C++ functions in python*/
class PyProblem : public ProblemBase //PyProblem? 
{
    public:
        /* Inherit the constructor*/
        using ProblemBase::ProblemBase;

        void ffct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override 
        {
            PYBIND11_OVERRIDE_PURE(
                void, 
                ProblemBase, 
                ffct,
                out, t, x, u, p);
        }

        void dfdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override
        {
            PYBIND11_OVERRIDE_PURE(
                void, 
                ProblemBase, 
                dfdx_vec,
                out, t, x, vec, u, p);
        }

        void dfdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dfdu_vec,
                out, t, x, vec, u, p);
        }

        void dfdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef vec, cVectorRef u, cVectorRef p) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dfdp_vec,
                out, t, x, vec, u, p);
        }

        void lfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                lfct,
                out, t, x, u, p, xdes, udes);
        }

        void dldx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dldx,
                out, t, x, u, p, xdes, udes);
        }

        void dldu(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dldu,
                out, t, x, u, p, xdes, udes);
        }
        
        void dldp(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef xdes, cVectorRef udes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dldp,
                out, t, x, u, p, xdes, udes);
        }

        void Vfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef xdes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                Vfct,
                out, t, x, p, xdes);
        }

        void dVdx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef xdes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dVdx,
                out, t, x, p, xdes);
        }

        void dVdp(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef xdes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dVdp,
                out, t, x, p, xdes);
        }
        
        void dVdT(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef xdes) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dVdT,
                out, t, x, p, xdes);
        }

		void gfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                gfct,
                out, t, x, u, p);
        }

		void dgdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dgdx_vec,
                out, t, x, u, p, vec);
        }

		void dgdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dgdu_vec,
                out, t, x, u, p, vec);
        }
		void dgdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dgdp_vec,
                out, t, x, u, p, vec);
        }

		void hfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                hfct,
                out, t, x, u, p);
        }

		void dhdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dhdx_vec,
                out, t, x, u, p, vec);
        }

		void dhdu_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dhdu_vec,
                out, t, x, u, p, vec);
        }
		void dhdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p, cVectorRef vec) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dhdp_vec,
                out, t, x, u, p, vec);
        }

		void gTfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                gTfct,
                out, t, x, p);
        }

		void dgTdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dgTdx_vec,
                out, t, x, p, vec);
        }

		void dgTdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dgTdp_vec,
                out, t, x, p, vec);
        }

		void dgTdT_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dgTdT_vec,
                out, t, x, p, vec);
        }

		void hTfct(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                hTfct,
                out, t, x, p);
        }

		void dhTdx_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dhTdx_vec,
                out, t, x, p, vec);
        }

		void dhTdp_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dhTdp_vec,
                out, t, x, p, vec);
        }

		void dhTdT_vec(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef p, cVectorRef vec) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dhTdT_vec,
                out, t, x, p, vec);
        }

		void dfdx(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dfdx,
                out, t, x, u, p);
        }
        
		void dfdxtrans(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dfdxtrans,
                out, t, x, u, p);
        }
        
		void dfdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef p) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dfdt,
                out, t, x, u, p);
        }
        
		void dHdxdt(VectorRef out, ctypeRNum t, cVectorRef x, cVectorRef u, cVectorRef vec, cVectorRef p) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                dHdxdt,
                out, t, x, u, vec, p);
        }
        
		void Mfct(VectorRef out) override
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                Mfct,
                out);
        }
		void Mtrans(VectorRef out) 
        {
            PYBIND11_OVERRIDE(
                void, 
                ProblemBase, 
                Mtrans,
                out);
        }
};

#endif
