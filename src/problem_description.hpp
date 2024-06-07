/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#ifndef PROBLEM_DESCRIPTION_HPP
#define PROBLEM_DESCRIPTION_HPP

extern "C"
{
    #include "grampc.h"
}

#include <Eigen/Dense>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace grampc
{
    typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Ref<Vector> VectorRef;
    typedef const Eigen::Ref<const Vector>& VectorConstRef;

    class PYBIND11_EXPORT ProblemDescription
    {
        public:
            typeInt Nx;
            typeInt Nu;
            typeInt Np;
            typeInt Ng;
            typeInt Nh;
            typeInt NgT;
            typeInt NhT;
            typeInt Rodas_Jac = 0;
            typeInt Rodas_M = 0;

        public:

            virtual ~ProblemDescription() {}

            /** System function f(t,x,u,p,userparam)
            ------------------------------------ **/
            virtual void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) = 0;
            /** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
            virtual void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) = 0;
            /** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
            virtual void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
            virtual void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) {}


            /** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam)
            -------------------------------------------------- **/
            virtual void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}
            /** Gradient dl/dx **/
            virtual void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}
            /** Gradient dl/du **/
            virtual void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}
            /** Gradient dl/dp **/
            virtual void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) {}


            /** Terminal cost V(T,x,p) */
            virtual void Vfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}
            /** Gradient dV/dx **/
            virtual void dVdx(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}
            /** Gradient dV/dp **/
            virtual void dVdp(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}
            /** Gradient dV/dT **/
            virtual void dVdT(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) {}


            /** Equality constraints g(t,x,u,p) = 0 */
            virtual void gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
            virtual void dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
            virtual void dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
            virtual void dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}


            /** Inequality constraints h(t,x,u,p) < 0 */
            virtual void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dh/dx) **/
            virtual void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dh/du) **/
            virtual void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dh/dp) **/
            virtual void dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) {}


            /** Terminal equality constraints gT(T,x,p) = 0 */
            virtual void gTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p) {}
            /** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
            virtual void dgTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
            virtual void dgTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
            virtual void dgTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}


            /** Terminal inequality constraints hT(T,x,p) < 0 */
            virtual void hTfct(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p) {}
            /** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
            virtual void dhTdx_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
            virtual void dhTdp_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}
            /** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
            virtual void dhTdT_vec(VectorRef out, ctypeRNum T, VectorConstRef x, VectorConstRef p, VectorConstRef vec) {}


            /** Additional functions required for semi-implicit systems
            M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS
            ------------------------------------------------------- **/
            /** Jacobian df/dx in vector form (column-wise) **/
            virtual void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian df/dx in vector form (column-wise) **/
            virtual void dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian df/dt **/
            virtual void dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian d(dH/dx)/dt  **/
            virtual void dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef vec, VectorConstRef p) {}
            /** Mass matrix in vector form (column-wise, either banded or full matrix) **/
            virtual void Mfct(VectorRef out) {}
            /** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
            virtual void Mtrans(VectorRef out) {}
            

            /** Additional functions required for Taylor-SMPC in GRAMPC-S */
            /*------------------------------------------------*/

            /** Jacobian df/dp in vector form  **/
            virtual void dfdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Hessian d(df/dx)/dx in vector form  */
            virtual void dfdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian d(df/dx)/dp in vector form  */
            virtual void dfdxdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Hessian d(df/dx)/du in vector form  */
            virtual void dfdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian d(df/dp)/du in vector form  */
            virtual void dfdpdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            
            /** Jacobian dh/dx in vector form  **/
            virtual void dhdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian dh/du in vector form  **/
            virtual void dhdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Hessian d(dh/dx)/dx in vector form  **/
            virtual void dhdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian d(dh/dx)/du in vector form **/
            virtual void dhdxdu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) {}
            /** Jacobian dhT/dx in vector form  **/
            virtual void dhTdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}
            /** Hessian d(dhT/dx)/dx in vector form  **/
            virtual void dhTdxdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}
            /** Hessian d(dhT/dx)/dT in vector form  **/
            virtual void dhTdxdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) {}
    };
    
    // Alias
    typedef std::shared_ptr<ProblemDescription> ProblemDescriptionPtr;
    typedef std::shared_ptr<const ProblemDescription> ProblemDescriptionConstPtr;

    /*Template class needed for overriting C++ functions in python*/
    class PyProblem : public ProblemDescription
    {
        public:
            /* Inherit the constructor*/
            using ProblemDescription::ProblemDescription;

            void ffct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override 
            {
                PYBIND11_OVERRIDE_PURE(
                    void, 
                    ProblemDescription, 
                    ffct,
                    out, t, x, u, p);
            }

            void dfdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override
            {
                PYBIND11_OVERRIDE_PURE(
                    void, 
                    ProblemDescription, 
                    dfdx_vec,
                    out, t, x, vec, u, p);
            }

            void dfdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdu_vec,
                    out, t, x, vec, u, p);
            }

            void dfdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef vec, VectorConstRef u, VectorConstRef p) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdp_vec,
                    out, t, x, vec, u, p);
            }

            void lfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    lfct,
                    out, t, x, u, p, xdes, udes);
            }

            void dldx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dldx,
                    out, t, x, u, p, xdes, udes);
            }

            void dldu(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dldu,
                    out, t, x, u, p, xdes, udes);
            }
            
            void dldp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef xdes, VectorConstRef udes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dldp,
                    out, t, x, u, p, xdes, udes);
            }

            void Vfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    Vfct,
                    out, t, x, p, xdes);
            }

            void dVdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dVdx,
                    out, t, x, p, xdes);
            }

            void dVdp(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dVdp,
                    out, t, x, p, xdes);
            }
            
            void dVdT(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef xdes) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dVdT,
                    out, t, x, p, xdes);
            }

            void gfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    gfct,
                    out, t, x, u, p);
            }

            void dgdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgdx_vec,
                    out, t, x, u, p, vec);
            }

            void dgdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgdu_vec,
                    out, t, x, u, p, vec);
            }
            void dgdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgdp_vec,
                    out, t, x, u, p, vec);
            }

            void hfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    hfct,
                    out, t, x, u, p);
            }

            void dhdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhdx_vec,
                    out, t, x, u, p, vec);
            }

            void dhdu_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhdu_vec,
                    out, t, x, u, p, vec);
            }
            void dhdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p, VectorConstRef vec) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhdp_vec,
                    out, t, x, u, p, vec);
            }

            void gTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    gTfct,
                    out, t, x, p);
            }

            void dgTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgTdx_vec,
                    out, t, x, p, vec);
            }

            void dgTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgTdp_vec,
                    out, t, x, p, vec);
            }

            void dgTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dgTdT_vec,
                    out, t, x, p, vec);
            }

            void hTfct(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    hTfct,
                    out, t, x, p);
            }

            void dhTdx_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhTdx_vec,
                    out, t, x, p, vec);
            }

            void dhTdp_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhTdp_vec,
                    out, t, x, p, vec);
            }

            void dhTdT_vec(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef p, VectorConstRef vec) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dhTdT_vec,
                    out, t, x, p, vec);
            }

            void dfdx(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdx,
                    out, t, x, u, p);
            }
            
            void dfdxtrans(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdxtrans,
                    out, t, x, u, p);
            }
            
            void dfdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef p) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dfdt,
                    out, t, x, u, p);
            }
            
            void dHdxdt(VectorRef out, ctypeRNum t, VectorConstRef x, VectorConstRef u, VectorConstRef vec, VectorConstRef p) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    dHdxdt,
                    out, t, x, u, vec, p);
            }
            
            void Mfct(VectorRef out) override
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    Mfct,
                    out);
            }
            void Mtrans(VectorRef out) 
            {
                PYBIND11_OVERRIDE(
                    void, 
                    ProblemDescription, 
                    Mtrans,
                    out);
            }
    };
}

#endif
