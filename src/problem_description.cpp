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

extern "C"
{
    void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        *Nx = problem->Nx;
        *Nu = problem->Nu;
        *Np = problem->Np;
        *Ng = problem->Ng;
        *Nh = problem->Nh;
        *NgT = problem->NgT;
        *NhT = problem->NhT;
    }

    void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);

        problem->ffct(outMap, t, xMap, uMap, pMap);
    }
    
    void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dfdx_vec(outMap, t, xMap, vecMap, uMap, pMap);
    }

    void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dfdu_vec(outMap, t, xMap, vecMap, uMap, pMap);
    }

    void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dfdp_vec(outMap, t, xMap, vecMap, uMap, pMap);
    }

    void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> udesMap(udes, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->lfct(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
    }

    void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> udesMap(udes, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dldx(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
    }

    void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> udesMap(udes, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dldu(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
    }

    void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> udesMap(udes, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dldp(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
    }

    void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->Vfct(outMap, T, xMap, pMap, xdesMap);
    }

    void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dVdx(outMap, T, xMap, pMap, xdesMap);
    }

    void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dVdp(outMap, T, xMap, pMap, xdesMap);
    }

    void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> xdesMap(xdes, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dVdT(outMap, T, xMap, pMap, xdesMap);
    }

    void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Ng);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->gfct(outMap, t, xMap, uMap, pMap);
    }

    void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Ng);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dgdx_vec(outMap, t, xMap, uMap, pMap, vecMap);
    }

    void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Ng);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dgdu_vec(outMap, t, xMap, uMap, pMap, vecMap);
    }

    void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Ng);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dgdp_vec(outMap, t, xMap, uMap, pMap, vecMap);
    }

    void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nh);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->hfct(outMap, t, xMap, uMap, pMap);
    }

    void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nh);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dhdx_vec(outMap, t, xMap, uMap, pMap, vecMap);
    }

    void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nu);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nh);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dhdu_vec(outMap, t, xMap, uMap, pMap, vecMap);
    }

    void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nh);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dhdp_vec(outMap, t, xMap, uMap, pMap, vecMap);
    }

    void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->NgT);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->gTfct(outMap, T, xMap, pMap);
    }

    void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NgT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dgTdx_vec(outMap, T, xMap, pMap, vecMap);
    }

    void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NgT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dgTdp_vec(outMap, T, xMap, pMap, vecMap);
    }

    void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NgT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dgTdT_vec(outMap, T, xMap, pMap, vecMap);
    }

    void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->NhT);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->hTfct(outMap, T, xMap, pMap);
    }

    void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NhT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dhTdx_vec(outMap, T, xMap, pMap, vecMap);
    }

    void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Np);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NhT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);

        problem->dhTdp_vec(outMap, T, xMap, pMap, vecMap);
    }

    void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, 1);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->NhT);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dhTdT_vec(outMap, T, xMap, pMap, vecMap);
    }

    void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_Jac); //Nx * (MLJAC + MUJAC + 1)
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dfdx(outMap, t, xMap, uMap, pMap);
    }

    void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_Jac); //Nx * (MLJAC + MUJAC + 1)
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dfdxtrans(outMap, t, xMap, uMap, pMap);
    }

    void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dfdt(outMap, t, xMap, uMap, pMap);
    }

    void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec, ctypeRNum *p, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Nx);
        Eigen::Map<const grampc::Vector> xMap(x, problem->Nx);
        Eigen::Map<const grampc::Vector> vecMap(vec, problem->Nx);
        Eigen::Map<const grampc::Vector> uMap(u, problem->Nu);
        Eigen::Map<const grampc::Vector> pMap(p, problem->Np);
        
        problem->dHdxdt(outMap, t, xMap, uMap, vecMap, pMap);
    }

    void Mfct(typeRNum *out, typeUSERPARAM *userparam) // Auf Python Seite eine Matrix?
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_M); //Nx * (MLMAS + MUMAS + 1)
        
        problem->Mfct(outMap);
    }

    void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
    {
        grampc::ProblemDescription* problem = (grampc::ProblemDescription*) userparam;
        Eigen::Map<grampc::Vector> outMap(out, problem->Rodas_M); //Nx * (MLMAS + MUMAS + 1)
        
        problem->Mtrans(outMap);
    }

}
