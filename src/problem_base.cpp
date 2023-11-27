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

extern "C"
{

	void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        *Nx = problem->Nx_;
        *Nu = problem->Nu_;
        *Np = problem->Np_;
        *Ng = problem->Ng_;
        *Nh = problem->Nh_;
        *NgT = problem->NgT_;
        *NhT = problem->NhT_;
	}

	void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        // mapping the C pointer to an Eigen::Vector, so that the data is exposed as an numpy array in Python
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);

		problem->ffct(outMap, t, xMap, uMap, pMap);
	}
    
	void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dfdx_vec(outMap, t, xMap, vecMap, uMap, pMap);
	}

	void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nu_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dfdu_vec(outMap, t, xMap, vecMap, uMap, pMap);
	}

	void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Np_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dfdp_vec(outMap, t, xMap, vecMap, uMap, pMap);
	}

	void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, 1);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> udesMap(udes, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->lfct(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
	}

	void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> udesMap(udes, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dldx(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
	}

	void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nu_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> udesMap(udes, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dldu(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
	}

	void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Np_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> udesMap(udes, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dldp(outMap, t, xMap, uMap, pMap, xdesMap, udesMap);
	}

	void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, 1);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->Vfct(outMap, T, xMap, pMap, xdesMap);
	}

	void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dVdx(outMap, T, xMap, pMap, xdesMap);
	}

	void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Np_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dVdp(outMap, T, xMap, pMap, xdesMap);
	}

	void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, 1);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> xdesMap(xdes, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dVdT(outMap, T, xMap, pMap, xdesMap);
	}

	void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Ng_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->gfct(outMap, t, xMap, uMap, pMap);
	}

	void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dgdx_vec(outMap, t, xMap, uMap, pMap, vecMap);
	}

	void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nu_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dgdu_vec(outMap, t, xMap, uMap, pMap, vecMap);
	}

	void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Np_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dgdp_vec(outMap, t, xMap, uMap, pMap, vecMap);
	}

	void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nh_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->hfct(outMap, t, xMap, uMap, pMap);
	}

	void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dhdx_vec(outMap, t, xMap, uMap, pMap, vecMap);
	}

	void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nu_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dhdu_vec(outMap, t, xMap, uMap, pMap, vecMap);
	}

	void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Np_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dhdp_vec(outMap, t, xMap, uMap, pMap, vecMap);
	}

	void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->NgT_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->gTfct(outMap, T, xMap, pMap);
	}

	void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dgTdx_vec(outMap, T, xMap, pMap, vecMap);
	}

	void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Np_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dgTdp_vec(outMap, T, xMap, pMap, vecMap);
	}

	void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, 1);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dgTdT_vec(outMap, T, xMap, pMap, vecMap);
	}

	void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->NhT_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->hTfct(outMap, T, xMap, pMap);
	}

	void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dhTdx_vec(outMap, T, xMap, pMap, vecMap);
	}

	void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Np_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);

		problem->dhTdp_vec(outMap, T, xMap, pMap, vecMap);
	}

	void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM* userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, 1);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dhTdT_vec(outMap, T, xMap, pMap, vecMap);
	}

	void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Rodas_Jac); //Nx * (MLJAC + MUJAC + 1)
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dfdx(outMap, t, xMap, uMap, pMap);
	}

	void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Rodas_Jac); //Nx * (MLJAC + MUJAC + 1)
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dfdxtrans(outMap, t, xMap, uMap, pMap);
	}

	void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dfdt(outMap, t, xMap, uMap, pMap);
	}

	void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec, ctypeRNum *p, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Nx_);
        Eigen::Map<const Vector> xMap(x, problem->Nx_);
        Eigen::Map<const Vector> vecMap(vec, problem->Nx_);
        Eigen::Map<const Vector> uMap(u, problem->Nu_);
        Eigen::Map<const Vector> pMap(p, problem->Np_);
        
		problem->dHdxdt(outMap, t, xMap, uMap, vecMap, pMap);
	}

	void Mfct(typeRNum *out, typeUSERPARAM *userparam) // Auf Python Seite eine Matrix?
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Rodas_M); //Nx * (MLMAS + MUMAS + 1)
        
		problem->Mfct(outMap);
	}

	void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
	{
        ProblemBase* problem = (ProblemBase*) userparam;
        Eigen::Map<Vector> outMap(out, problem->Rodas_M); //Nx * (MLMAS + MUMAS + 1)
        
		problem->Mtrans(outMap);
	}

}
