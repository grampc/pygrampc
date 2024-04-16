/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#include "grampc_interface.hpp"

namespace grampc
{
    // Eigen::Map needs to be explicitly constructed
    GrampcBinding::grampc_param::grampc_param()
        : x0(NULL, 0), xdes(NULL, 0), u0(NULL, 0), udes(NULL, 0),
        umax(NULL, 0), umin(NULL, 0), p0(NULL, 0), pmax(NULL, 0), pmin(NULL, 0)
    {
    }

    void GrampcBinding::grampc_param::remap_memory(const typeGRAMPC *grampc)
    {
        Nx = &grampc->param->Nx;
        Nu = &grampc->param->Nu;
        Np = &grampc->param->Np;
        Ng = &grampc->param->Ng;
        Nh = &grampc->param->Nh;
        NgT = &grampc->param->NgT;
        NhT = &grampc->param->NhT;
        Nc = &grampc->param->Nc;

        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
        new (&x0) Eigen::Map<Vector> (grampc->param->x0, grampc->param->Nx);
        new (&xdes) Eigen::Map<Vector>(grampc->param->xdes, grampc->param->Nx);

        new (&u0) Eigen::Map<Vector>(grampc->param->u0, grampc->param->Nu);
        new (&udes) Eigen::Map<Vector>(grampc->param->udes, grampc->param->Nu);
        new (&umax) Eigen::Map<Vector>(grampc->param->umax, grampc->param->Nu);
        new (&umin) Eigen::Map<Vector>(grampc->param->umin, grampc->param->Nu);

        new (&p0) Eigen::Map<Vector>(grampc->param->p0, grampc->param->Np);
        new (&pmax) Eigen::Map<Vector>(grampc->param->pmax, grampc->param->Np);
        new (&pmin) Eigen::Map<Vector>(grampc->param->pmin, grampc->param->Np);

        Thor = &grampc->param->Thor;
        Tmax = &grampc->param->Tmax;
        Tmin = &grampc->param->Tmin;

        dt = &grampc->param->dt;
        t0 = &grampc->param->t0;
    }

    // Eigen::Map needs to be explicitly constructed
    GrampcBinding::grampc_opt::grampc_opt()
        : xScale(NULL, 0), xOffset(NULL, 0), uScale(NULL, 0), uOffset(NULL, 0),
        pScale(NULL, 0), pOffset(NULL, 0), cScale(NULL, 0), ConstraintsAbsTol(NULL, 0)
    {
    }

    void GrampcBinding::grampc_opt::remap_memory(const typeGRAMPC *grampc)
    {
        Nhor = &grampc->opt->Nhor;
        MaxGradIter = &grampc->opt->MaxGradIter;
        MaxMultIter = &grampc->opt->MaxMultIter;
        ShiftControl = &grampc->opt->ShiftControl;

        TimeDiscretization = &grampc->opt->TimeDiscretization;

        IntegralCost = &grampc->opt->IntegralCost;
        TerminalCost = &grampc->opt->TerminalCost;
        IntegratorCost = &grampc->opt->IntegratorCost;
        Integrator = &grampc->opt->Integrator;
        IntegratorRelTol = &grampc->opt->IntegratorRelTol;
        IntegratorAbsTol = &grampc->opt->IntegratorAbsTol;
        IntegratorMinStepSize = &grampc->opt->IntegratorMinStepSize;
        IntegratorMaxSteps = &grampc->opt->IntegratorMaxSteps;
        FlagsRodas = std::vector<int>(grampc->opt->FlagsRodas, grampc->opt->FlagsRodas + 8);

        LineSearchType = &grampc->opt->LineSearchType;
        LineSearchExpAutoFallback = &grampc->opt->LineSearchExpAutoFallback;
        LineSearchMax = &grampc->opt->LineSearchMax;
        LineSearchMin = &grampc->opt->LineSearchMin;
        LineSearchInit = &grampc->opt->LineSearchInit;
        LineSearchAdaptAbsTol = &grampc->opt->LineSearchAdaptAbsTol;
        LineSearchAdaptFactor = &grampc->opt->LineSearchAdaptFactor;
        LineSearchIntervalTol = &grampc->opt->LineSearchIntervalTol;
        LineSearchIntervalFactor = &grampc->opt->LineSearchIntervalFactor;

        OptimControl = &grampc->opt->OptimControl;
        OptimParam = &grampc->opt->OptimParam;
        OptimParamLineSearchFactor = &grampc->opt->OptimParamLineSearchFactor;
        OptimTime = &grampc->opt->OptimTime;
        OptimTimeLineSearchFactor = &grampc->opt->OptimTimeLineSearchFactor;

        ScaleProblem = &grampc->opt->ScaleProblem;
        
        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
        new (&xScale) Eigen::Map<Vector>(grampc->opt->xScale, grampc->param->Nx);
        new (&xOffset) Eigen::Map<Vector>(grampc->opt->xOffset, grampc->param->Nx);
        new (&uScale) Eigen::Map<Vector>(grampc->opt->uScale, grampc->param->Nu);
        new (&uOffset) Eigen::Map<Vector>(grampc->opt->uOffset, grampc->param->Nu);
        new (&pScale) Eigen::Map<Vector>(grampc->opt->pScale, grampc->param->Np);
        new (&pOffset) Eigen::Map<Vector>(grampc->opt->pOffset, grampc->param->Np);
        TScale = &grampc->opt->TScale;
        TOffset = &grampc->opt->TOffset;
        JScale = &grampc->opt->JScale;
        new (&cScale) Eigen::Map<Vector>(grampc->opt->cScale, grampc->param->Nc);

        EqualityConstraints = &grampc->opt->EqualityConstraints;
        InequalityConstraints = &grampc->opt->InequalityConstraints;
        TerminalEqualityConstraints = &grampc->opt->TerminalEqualityConstraints;
        TerminalInequalityConstraints = &grampc->opt->TerminalInequalityConstraints;
        ConstraintsHandling = &grampc->opt->ConstraintsHandling;
        new (&ConstraintsAbsTol) Eigen::Map<Vector>(grampc->opt->ConstraintsAbsTol, grampc->param->Nc);

        MultiplierMax = &grampc->opt->MultiplierMax;
        MultiplierDampingFactor = &grampc->opt->MultiplierDampingFactor;
        PenaltyMax = &grampc->opt->PenaltyMax;
        PenaltyMin = &grampc->opt->PenaltyMin;
        PenaltyIncreaseFactor = &grampc->opt->PenaltyIncreaseFactor;
        PenaltyDecreaseFactor = &grampc->opt->PenaltyDecreaseFactor;
        PenaltyIncreaseThreshold = &grampc->opt->PenaltyIncreaseThreshold;
        AugLagUpdateGradientRelTol = &grampc->opt->AugLagUpdateGradientRelTol;

        ConvergenceCheck = &grampc->opt->ConvergenceCheck;
        ConvergenceGradientRelTol = &grampc->opt->ConvergenceGradientRelTol;
    }

    // Eigen::Map needs to be explicitly constructed
    GrampcBinding::grampc_sol::grampc_sol()
        : xnext(NULL, 0), unext(NULL, 0), pnext(NULL, 0), J(NULL, 0)
    {
    }

    void GrampcBinding::grampc_sol::remap_memory(const typeGRAMPC *grampc)
    {
        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
        new (&xnext) Eigen::Map<Vector>(grampc->sol->xnext, grampc->param->Nx);
        new (&unext) Eigen::Map<Vector>(grampc->sol->unext, grampc->param->Nu);
        new (&pnext) Eigen::Map<Vector>(grampc->sol->pnext, grampc->param->Np);
        Tnext = &grampc->sol->Tnext;
        new (&J) Eigen::Map<Vector>((typeRNum*)&grampc->sol->J[0], 2);
        cfct = &grampc->sol->cfct;
        pen = &grampc->sol->pen;
        iter = std::vector<int>(grampc->sol->iter, (grampc->sol->iter) + (grampc->opt->MaxMultIter));
        status = grampc->sol->status;
    }

    // Eigen::Map needs to be explicitly constructed
    GrampcBinding::grampc_rws::grampc_rws()
        : t(NULL, 0), tls(NULL, 0), x(NULL, 0, 0), adj(NULL, 0, 0), dcdx(NULL, 0, 0), u(NULL, 0, 0), uls(NULL, 0, 0), 
        uprev(NULL, 0, 0), gradu(NULL, 0, 0), graduprev(NULL, 0, 0), dcdu(NULL, 0, 0), p(NULL, 0), pls(NULL, 0), 
        pprev(NULL, 0), gradp(NULL, 0), gradpprev(NULL, 0), dcdp(NULL, 0, 0), mult(NULL, 0, 0), pen(NULL, 0, 0), 
        cfct(NULL, 0, 0), cfctprev(NULL, 0, 0), cfctAbsTol(NULL, 0), lsAdapt(NULL, 0), lsExplicit(NULL, 0), 
        rwsScale(NULL, 0), rwsGeneral(NULL, 0), rparRodas(NULL, 0), workRodas(NULL, 0)
    {
    }

    void GrampcBinding::grampc_rws::remap_memory(const typeGRAMPC *grampc)
    {
        // placement new uses the preallocated memory of the Eigen::Map types on the stack, so no delete has to be called
        // this is the way to go according to the Eigen documentation: https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html

        // dim = Nhor;
        new (&t) Eigen::Map<Vector>(grampc->rws->t, grampc->opt->Nhor);
        new (&tls) Eigen::Map<Vector>(grampc->rws->tls, grampc->opt->Nhor);

        // dim = Nx x Nhor;
        new (&x) Eigen::Map<Matrix>(grampc->rws->x, grampc->param->Nx, grampc->opt->Nhor);
        new (&adj) Eigen::Map<Matrix>(grampc->rws->adj, grampc->param->Nx, grampc->opt->Nhor);
        new (&dcdx) Eigen::Map<Matrix>(grampc->rws->dcdx, grampc->param->Nx, grampc->opt->Nhor + 1);

        // dim = Nu x Nhor;
        new (&u) Eigen::Map<Matrix>(grampc->rws->u, grampc->param->Nu, grampc->opt->Nhor);
        new (&uls) Eigen::Map<Matrix>(grampc->rws->uls, grampc->param->Nu, grampc->opt->Nhor);
        new (&uprev) Eigen::Map<Matrix>(grampc->rws->uprev, grampc->param->Nu, grampc->opt->Nhor);
        new (&gradu) Eigen::Map<Matrix>(grampc->rws->gradu, grampc->param->Nu, grampc->opt->Nhor);
        new (&graduprev) Eigen::Map<Matrix>(grampc->rws->graduprev, grampc->param->Nu, grampc->opt->Nhor);
        new (&dcdu) Eigen::Map<Matrix>(grampc->rws->dcdu, grampc->param->Nu, grampc->opt->Nhor);

        // dim = Np
        new (&p) Eigen::Map<Vector>(grampc->rws->p, grampc->param->Np);
        new (&pls) Eigen::Map<Vector>(grampc->rws->pls, grampc->param->Np);
        new (&pprev) Eigen::Map<Vector>(grampc->rws->pprev, grampc->param->Np);
        new (&gradp) Eigen::Map<Vector>(grampc->rws->gradp, grampc->param->Np);
        new (&gradpprev) Eigen::Map<Vector>(grampc->rws->gradpprev, grampc->param->Np);
        // dim = Np x (Nhor + 1)
        new (&dcdp) Eigen::Map<Matrix>(grampc->rws->dcdp, grampc->param->Np, grampc->opt->Nhor + 1);

        // dim = 1
        T = &grampc->rws->T;
        Tprev = &grampc->rws->Tprev;
        gradT = &grampc->rws->gradT;
        gradTprev = &grampc->rws->gradTprev;
        dcdt = &grampc->rws->dcdt;

        // dim = Nc x Nhor
        new (&mult) Eigen::Map<Matrix>(grampc->rws->mult, grampc->param->Nc, grampc->opt->Nhor);
        new (&pen) Eigen::Map<Matrix>(grampc->rws->pen, grampc->param->Nc, grampc->opt->Nhor);
        new (&cfct) Eigen::Map<Matrix>(grampc->rws->cfct, grampc->param->Nc, grampc->opt->Nhor);
        new (&cfctprev) Eigen::Map<Matrix>(grampc->rws->cfctprev, grampc->param->Nc, grampc->opt->Nhor);

        // dim = 1 x Nc
        new (&cfctAbsTol) Eigen::Map<Vector>(grampc->rws->cfctAbsTol, grampc->param->Nc);

        // dim = DYN!
        // Only memory for the selected linesearch is allocated
        if (grampc->opt->LineSearchType == INT_ADAPTIVELS) 
        {
            const int dim = 2 * (NALS + 1) * (1 + grampc->opt->MaxGradIter);
            new (&lsAdapt) Eigen::Map<Vector>(grampc->rws->lsAdapt, dim);
        }
        else // if(grampc->opt->LineSearchType == INT_EXPLICIT)
        {
            new (&lsExplicit) Eigen::Map<Vector>(grampc->rws->lsExplicit, NELS);
        }

        // dim = 2*(Nx+Nu+Np)
        const int dim = 2 * (grampc->param->Nx + grampc->param->Nu + grampc->param->Np);
        new (&rwsScale) Eigen::Map<Vector>(grampc->rws->rwsScale, dim);
        lrwsGeneral = &grampc->rws->lrwsGeneral;
        // dim = lrwsGeneral
        new (&rwsGeneral) Eigen::Map<Vector>(grampc->rws->rwsGeneral, grampc->rws->lrwsGeneral);

        lworkRodas = &grampc->rws->lworkRodas;
        liworkRodas = &grampc->rws->liworkRodas;

        if (grampc->opt->Integrator == INT_RODAS)
        {
            // dim = Nhor;
            new (&rparRodas) Eigen::Map<Vector>(grampc->rws->rparRodas, grampc->opt->Nhor);
            // dim = 20;
            iparRodas = std::vector<int>(grampc->rws->iparRodas, (grampc->rws->iparRodas)+20);
            // dim = lworkRodas
            new (&workRodas) Eigen::Map<Vector>(grampc->rws->workRodas, grampc->rws->lworkRodas);
            // dim = liworkRodas
            iworkRodas = std::vector<int>(grampc->rws->iworkRodas, (grampc->rws->iworkRodas) + (grampc->rws->liworkRodas));
        }
    }

    GrampcBinding::GrampcBinding(ProblemDescriptionPtr problem)
        : problem_description(problem)
    {
        grampc_init(&grampc_, problem_description.get());
        param.remap_memory(grampc_);
        opt.remap_memory(grampc_);
        rws.remap_memory(grampc_);
        sol.remap_memory(grampc_);
    }

    GrampcBinding::~GrampcBinding()
    {
        grampc_free(&grampc_);
    }

    typeRNum GrampcBinding::run()
    {
        // check if dt and Thor are valid. This redundancy is necessary to throw an error for Python
        if (*param.dt <= 0.0)
        {
            throw std::runtime_error("Sampling time dt is not valid. Must be greater than zero");
        } 
        else if (*param.Thor < *param.dt)
        {
            throw std::runtime_error("Horizon Thor is not valid. Must be greater than sampling time");
        }

        auto begin = std::chrono::high_resolution_clock::now();
        grampc_run(grampc_);
        auto end = std::chrono::high_resolution_clock::now();
        // time in milliseconds
        return (typeRNum)std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() * 1e-6; 
    }

    void GrampcBinding::estim_penmin(bool run_grampc)
    {
        // check if dt and Thor are valid. This redundancy is necessary to throw an error for Python
        if (*param.dt <= 0.0 && run_grampc)
        {
            throw std::runtime_error("Sampling time dt is not valid. Must be greater than zero");
        } 
        else if (*param.Thor < *param.dt && run_grampc)
        {
            throw std::runtime_error("Horizon Thor is not valid. Must be greater than sampling time");
        }

        grampc_estim_penmin(grampc_, run_grampc);
    }

    Vector GrampcBinding::ffct(typeRNum t, Vector x, Vector u, Vector p)
    {
        Vector out(problem_description->Nx);
        problem_description->ffct(out, t, x, u, p);
        return out;
    }

    typeRNum GrampcBinding::lfct(typeRNum t, Vector x, Vector u, Vector p, Vector xdes, Vector udes)
    {
        Vector out(1);
        problem_description->lfct(out, t, x, u, p, xdes, udes);
        return out[0];
    }

    typeRNum GrampcBinding::Vfct(typeRNum t, Vector x, Vector p, Vector xdes)
    {
        Vector out(1);
        problem_description->Vfct(out, t, x, p, xdes);
        return out[0];
    }

    Vector GrampcBinding::gfct(typeRNum t, Vector x, Vector u, Vector p)
    {
        Vector out(problem_description->Ng);
        problem_description->gfct(out, t, x, u, p);
        return out;
    }

    Vector GrampcBinding::hfct(typeRNum t, Vector x, Vector u, Vector p)
    {
        Vector out(problem_description->Nh);
        problem_description->hfct(out, t, x, u, p);
        return out;
    }

    Vector GrampcBinding::gTfct(typeRNum t, Vector x, Vector p)
    {
        Vector out(problem_description->NgT);
        problem_description->gTfct(out, t, x, p);
        return out;
    }

    Vector GrampcBinding::hTfct(typeRNum t, Vector x, Vector p)
    {
        Vector out(problem_description->NhT);
        problem_description->hTfct(out, t, x, p);
        return out;
    }

    void GrampcBinding::set_param_real(const std::string &key, typeRNum value)
    {
        grampc_setparam_real(grampc_, key.c_str(), value);
    }

    void GrampcBinding::set_param_real_vec(const std::string &key, const std::vector<typeRNum> &values)
    {
        grampc_setparam_real_vector(grampc_, key.c_str(), values.data());
    }

    void GrampcBinding::set_opt_str(const std::string &key, const std::string &vstr)
    {
        grampc_setopt_string(grampc_, key.c_str(), vstr.c_str());
        if (key == "Integrator" || key == "LineSearchType" || key == "IntegratorCost")
        {    
            rws.remap_memory(grampc_);
        }
    }

    void GrampcBinding::set_opt_int(const std::string &key, int value)
    {
        grampc_setopt_int(grampc_, key.c_str(), value);
        if (key == "Nhor" || key == "MaxMultIter" || key == "MaxGradIter")
        {
            rws.remap_memory(grampc_);
            sol.remap_memory(grampc_);
        }
    }

    void GrampcBinding::set_opt_int_vec(const std::string &key, const std::vector<int> &values)
    {
        grampc_setopt_int_vector(grampc_, key.c_str(), values.data());
        if (values[4] < problem_description->Nx) // Check for MLJAC, if a banded structure is used
        {
            problem_description->Rodas_Jac = problem_description->Nx * (values[4] + values[5] + 1); // Nx * (MLJAC + MUJAC + 1)
        }
        else if (values[4] == problem_description->Nx)
        {
            problem_description->Rodas_Jac = problem_description->Nx * problem_description->Nx;
        }

        if (values[6] < problem_description->Nx) // Check for MLMAS, if a banded structure is used
        {
            problem_description->Rodas_M = problem_description->Nx * (values[6] + values[7] + 1); // Nx * (MLMAS + MUMAS + 1)
        }
        else if (values[6] == problem_description->Nx)
        {
            problem_description->Rodas_M = problem_description->Nx * problem_description->Nx;
        }
    }

    void GrampcBinding::set_opt_real(const std::string &key, typeRNum value)
    {
        grampc_setopt_real(grampc_, key.c_str(), value);
    }

    void GrampcBinding::set_opt_real_vec(const std::string &key, const std::vector<typeRNum> &values)
    {
        grampc_setopt_real_vector(grampc_, key.c_str(), values.data());
    }

    void GrampcBinding::fill_rws_memory(Eigen::Ref<Matrix> rws_matrix, const Eigen::Ref<const Matrix>& new_data)
    {
        if (rws_matrix.data() == nullptr)
        {
            throw std::runtime_error("Memory for matrix is not allocated");
        }
        if (rws_matrix.rows() != new_data.rows() || rws_matrix.cols() != new_data.cols())
        {
            // Formats an error message which outputs the expected dimensions and the dimensions of 
            // new_data. Results in an message like: 
            // "Wrong dimensions detected. Expected (2, 10), got (2, 15)."
            std::ostringstream error_message; 
            error_message << "Wrong dimensions detected. Expected (" << rws_matrix.rows() << ", " << rws_matrix.cols() << "), got ("
                << new_data.rows() << ", " << new_data.cols() << ")";
            throw std::length_error(error_message.str());
        }

        // copy the new data to the rws_matrix.
        rws_matrix = new_data;
    }

    void GrampcBinding::set_rws_u(const Eigen::Ref<const Matrix>& u_new)
    {
        fill_rws_memory(rws.u, u_new);
    }

    void GrampcBinding::set_rws_multiplier(const Eigen::Ref<const Matrix>& multiplier_new)
    {
        fill_rws_memory(rws.mult, multiplier_new);
    }

    void GrampcBinding::set_rws_penalty(const Eigen::Ref<const Matrix>& penalty_new)
    {
        fill_rws_memory(rws.pen, penalty_new);
    }

    void GrampcBinding::print_params()
    {
        grampc_printparam(grampc_);
    }

    void GrampcBinding::print_opts()
    {
        grampc_printopt(grampc_);
    }

    void GrampcBinding::print_status()
    {
        grampc_printstatus(sol.status, STATUS_LEVEL_DEBUG);
    }
}