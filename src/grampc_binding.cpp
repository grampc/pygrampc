/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#include "grampc_binding.hpp"

// Eigen::Map needs to be explicitly constructed
GrampcBinding::grampc_param::grampc_param()
    : x0(NULL, 0), xdes(NULL, 0), u0(NULL, 0), udes(NULL, 0),
      umax(NULL, 0), umin(NULL, 0), p0(NULL, 0), pmax(NULL, 0), pmin(NULL, 0)
{
}

void GrampcBinding::grampc_param::reMapMemory(const typeGRAMPC *grampc)
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

void GrampcBinding::grampc_opt::reMapMemory(const typeGRAMPC *grampc)
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

void GrampcBinding::grampc_sol::reMapMemory(const typeGRAMPC *grampc)
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
      uprev(NULL, 0, 0), gradu(NULL, 0, 0), graduprev(NULL, 0, 0), dcdu(NULL, 0, 0), p(NULL, 0, 0), pls(NULL, 0, 0), 
      pprev(NULL, 0, 0), gradp(NULL, 0, 0), gradpprev(NULL, 0, 0), dcdp(NULL, 0, 0), mult(NULL, 0, 0), pen(NULL, 0, 0), 
      cfct(NULL, 0, 0), cfctprev(NULL, 0, 0), cfctAbsTol(NULL, 0), lsAdapt(NULL, 0), lsExplicit(NULL, 0), 
      rwsScale(NULL, 0), rwsGeneral(NULL, 0), rparRodas(NULL, 0), workRodas(NULL, 0)
{
}

void GrampcBinding::grampc_rws::reMapMemory(const typeGRAMPC *grampc)
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

    // dim = Np x Nhor
    new (&p) Eigen::Map<Matrix>(grampc->rws->p, grampc->param->Np, grampc->opt->Nhor);
    new (&pls) Eigen::Map<Matrix>(grampc->rws->pls, grampc->param->Np, grampc->opt->Nhor);
    new (&pprev) Eigen::Map<Matrix>(grampc->rws->pprev, grampc->param->Np, grampc->opt->Nhor);
    new (&gradp) Eigen::Map<Matrix>(grampc->rws->gradp, grampc->param->Np, grampc->opt->Nhor);
    new (&gradpprev) Eigen::Map<Matrix>(grampc->rws->gradpprev, grampc->param->Np, grampc->opt->Nhor);
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

GrampcBinding::GrampcBinding(ProblemBase *problem)
    : problem_binding(problem)
{
    grampc_init(&grampc_, problem_binding);
    param.reMapMemory(grampc_);
    opt.reMapMemory(grampc_);
    rws.reMapMemory(grampc_);
    sol.reMapMemory(grampc_);
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
    Vector out(problem_binding->Nx_);
    problem_binding->ffct(out, t, x, u, p);
    return out;
}

typeRNum GrampcBinding::lfct(typeRNum t, Vector x, Vector u, Vector p, Vector xdes, Vector udes)
{
    Vector out(1);
    problem_binding->lfct(out, t, x, u, p, xdes, udes);
    return out[0];
}

typeRNum GrampcBinding::Vfct(typeRNum t, Vector x, Vector p, Vector xdes)
{
    Vector out(1);
    problem_binding->Vfct(out, t, x, p, xdes);
    return out[0];
}

Vector GrampcBinding::gfct(typeRNum t, Vector x, Vector u, Vector p)
{
    Vector out(problem_binding->Ng_);
    problem_binding->gfct(out, t, x, u, p);
    return out;
}

Vector GrampcBinding::hfct(typeRNum t, Vector x, Vector u, Vector p)
{
    Vector out(problem_binding->Nh_);
    problem_binding->hfct(out, t, x, u, p);
    return out;
}

Vector GrampcBinding::gTfct(typeRNum t, Vector x, Vector p)
{
    Vector out(problem_binding->NgT_);
    problem_binding->gTfct(out, t, x, p);
    return out;
}

Vector GrampcBinding::hTfct(typeRNum t, Vector x, Vector p)
{
    Vector out(problem_binding->NhT_);
    problem_binding->hTfct(out, t, x, p);
    return out;
}

void GrampcBinding::set_param_real(const std::string &key, typeRNum value)
{
    grampc_setparam_real(grampc_, key.c_str(), value);
}

void GrampcBinding::set_param_real_vec(const std::string &key,
                               const std::vector<typeRNum> &values)
{
    grampc_setparam_real_vector(grampc_, key.c_str(), values.data());
}

void GrampcBinding::set_opt_str(const std::string &key, const std::string &vstr)
{
    grampc_setopt_string(grampc_, key.c_str(), vstr.c_str());
    if (key == "Integrator" || key == "LineSearchType" || key == "IntegratorCost")
    {    
        rws.reMapMemory(grampc_);
    }
}

void GrampcBinding::set_opt_int(const std::string &key, int value)
{
    grampc_setopt_int(grampc_, key.c_str(), value);
    if (key == "Nhor" || key == "MaxMultIter" || key == "MaxGradIter")
    {
        rws.reMapMemory(grampc_);
        sol.reMapMemory(grampc_);
    }
}

void GrampcBinding::set_opt_int_vec(const std::string &key, const std::vector<int> &values)
{
    grampc_setopt_int_vector(grampc_, key.c_str(), values.data());
    if (values[4] < problem_binding->Nx_) // Check for MLJAC, if a banded structure is used
    {
        problem_binding->Rodas_Jac = problem_binding->Nx_ * (values[4] + values[5] + 1); // Nx * (MLJAC + MUJAC + 1)
    }
    else if (values[4] == problem_binding->Nx_)
    {
        problem_binding->Rodas_Jac = problem_binding->Nx_ * problem_binding->Nx_;
    }

    if (values[6] < problem_binding->Nx_) // Check for MLMAS, if a banded structure is used
    {
        problem_binding->Rodas_M = problem_binding->Nx_ * (values[6] + values[7] + 1); // Nx * (MLMAS + MUMAS + 1)
    }
    else if (values[6] == problem_binding->Nx_)
    {
        problem_binding->Rodas_M = problem_binding->Nx_ * problem_binding->Nx_;
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

PYBIND11_MODULE(_core, m)
{
    pybind11::class_<GrampcBinding> binding(m, "GrampcBinding");
    binding.def(pybind11::init<ProblemBase *>())
        .def_readonly("param", &GrampcBinding::param)
        .def_readonly("opt", &GrampcBinding::opt)
        .def_readonly("sol", &GrampcBinding::sol)
        .def_readonly("rws", &GrampcBinding::rws)
        .def_readwrite("problem", &GrampcBinding::problem_binding)
        .def("run", &GrampcBinding::run)
        .def("_set_param_real", &GrampcBinding::set_param_real)
        .def("_set_param_real_vec", &GrampcBinding::set_param_real_vec)
        .def("_set_opt_str", &GrampcBinding::set_opt_str)
        .def("_set_opt_int", &GrampcBinding::set_opt_int)
        .def("_set_opt_int_vec", &GrampcBinding::set_opt_int_vec)
        .def("_set_opt_real", &GrampcBinding::set_opt_real)
        .def("_set_opt_real_vec", &GrampcBinding::set_opt_real_vec)
        .def("print_opts", &GrampcBinding::print_opts)
        .def("print_params", &GrampcBinding::print_params)
        .def("print_status", &GrampcBinding::print_status)

        .def("ffct", &GrampcBinding::ffct)
        .def("lfct", &GrampcBinding::lfct)
        .def("Vfct", &GrampcBinding::Vfct)
        .def("gfct", &GrampcBinding::gfct)
        .def("hfct", &GrampcBinding::hfct)
        .def("gTfct", &GrampcBinding::gTfct)
        .def("hTfct", &GrampcBinding::hTfct)
        .def("estim_penmin", &GrampcBinding::estim_penmin);

    typedef GrampcBinding::grampc_param prefix_param;
    pybind11::class_<GrampcBinding::grampc_param>(binding, "_grampc_param")
        .def(pybind11::init<>())
        .def_readonly("Nx", &prefix_param::Nx)
        .def_readonly("Nu", &prefix_param::Nu)
        .def_readonly("Np", &prefix_param::Np)
        .def_readonly("Ng", &prefix_param::Ng)
        .def_readonly("Nh", &prefix_param::Nh)
        .def_readonly("NgT", &prefix_param::NgT)
        .def_readonly("NhT", &prefix_param::NhT)
        .def_readonly("Nc", &prefix_param::Nc)

        .def_readonly("x0", &prefix_param::x0)
        .def_readonly("xdes", &prefix_param::xdes)

        .def_readonly("u0", &prefix_param::u0)
        .def_readonly("udes", &prefix_param::udes)
        .def_readonly("umax", &prefix_param::umax)
        .def_readonly("umin", &prefix_param::umin)

        .def_readonly("p0", &prefix_param::p0)
        .def_readonly("pmax", &prefix_param::pmax)
        .def_readonly("pmin", &prefix_param::pmin)

        .def_readonly("Thor", &prefix_param::Thor)
        .def_readonly("Tmax", &prefix_param::Tmax)
        .def_readonly("Tmin", &prefix_param::Tmin)
        .def_readonly("dt", &prefix_param::dt)
        .def_readonly("t0", &prefix_param::t0);

    typedef GrampcBinding::grampc_opt prefix_opt;
    pybind11::class_<GrampcBinding::grampc_opt>(binding, "_grampc_opt")
        .def(pybind11::init<>())
        .def_readonly("Nhor", &prefix_opt::Nhor)
        .def_readonly("MaxGradIter", &prefix_opt::MaxGradIter)
        .def_readonly("MaxMultIter", &prefix_opt::MaxMultIter)
        .def_readonly("ShiftControl", &prefix_opt::ShiftControl)

        .def_readonly("TimeDiscretization", &prefix_opt::TimeDiscretization)

        .def_readonly("IntegralCost", &prefix_opt::IntegralCost)
        .def_readonly("TerminalCost", &prefix_opt::TerminalCost)
        .def_readonly("IntegratorCost", &prefix_opt::IntegratorCost)

        .def_readonly("Integrator", &prefix_opt::Integrator)
        .def_readonly("IntegratorRelTol", &prefix_opt::IntegratorRelTol)
        .def_readonly("IntegratorAbsTol", &prefix_opt::IntegratorAbsTol)
        .def_readonly("IntegratorMinStepSize", &prefix_opt::IntegratorMinStepSize)
        .def_readonly("IntegratorMaxSteps", &prefix_opt::IntegratorMaxSteps)
        .def_readonly("FlagsRodas", &prefix_opt::FlagsRodas)

        .def_readonly("LineSearchType", &prefix_opt::LineSearchType)
        .def_readonly("LineSearchExpAutoFallback", &prefix_opt::LineSearchExpAutoFallback)
        .def_readonly("LineSearchMax", &prefix_opt::LineSearchMax)
        .def_readonly("LineSearchMin", &prefix_opt::LineSearchMin)
        .def_readonly("LineSearchInit", &prefix_opt::LineSearchInit)
        .def_readonly("LineSearchAdaptAbsTol", &prefix_opt::LineSearchAdaptAbsTol)
        .def_readonly("LineSearchAdaptFactor", &prefix_opt::LineSearchAdaptFactor)
        .def_readonly("LineSearchIntervalTol", &prefix_opt::LineSearchIntervalTol)
        .def_readonly("LineSearchIntervalFactor", &prefix_opt::LineSearchIntervalFactor)

        .def_readonly("OptimControl", &prefix_opt::OptimControl)
        .def_readonly("OptimParam", &prefix_opt::OptimParam)
        .def_readonly("OptimParamLineSearchFactor", &prefix_opt::OptimParamLineSearchFactor)
        .def_readonly("OptimTime", &prefix_opt::OptimTime)
        .def_readonly("OptimTimeLineSearchFactor", &prefix_opt::OptimTimeLineSearchFactor)

        .def_readonly("ScaleProblem", &prefix_opt::ScaleProblem)
        .def_readonly("xScale", &prefix_opt::xScale)
        .def_readonly("xOffset", &prefix_opt::xOffset)
        .def_readonly("uScale", &prefix_opt::uScale)
        .def_readonly("uOffset", &prefix_opt::uOffset)
        .def_readonly("pScale", &prefix_opt::pScale)
        .def_readonly("pOffset", &prefix_opt::pOffset)
        .def_readonly("TScale", &prefix_opt::TScale)
        .def_readonly("TOffset", &prefix_opt::TOffset)
        .def_readonly("JScale", &prefix_opt::JScale)
        .def_readonly("cScale", &prefix_opt::cScale)

        .def_readonly("EqualityConstraints", &prefix_opt::EqualityConstraints)
        .def_readonly("InequalityConstraints", &prefix_opt::InequalityConstraints)
        .def_readonly("TerminalEqualityConstraints", &prefix_opt::TerminalEqualityConstraints)
        .def_readonly("TerminalInequalityConstraints", &prefix_opt::TerminalInequalityConstraints)
        .def_readonly("ConstraintsHandling", &prefix_opt::ConstraintsHandling)
        .def_readonly("ConstraintsAbsTol", &prefix_opt::ConstraintsAbsTol)

        .def_readonly("MultiplierMax", &prefix_opt::MultiplierMax)
        .def_readonly("MultiplierDampingFactor", &prefix_opt::MultiplierDampingFactor)
        .def_readonly("PenaltyMax", &prefix_opt::PenaltyMax)
        .def_readonly("PenaltyMin", &prefix_opt::PenaltyMin)
        .def_readonly("PenaltyIncreaseFactor", &prefix_opt::PenaltyIncreaseFactor)
        .def_readonly("PenaltyDecreaseFactor", &prefix_opt::PenaltyDecreaseFactor)
        .def_readonly("PenaltyIncreaseThreshold", &prefix_opt::PenaltyIncreaseThreshold)
        .def_readonly("AugLagUpdateGradientRelTol", &prefix_opt::AugLagUpdateGradientRelTol)

        .def_readonly("ConvergenceCheck", &prefix_opt::ConvergenceCheck)
        .def_readonly("ConvergenceGradientRelTol", &prefix_opt::ConvergenceGradientRelTol);

    typedef GrampcBinding::grampc_sol prefix_sol;
    pybind11::class_<GrampcBinding::grampc_sol>(binding, "_grampc_sol")
        .def(pybind11::init<>())
        .def_readonly("xnext", &prefix_sol::xnext)
        .def_readonly("unext", &prefix_sol::unext)
        .def_readonly("pnext", &prefix_sol::pnext)
        .def_readonly("Tnext", &prefix_sol::Tnext)
        .def_readonly("J", &prefix_sol::J)
        .def_readonly("cfct", &prefix_sol::cfct)
        .def_readonly("pen", &prefix_sol::pen)
        .def_readonly("iter", &prefix_sol::iter)
        .def_readonly("status", &prefix_sol::status);

    typedef GrampcBinding::grampc_rws prefix_rws;
    pybind11::class_<GrampcBinding::grampc_rws>(binding, "_grampc_rws")
        .def(pybind11::init<>())
        .def_readonly("t", &prefix_rws::t)
        .def_readonly("tls", &prefix_rws::tls)

        .def_readonly("x", &prefix_rws::x)
        .def_readonly("adj", &prefix_rws::adj)
        .def_readonly("dcdx", &prefix_rws::dcdx)

        .def_readonly("u", &prefix_rws::u)
        .def_readonly("uls", &prefix_rws::uls)
        .def_readonly("uprev", &prefix_rws::uprev)
        .def_readonly("gradu", &prefix_rws::gradu)
        .def_readonly("graduprev", &prefix_rws::graduprev)
        .def_readonly("dcdu", &prefix_rws::dcdu)

        .def_readonly("p", &prefix_rws::p)
        .def_readonly("pls", &prefix_rws::pls)
        .def_readonly("pprev", &prefix_rws::pprev)
        .def_readonly("gradp", &prefix_rws::gradp)
        .def_readonly("gradpprev", &prefix_rws::gradpprev)
        .def_readonly("dcdp", &prefix_rws::dcdp)

        .def_readonly("T", &prefix_rws::T)
        .def_readonly("Tprev", &prefix_rws::Tprev)
        .def_readonly("gradT", &prefix_rws::gradT)
        .def_readonly("dcdt", &prefix_rws::dcdt)

        .def_readonly("mult", &prefix_rws::mult)
        .def_readonly("pen", &prefix_rws::pen)
        .def_readonly("cfct", &prefix_rws::cfct)
        .def_readonly("cfctprev", &prefix_rws::cfctprev)
        .def_readonly("cfctAbsTol", &prefix_rws::cfctAbsTol)

        .def_readonly("lsAdapt", &prefix_rws::lsAdapt)
        .def_readonly("lsExplicit", &prefix_rws::lsExplicit)
        .def_readonly("rwsScale", &prefix_rws::rwsScale)
        .def_readonly("lrwsGeneral", &prefix_rws::lrwsGeneral)
        .def_readonly("rwsGeneral", &prefix_rws::rwsGeneral)

        .def_readonly("lworkRodas", &prefix_rws::lworkRodas)
        .def_readonly("liworkRodas", &prefix_rws::liworkRodas)
        .def_readonly("rparRodas", &prefix_rws::rparRodas)
        .def_readonly("iparRodas", &prefix_rws::iparRodas)
        .def_readonly("workRodas", &prefix_rws::workRodas)
        .def_readonly("iworkRodas", &prefix_rws::iworkRodas);

    pybind11::class_<ProblemBase, PyProblem>(m, "ProblemBase")
        .def(pybind11::init<>())
        .def_readwrite("Nx", &ProblemBase::Nx_)
        .def_readwrite("Nu", &ProblemBase::Nu_)
        .def_readwrite("Np", &ProblemBase::Np_)
        .def_readwrite("Ng", &ProblemBase::Ng_)
        .def_readwrite("Nh", &ProblemBase::Nh_)
        .def_readwrite("NgT", &ProblemBase::NgT_)
        .def_readwrite("NhT", &ProblemBase::NhT_)

        .def("ffct", &ProblemBase::ffct)
        .def("dfdx_vec", &ProblemBase::dfdx_vec)
        .def("dfdu_vec", &ProblemBase::dfdu_vec)
        .def("dfdp_vec", &ProblemBase::dfdp_vec)

        .def("lfct", &ProblemBase::lfct)
        .def("dldx", &ProblemBase::dldx)
        .def("dldu", &ProblemBase::dldu)
        .def("dldp", &ProblemBase::dldp)

        .def("Vfct", &ProblemBase::Vfct)
        .def("dVdx", &ProblemBase::dVdx)
        .def("dVdp", &ProblemBase::dVdp)
        .def("dVdT", &ProblemBase::dVdT)

        .def("gfct", &ProblemBase::gfct)
        .def("dgdx_vec", &ProblemBase::dgdx_vec)
        .def("dgdu_vec", &ProblemBase::dgdu_vec)
        .def("dgdp_vec", &ProblemBase::dgdp_vec)

        .def("hfct", &ProblemBase::hfct)
        .def("dhdx_vec", &ProblemBase::dhdx_vec)
        .def("dhdu_vec", &ProblemBase::dhdu_vec)
        .def("dhdp_vec", &ProblemBase::dhdp_vec)

        .def("gTfct", &ProblemBase::gTfct)
        .def("dgTdx_vec", &ProblemBase::dgTdx_vec)
        .def("dgTdp_vec", &ProblemBase::dgTdp_vec)
        .def("dgTdT_vec", &ProblemBase::dgTdT_vec)

        .def("hTfct", &ProblemBase::hTfct)
        .def("dhTdx_vec", &ProblemBase::dhTdx_vec)
        .def("dhTdp_vec", &ProblemBase::dhTdp_vec)
        .def("dhTdT_vec", &ProblemBase::dhTdT_vec)

        .def("dfdx", &ProblemBase::dfdx)
        .def("dfdxtrans", &ProblemBase::dfdxtrans)
        .def("dfdt", &ProblemBase::dfdt)
        .def("dHdxdt", &ProblemBase::dHdxdt)
        .def("Mfct", &ProblemBase::Mfct)
        .def("Mtrans", &ProblemBase::Mtrans);
}
