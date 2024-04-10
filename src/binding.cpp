/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#include <grampc_interface.hpp>
#include <problem_description.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

PYBIND11_MODULE(_core, m)
{
    pybind11::class_<GrampcBinding> binding(m, "GrampcBinding");
    binding.def(pybind11::init<ProblemDescription *>())
        .def_readonly("param", &GrampcBinding::param)
        .def_readonly("opt", &GrampcBinding::opt)
        .def_readonly("sol", &GrampcBinding::sol)
        .def_readonly("rws", &GrampcBinding::rws)
        .def_readwrite("problem", &GrampcBinding::problem_description)
        .def("run", &GrampcBinding::run)
        .def("_set_param_real", &GrampcBinding::set_param_real)
        .def("_set_param_real_vec", &GrampcBinding::set_param_real_vec)
        .def("_set_opt_str", &GrampcBinding::set_opt_str)
        .def("_set_opt_int", &GrampcBinding::set_opt_int)
        .def("_set_opt_int_vec", &GrampcBinding::set_opt_int_vec)
        .def("_set_opt_real", &GrampcBinding::set_opt_real)
        .def("_set_opt_real_vec", &GrampcBinding::set_opt_real_vec)
        .def("set_rws_u", &GrampcBinding::set_rws_u)
        .def("set_rws_multiplier", &GrampcBinding::set_rws_multiplier)
        .def("set_rws_penalty", &GrampcBinding::set_rws_penalty)
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

    pybind11::class_<ProblemDescription, PyProblem>(m, "ProblemDescription")
        .def(pybind11::init<>())
        .def_readwrite("Nx", &ProblemDescription::Nx)
        .def_readwrite("Nu", &ProblemDescription::Nu)
        .def_readwrite("Np", &ProblemDescription::Np)
        .def_readwrite("Ng", &ProblemDescription::Ng)
        .def_readwrite("Nh", &ProblemDescription::Nh)
        .def_readwrite("NgT", &ProblemDescription::NgT)
        .def_readwrite("NhT", &ProblemDescription::NhT)

        .def("ffct", &ProblemDescription::ffct)
        .def("dfdx_vec", &ProblemDescription::dfdx_vec)
        .def("dfdu_vec", &ProblemDescription::dfdu_vec)
        .def("dfdp_vec", &ProblemDescription::dfdp_vec)

        .def("lfct", &ProblemDescription::lfct)
        .def("dldx", &ProblemDescription::dldx)
        .def("dldu", &ProblemDescription::dldu)
        .def("dldp", &ProblemDescription::dldp)

        .def("Vfct", &ProblemDescription::Vfct)
        .def("dVdx", &ProblemDescription::dVdx)
        .def("dVdp", &ProblemDescription::dVdp)
        .def("dVdT", &ProblemDescription::dVdT)

        .def("gfct", &ProblemDescription::gfct)
        .def("dgdx_vec", &ProblemDescription::dgdx_vec)
        .def("dgdu_vec", &ProblemDescription::dgdu_vec)
        .def("dgdp_vec", &ProblemDescription::dgdp_vec)

        .def("hfct", &ProblemDescription::hfct)
        .def("dhdx_vec", &ProblemDescription::dhdx_vec)
        .def("dhdu_vec", &ProblemDescription::dhdu_vec)
        .def("dhdp_vec", &ProblemDescription::dhdp_vec)

        .def("gTfct", &ProblemDescription::gTfct)
        .def("dgTdx_vec", &ProblemDescription::dgTdx_vec)
        .def("dgTdp_vec", &ProblemDescription::dgTdp_vec)
        .def("dgTdT_vec", &ProblemDescription::dgTdT_vec)

        .def("hTfct", &ProblemDescription::hTfct)
        .def("dhTdx_vec", &ProblemDescription::dhTdx_vec)
        .def("dhTdp_vec", &ProblemDescription::dhTdp_vec)
        .def("dhTdT_vec", &ProblemDescription::dhTdT_vec)

        .def("dfdx", &ProblemDescription::dfdx)
        .def("dfdxtrans", &ProblemDescription::dfdxtrans)
        .def("dfdt", &ProblemDescription::dfdt)
        .def("dHdxdt", &ProblemDescription::dHdxdt)
        .def("Mfct", &ProblemDescription::Mfct)
        .def("Mtrans", &ProblemDescription::Mtrans);
}