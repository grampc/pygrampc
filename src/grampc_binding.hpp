/* This file is part of PyGRAMPC - (https://github.com/grampc/pygrampc)
 *
 * PyGRAMPC -- A Python interface for the GRAMPC solver
 *
 * Copyright 2023 by Thore Wietzke and Andreas Voelz
 * All rights reserved.
 *
 * PyGRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 */

#ifndef GRAMPC_BINDINGS_HPP
#define GRAMPC_BINDINGS_HPP

extern "C" 
{
    #include "grampc.h"
}
#include "problem_base.hpp"
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <chrono>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, 1> Vector;
typedef Eigen::Matrix<typeRNum, Eigen::Dynamic, Eigen::Dynamic> Matrix;

class GrampcBinding
{
    public:

    struct grampc_param
    {
    public:
        grampc_param();
        void reMapMemory(const typeGRAMPC* grampc);
    public:
        int* Nx;
        int* Nu;
        int* Np;
        int* Ng;
        int* Nh;
        int* NgT;
        int* NhT;
        int* Nc;

        Eigen::Map<Vector> x0;
        Eigen::Map<Vector> xdes;

        Eigen::Map<Vector> u0;
        Eigen::Map<Vector> udes;
        Eigen::Map<Vector> umax;
        Eigen::Map<Vector> umin;

        Eigen::Map<Vector> p0;
        Eigen::Map<Vector> pmax;
        Eigen::Map<Vector> pmin;

        typeRNum* Thor;
        typeRNum* Tmax;
        typeRNum* Tmin;

        typeRNum* dt;
        typeRNum* t0;
    };
    
    /* OPTIONS STRUCTURE ******************************************************/
    // Options for GRAMPC computations, which are part of the GRAMPC data struct.
    struct grampc_opt
    {
    public:
        grampc_opt();
        void reMapMemory(const typeGRAMPC* grampc);
    public:
        int* Nhor;
        int* MaxGradIter;
        int* MaxMultIter;
        int* ShiftControl;

        int* TimeDiscretization;

        int* IntegralCost;
        int* TerminalCost;
        int* IntegratorCost;

        int* Integrator;
        typeRNum* IntegratorRelTol;
        typeRNum* IntegratorAbsTol;
        typeRNum* IntegratorMinStepSize;
        int*    IntegratorMaxSteps;
        std::vector<int> FlagsRodas;

        int* LineSearchType;
        int* LineSearchExpAutoFallback;
        typeRNum* LineSearchMax;
        typeRNum* LineSearchMin;
        typeRNum* LineSearchInit;
        typeRNum* LineSearchAdaptAbsTol;
        typeRNum* LineSearchAdaptFactor;
        typeRNum* LineSearchIntervalTol;
        typeRNum* LineSearchIntervalFactor;

        int* OptimControl;
        int* OptimParam;
        typeRNum* OptimParamLineSearchFactor;
        int* OptimTime;
        typeRNum* OptimTimeLineSearchFactor;

        int* ScaleProblem;
        Eigen::Map<Vector> xScale;
        Eigen::Map<Vector> xOffset;
        Eigen::Map<Vector> uScale;
        Eigen::Map<Vector> uOffset;
        Eigen::Map<Vector> pScale;
        Eigen::Map<Vector> pOffset;
        typeRNum* TScale;
        typeRNum* TOffset;
        typeRNum* JScale;
        Eigen::Map<Vector> cScale;

        int* EqualityConstraints;
        int* InequalityConstraints;
        int* TerminalEqualityConstraints;
        int* TerminalInequalityConstraints;
        int* ConstraintsHandling;
        Eigen::Map<Vector> ConstraintsAbsTol;

        typeRNum* MultiplierMax;
        typeRNum* MultiplierDampingFactor;
        typeRNum* PenaltyMax;
        typeRNum* PenaltyMin;
        typeRNum* PenaltyIncreaseFactor;
        typeRNum* PenaltyDecreaseFactor;
        typeRNum* PenaltyIncreaseThreshold;
        typeRNum* AugLagUpdateGradientRelTol;

        int* ConvergenceCheck;
        typeRNum* ConvergenceGradientRelTol;
    };

    /* SOLUTION STRUCTURE ******************************************************/
    // Solution structure of GRAMPC computations, which is part of the GRAMPC data struct.
    struct grampc_sol
    {
    public:
        grampc_sol();
        void reMapMemory(const typeGRAMPC* grampc);
    public:
        Eigen::Map<Vector> xnext;
        Eigen::Map<Vector> unext;
        Eigen::Map<Vector> pnext;
        typeRNum* Tnext;
        Eigen::Map<Vector> J;
        typeRNum* cfct;
        typeRNum* pen;
        std::vector<int> iter;
        int status;
    };

    /* RWS STRUCTURE **********************************************************/
    // Real workspace structure of GRAMPC that holds intermediate results
    // and computation results.
    struct grampc_rws
    {
    public:
        grampc_rws();
        void reMapMemory(const typeGRAMPC* grampc);
    public:
        Eigen::Map<Vector> t;
        Eigen::Map<Vector> tls;

        Eigen::Map<Matrix> x;
        Eigen::Map<Matrix> adj;
        Eigen::Map<Matrix> dcdx;

        Eigen::Map<Matrix> u;
        Eigen::Map<Matrix> uls;
        Eigen::Map<Matrix> uprev;
        Eigen::Map<Matrix> gradu;
        Eigen::Map<Matrix> graduprev;
        Eigen::Map<Matrix> dcdu;

        Eigen::Map<Matrix> p;
        Eigen::Map<Matrix> pls;
        Eigen::Map<Matrix> pprev;
        Eigen::Map<Matrix> gradp;
        Eigen::Map<Matrix> gradpprev;
        Eigen::Map<Matrix> dcdp;

        typeRNum* T;
        typeRNum* Tprev;
        typeRNum* gradT;
        typeRNum* gradTprev;
        typeRNum* dcdt;

        Eigen::Map<Matrix> mult;
        Eigen::Map<Matrix> pen;
        Eigen::Map<Matrix> cfct;
        Eigen::Map<Matrix> cfctprev;
        Eigen::Map<Vector> cfctAbsTol;

        Eigen::Map<Vector> lsAdapt;
        Eigen::Map<Vector> lsExplicit;
        Eigen::Map<Vector> rwsScale;
        int* lrwsGeneral;
        Eigen::Map<Vector> rwsGeneral;

        int* lworkRodas;
        int* liworkRodas;
        Eigen::Map<Vector> rparRodas;
        std::vector<int> iparRodas;
        Eigen::Map<Vector> workRodas;
        std::vector<int> iworkRodas;
    };

    public:
        grampc_param param;
        grampc_opt opt;
        grampc_sol sol;
        grampc_rws rws;

        // id of the concrete problem
        ProblemBase* problem_binding;

    public:
        // Create Python interface to the GRAMPC solver for problem description
        GrampcBinding(ProblemBase* problem);
        ~GrampcBinding();

        // Calls grampc_run and returns the wall clock time of one call in milliseconds.
        typeRNum run();

        // Wraps the grampc c-function grampc_estim_penmin.
        void estim_penmin(bool run_grampc);

        // ---- PARAMETER SETTINGS ----
        // Sets the parameter of type typeRNum value given by the literal string key.
        void set_param_real(const std::string& key, typeRNum value);
        // Sets the parameter of type typeRNum-array given by the literal string key.
        void set_param_real_vec(const std::string& key, const std::vector<typeRNum>& values);

        // ---- OPTIONS SETTINGS ----
        // Sets the option of type string given by the literal string key.
        void set_opt_str(const std::string& key, const std::string& vstr);
        // Sets the option of type integer given by the literal string key.
        void set_opt_int(const std::string& key, int value);
        // Sets the option of type integer-array given by the literal string key.
        void set_opt_int_vec(const std::string& key, const std::vector<int>& values);
        // Sets the option of type typeRNum value given by the literal string key.
        void set_opt_real(const std::string& key, typeRNum value);
        // Sets the option of type typeRNum array given by the literal string key.
        void set_opt_real_vec(const std::string& key, const std::vector<typeRNum>& values);

        // access the probFct like the Matlab interface
        Vector ffct(typeRNum t, Vector x, Vector u, Vector p);

        typeRNum lfct(typeRNum t, Vector x, Vector u, Vector p, Vector xdes, Vector udes);

        typeRNum Vfct(typeRNum t, Vector x, Vector p, Vector xdes);

        Vector gfct(typeRNum t, Vector x, Vector u, Vector p);
        Vector hfct(typeRNum t, Vector x, Vector u, Vector p);
        
        Vector gTfct(typeRNum t, Vector x, Vector p);
        Vector hTfct(typeRNum t, Vector x, Vector p);

        // Prints the parameter settings to standard console.
        void print_params();
        // Prints the options settings to standard console.
        void print_opts();
        // Prints the most recent status variable values to standard console.
        void print_status();

    private: 
        // grampc handler
        typeGRAMPC* grampc_;

};

#endif