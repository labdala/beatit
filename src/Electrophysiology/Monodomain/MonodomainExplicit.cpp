/*
 ============================================================================

 .______    _______     ___   .___________.    __  .___________.
 |   _  \  |   ____|   /   \  |           |   |  | |           |
 |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
 |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |     
 |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |     
 |______/  |_______/__/     \__\  |__|        |__|     |__|     
 
 BeatIt - code for cardiovascular simulations
 Copyright (C) 2016 Simone Rossi

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ============================================================================
 */

/**
 * \file MonodomainExplicit.cpp
 *
 * \class MonodomainExplicit
 *
 * \brief This class provides a simple factory implementation
 *
 * For details on how to use it check the test_factory in the testsuite folder
 *
 *
 * \author srossi
 *
 * \version 0.0
 *
 *
 * Contact: srossi@gmail.com
 *
 * Created on: Aug 11, 2016
 *
 */
#include <exception>
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/type_tensor.h"
#include "PoissonSolver/Poisson.hpp"

// Include files that define a simple steady system
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
//
#include "libmesh/vector_value.h"
#include "libmesh/linear_solver.h"
// Define the Finite Element object.
#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"

// Define useful datatypes for finite element
// matrix and vector components.
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

#include "libmesh/exodusII_io.h"
//#include "libmesh/exodusII_io_helper.h"
#include "libmesh/gmv_io.h"

#include "libmesh/perf_log.h"

#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/fourth_error_estimators.h"
#include "libmesh/discontinuity_measure.h"
#include <sys/stat.h>

#include "Electrophysiology/IonicModels/NashPanfilov.hpp"
#include "Electrophysiology/IonicModels/Grandi11.hpp"
#include "Electrophysiology/IonicModels/ORd.hpp"
#include "Electrophysiology/IonicModels/TP06.hpp"
#include "Electrophysiology/IonicModels/Cubic.hpp"

#include "Electrophysiology/Monodomain/MonodomainExplicit.hpp"
#include "Util/SpiritFunction.hpp"

#include "libmesh/discontinuity_measure.h"
#include "libmesh/fe_interface.h"

#include "Electrophysiology/Pacing/PacingProtocolSpirit.hpp"

#include <cstdlib>

namespace BeatIt
{

// ///////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////
typedef libMesh::TransientLinearImplicitSystem MonodomainSystem;
typedef libMesh::TransientLinearImplicitSystem ElectroSystem;
typedef libMesh::TransientExplicitSystem IonicModelSystem;
typedef libMesh::ExplicitSystem ParameterSystem;

ElectroSolver* createMonodomainExplicit(libMesh::EquationSystems &es)
{
    return new MonodomainExplicit(es);
}

MonodomainExplicit::MonodomainExplicit(libMesh::EquationSystems &es) :
        ElectroSolver(es, "monodomain_explicit")
{

}

MonodomainExplicit::~MonodomainExplicit()
{
}

void MonodomainExplicit::setup_systems(GetPot &data, std::string section)
{
    std::cout << "* MONODOMAIN EXP: Setting up data from section " << section << std::endl;
    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Starts by creating the equation systems
    // 1) ADR
    std::cout << "* MONODOMAIN EXP: Creating new System for the monodomain diffusion reaction equation" << std::endl;
    ElectroSystem &monodomain_system = M_equationSystems.add_system < ElectroSystem > (M_model);
    // TO DO: Generalize to higher order
    monodomain_system.add_variable("Q", M_order, M_FEFamily);
    // Add 3 matrices
    monodomain_system.add_matrix("lumped_mass");
    monodomain_system.add_matrix("high_order_mass");
    // Add lumped mass matrix
    monodomain_system.add_vector("lumped_mass_vector");
    monodomain_system.add_matrix("mass");
    monodomain_system.add_matrix("stiffness");
    monodomain_system.add_vector("aux1");
    monodomain_system.add_vector("aux2");

    M_exporterNames.insert(M_model);

    monodomain_system.init();

    // WAVE
    ElectroSystem &wave_system = M_equationSystems.add_system < ElectroSystem > ("wave");
    wave_system.add_variable("V", M_order, M_FEFamily);
    M_exporterNames.insert("wave");
    wave_system.init();

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // 2) ODEs
    setup_ODE_systems(data, section);
    // Add the applied current to this system
    std::cout << "* MONODOMAIN EXP: Creating auxiliary explicit systems " << std::endl;

//    M_ionicModelExporterNames.insert("istim");
    IonicModelSystem &cut_system = M_equationSystems.add_system < IonicModelSystem > ("cut");
    cut_system.add_variable("cut", M_order, M_FEFamily);
    cut_system.init();

    // ///////////////////////////////////////////////////////////////////////
    // ///////////////////////////////////////////////////////////////////////
    // Distributed Parameters
    std::cout << "* MONODOMAIN EXP: Creating parameters spaces " << std::endl;
    ParameterSystem &conductivity_system = M_equationSystems.add_system < ParameterSystem > ("conductivity");
    conductivity_system.add_variable("Dff", libMesh::CONSTANT, libMesh::MONOMIAL);
    conductivity_system.add_variable("Dss", libMesh::CONSTANT, libMesh::MONOMIAL);
    conductivity_system.add_variable("Dnn", libMesh::CONSTANT, libMesh::MONOMIAL);
    M_parametersExporterNames.insert("conductivity");
    std::cout << "* MONODOMAIN EXP: Initializing equation systems " << std::endl;
    // Initializing
    conductivity_system.init();
//    M_equationSystems.init();
    M_equationSystems.print_info();

    // Conductivity Tensor in local coordinates
    // Fiber direction
    // Default Value = 1.3342  kOhm^-1 cm^-1
    std::cout << section + "/Dff" << std::endl;
    double Dff = M_datafile(section + "/Dff", 1.3342);
    // Sheet direction
    // Default Value = 1.3342  kOhm^-1 cm^-1
    double Dss = M_datafile(section + "/Dss", 0.17606);
    // Cross fiber direction
    // Default Value = 1.3342  kOhm^-1 cm^-1
    double Dnn = M_datafile(section + "/Dnn", 0.17606);
    double Chi = M_datafile(section + "/Chi", 1400.0);
    M_equationSystems.parameters.set<double>("Chi") = Chi;
    double tau = M_datafile(section + "/tau", 0.0);
    M_equationSystems.parameters.set<double>("tau") = tau;
    if (tau >= 0)
        M_equationType = EquationType::Wave;
    std::string anisotropy = M_datafile(section + "/anisotropy", "orhotropic");
    std::map<std::string, Anisotropy> aniso_map;

    aniso_map["orthotropic"] = Anisotropy::Orthotropic;
    aniso_map["isotropic"] = Anisotropy::Isotropic;
    aniso_map["transverse"] = Anisotropy::TransverselyIsotropic;

    auto it = aniso_map.find(anisotropy);
    if( it != aniso_map.end())
    {
        M_anisotropy = it->second;
    }
    std::cout << "* MONODOMAIN EXP: Parameters: " << std::endl;
    std::cout << "              Chi = " << Chi << std::endl;
    std::cout << "              Dff = " << Dff << std::endl;
    std::cout << "              Dss = " << Dss << std::endl;
    std::cout << "              Dnn = " << Dnn << std::endl;
    std::cout << "              tau = " << tau << std::endl;
    std::cout << "              anisotropy = " << anisotropy << std::endl;
}

void MonodomainExplicit::init_systems(double time)
{
}

void MonodomainExplicit::cut(double time, std::string f)
{
    std::cout << "* MONODOMAIN EXP: Cutting potential" << std::endl;

    SpiritFunction func;
    func.read(f);
    // Add the applied current to this system
    IonicModelSystem &cut_system = M_equationSystems.get_system < IonicModelSystem > ("cut");
    cut_system.time = time;
    cut_system.project_solution(&func);
    cut_system.solution->close();
//    cut_system.update();
    IonicModelSystem &ionic_model_system = M_equationSystems.get_system < IonicModelSystem > ("ionic_model");
    // WAVE
    ElectroSystem &wave_system = M_equationSystems.add_system < ElectroSystem > ("wave");
    ElectroSystem &monodomain_system = M_equationSystems.get_system < ElectroSystem > (M_model);

    auto first = cut_system.solution->first_local_index();
    auto last = cut_system.solution->last_local_index();
    int n_vars = ionic_model_system.n_vars();
    std::vector<double> init_val(n_vars + 1, 0.0);
    double cut_value = 0.0;

    // TODO: this is not going to work anymore
    //M_ionicModelMap->initialize(init_val);
    init_val[0] = -85.0;

    ionic_model_system.solution->close();
    wave_system.solution->close();
    monodomain_system.solution->close();
    cut_system.solution->close();

    for (int index = first; index < last; index++)
    {
        cut_value = (*cut_system.solution)(index);
        if (cut_value < 0.1)
        {
            if (M_equationType == EquationType::ReactionDiffusion)
            {
                wave_system.solution->set(index, init_val[0]);
                monodomain_system.solution->set(index, init_val[0]);
            }
            else
            {
                wave_system.solution->set(index, init_val[0]);
                monodomain_system.solution->set(index, 0.0);
            }
//       std::cout << "ion value" << std::endl;
            for (int m = 0; m < n_vars; m++)
            {
                int var_index = index * n_vars + m;
                ionic_model_system.solution->set(var_index, init_val[m + 1]);
            }
        }
    }

    ionic_model_system.solution->close();
    wave_system.solution->close();
    monodomain_system.solution->close();
    cut_system.solution->close();

}

void MonodomainExplicit::generate_fibers(const GetPot &data, const std::string &section)
{
    std::cout << "* MONODOMAIN EXP: Creating fiber fields" << std::endl;

    std::cout << "* MONODOMAIN EXP: Solving Poisson porblem" << std::endl;

    libMesh::Mesh new_mesh(dynamic_cast<libMesh::Mesh&>(M_equationSystems.get_mesh()));
    libMesh::EquationSystems es(new_mesh);
    Poisson p(es);
    p.setup(data, section);
    p.assemble_system();
    p.solve_system();
    p.compute_elemental_solution_gradient();
    p.save_exo();

    const auto &sol_ptr = p.get_P0_solution();

    ParameterSystem &fiber_system = M_equationSystems.get_system < ParameterSystem > ("fibers");
    ParameterSystem &sheets_system = M_equationSystems.get_system < ParameterSystem > ("sheets");
    ParameterSystem &xfiber_system = M_equationSystems.get_system < ParameterSystem > ("xfibers");

    *sheets_system.solution = *p.get_gradient();

    auto first_local_index = fiber_system.solution->first_local_index();
    auto last_local_index = fiber_system.solution->last_local_index();
    auto first_local_index_sol = sol_ptr->first_local_index();
    auto last_local_index_sol = sol_ptr->last_local_index();

    double norm = 0.0;
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    double cx = data(section + "/centerline_x", 0.0);
    double cy = data(section + "/centerline_y", 0.0);
    double cz = data(section + "/centerline_z", 1.0);
    double cdot = 0.0;
    double sx = 0.0;
    double sy = 0.0;
    double sz = 0.0;
    double xfx = 0.0;
    double xfy = 0.0;
    double xfz = 0.0;
    double epi_angle = data(section + "/epi_angle", -60.0);
    double endo_angle = data(section + "/endo_angle", 60.0);
    double potential = 0.0;
    double W11 = 0.0;
    double W12 = 0.0;
    double W13 = 0.0;
    double W21 = 0.0;
    double W22 = 0.0;
    double W23 = 0.0;
    double W31 = 0.0;
    double W32 = 0.0;
    double W33 = 0.0;
    //
    double R11 = 0.0;
    double R12 = 0.0;
    double R13 = 0.0;
    double R21 = 0.0;
    double R22 = 0.0;
    double R23 = 0.0;
    double R31 = 0.0;
    double R32 = 0.0;
    double R33 = 0.0;
    double sa = 0.0;
    double sa2 = 0.0;
    double teta1 = 0.0;
    double teta2 = 0.0;
    double teta = 0.0;
    double m = 0.0;
    double q = 0.0;
    double f0x = 0.0;
    double f0y = 0.0;
    double f0z = 0.0;

    auto normalize = [](double &x, double &y, double &z, double X, double Y, double Z)
    {
        double norm = std::sqrt( x * x + y * y + z * z);
        if(norm >= 1e-12 )
        {
            x /= norm;
            y /= norm;
            z /= norm;
        }
        else
        {
            x = X;
            y = Y;
            z = Z;
        }
    };

    std::cout << "* MONODOMAIN EXP: Computing rule-based fiber fields" << std::endl;

//    sol_ptr->print();

    auto j = first_local_index_sol;

    for (auto i = first_local_index; i < last_local_index;)
    {

//        std::cout << "* MONODOMAIN: getting sol ... " << std::flush;
        potential = (*sol_ptr)(j);
        j++;
//        std::cout << " done" << std::endl;

        sx = (*sheets_system.solution)(i);
        sy = (*sheets_system.solution)(i + 1);
        sz = (*sheets_system.solution)(i + 2);
        normalize(sx, sy, sz, 0.0, 1.0, 0.0);

//        norm = std::sqrt( sx * sx + sy * sy + sz * sz);
//        if(norm >= 1e-12 )
//        {
//            sx /= norm;
//            sy /= norm;
//            sz /= norm;
//        }
//        else
//        {
//            sx = 1.0;
//            sy = 0.0;
//            sz = 0.0;
//        }

        cdot = cx * sx + cy * sy + cz * sz;

        xfx = cx - cdot * sx;
        xfy = cy - cdot * sy;
        xfz = cz - cdot * sz;
        normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

        fx = sy * xfz - sz * xfy;
        fy = sz * xfx - sx * xfz;
        fz = sx * xfy - sy * xfx;
        normalize(fx, fy, fz, 1.0, 0.0, 0.0);

        teta1 = M_PI * epi_angle / 180.0;
        teta2 = M_PI * endo_angle / 180.0;
        m = (teta1 - teta2);
        q = teta2;
        teta = m * potential + q;

        //*************************************************************//
        // The fiber field F is a rotation of the flat fiber field f
        // F = R f
        // where R is the rotation matrix.
        // To compute R we need the sin(teta) and
        // the sin(teta)^2 and the cross-product matrix W (check
        // rodrigues formula on wikipedia :) )
        //*************************************************************//
        sa = std::sin(teta);
        sa2 = 2.0 * std::sin(0.5 * teta) * std::sin(0.5 * teta);

        W11 = 0.0;
        W12 = -sz;
        W13 = sy;
        W21 = sz;
        W22 = 0.0;
        W23 = -sx;
        W31 = -sy;
        W32 = sx;
        W33 = 0.0;
        //
        R11 = 1.0 + sa * W11 + sa2 * (sx * sx - 1.0);
        R12 = 0.0 + sa * W12 + sa2 * (sx * sy);
        R13 = 0.0 + sa * W13 + sa2 * (sx * sz);
        R21 = 0.0 + sa * W21 + sa2 * (sy * sx);
        R22 = 1.0 + sa * W22 + sa2 * (sy * sy - 1.0);
        R23 = 0.0 + sa * W23 + sa2 * (sy * sz);
        R31 = 0.0 + sa * W31 + sa2 * (sz * sx);
        R32 = 0.0 + sa * W32 + sa2 * (sz * sy);
        R33 = 1.0 + sa * W33 + sa2 * (sz * sz - 1.0);

        f0x = R11 * fx + R12 * fy + R13 * fz;
        f0y = R21 * fx + R22 * fy + R23 * fz;
        f0z = R31 * fx + R32 * fy + R33 * fz;
        normalize(f0x, f0y, f0z, 1.0, 0.0, 0.0);

        xfx = f0y * sz - f0z * sy;
        xfy = f0z * sx - f0x * sz;
        xfz = f0x * sy - f0y * sx;
        normalize(xfx, xfy, xfz, 0.0, 0.0, 1.0);

        sheets_system.solution->set(i, sx);
        sheets_system.solution->set(i + 1, sy);
        sheets_system.solution->set(i + 2, sz);
        xfiber_system.solution->set(i, xfx);
        xfiber_system.solution->set(i + 1, xfy);
        xfiber_system.solution->set(i + 2, xfz);
//        std::cout << "* MONODOMAIN: setting sol ... " << std::flush;

        fiber_system.solution->set(i, f0x);
        fiber_system.solution->set(i + 1, f0y);
        fiber_system.solution->set(i + 2, f0z);
//        std::cout << " done" << std::endl;

        i += 3;
    }

}

void MonodomainExplicit::amr(libMesh::MeshRefinement &mesh_refinement, const std::string &type)
{
//	std::cout << "* MONODOMAIN: starting AMR ... " << std::flush;
//	BeatIt:Timer timer;
//	std::cout << "Creating Error Vector " << std::endl;
//	timer.start();
    libMesh::ErrorVector error;
//	timer.stop();
//	timer.print(std::cout);
//	std::cout << "Creating Error Estimator " << std::endl;
//	timer.restart();
    libMesh::ErrorEstimator *p_error_estimator;
    if ("kelly" == type)
        p_error_estimator = new libMesh::KellyErrorEstimator;
    else if ("disc" == type)
        p_error_estimator = new libMesh::DiscontinuityMeasure;
    else
        p_error_estimator = new libMesh::LaplacianErrorEstimator;
//	timer.stop();
//	timer.print(std::cout);
//	 libMesh::KellyErrorEstimator error_estimator;
//	 libMesh::LaplacianErrorEstimator error_estimator;
    ElectroSystem &monodomain_system = M_equationSystems.get_system < ElectroSystem > (M_model);
    // WAVE
    ElectroSystem &wave_system = M_equationSystems.add_system < ElectroSystem > ("wave");

//	std::cout << "Estimating Monodomain Error  " << std::endl;
//	timer.restart();
    p_error_estimator->estimate_error(wave_system, error);
//	timer.stop();
//	timer.print(std::cout);
    // Flag elements to be refined and coarsened
//	std::cout << "Flagging Elements  " << std::endl;
//	timer.restart();
    mesh_refinement.flag_elements_by_error_fraction(error);
//	timer.stop();
//	timer.print(std::cout);
//
//	std::cout << "Refine and Coarsen  " << std::endl;
//	std::cout << " coarsen and refine ...  " << std::flush;
//	timer.restart();
    mesh_refinement.refine_and_coarsen_elements();
//	timer.stop();
//	timer.print(std::cout);
//	std::cout << " reinit ...  " << std::flush;
//	std::cout << "Reinit system  " << std::endl;
//	timer.restart();
    M_equationSystems.reinit();
//	timer.stop();
//	timer.print(std::cout);
//	timer.restart();
//	std::cout << " done  " << std::endl;
    delete p_error_estimator;
}

void MonodomainExplicit::assemble_matrices(double dt)
{
    if (M_FEFamily == libMesh::MONOMIAL || M_FEFamily == libMesh::L2_LAGRANGE)
        assemble_dg_matrices(dt);
    else
    {
        assemble_cg_matrices(dt);
    }

}
void MonodomainExplicit::assemble_cg_matrices(double dt)
{
    std::cout << "* MONODOMAIN: Assembling CG matrices ... " << std::endl;
    using std::unique_ptr;

    const libMesh::MeshBase &mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int max_dim = 3;
    const libMesh::Real Chi = M_equationSystems.parameters.get < libMesh::Real > ("Chi");

    // Get a reference to the LinearImplicitSystem we are solving
    ElectroSystem &monodomain_system = M_equationSystems.get_system < ElectroSystem > (M_model);
    //IonicModelSystem &ionic_model_system = M_equationSystems.get_system < IonicModelSystem > ("ionic_model");

    monodomain_system.get_matrix("mass").zero();
    monodomain_system.get_matrix("lumped_mass").zero();
    monodomain_system.get_matrix("high_order_mass").zero();
    monodomain_system.get_matrix("stiffness").zero();
    monodomain_system.get_vector("lumped_mass_vector").zero();

//     MatSetOption( (dynamic_cast<libMesh::PetscMatrix<libMesh::Number> >(monodomain_system.get_matrix("stiffness"))).mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//     MatSetOption( (dynamic_cast<libMesh::PetscMatrix<libMesh::Number> * >(monodomain_system.matrix))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    ParameterSystem &fiber_system = M_equationSystems.get_system < ParameterSystem > ("fibers");
    ParameterSystem &sheets_system = M_equationSystems.get_system < ParameterSystem > ("sheets");
    ParameterSystem &xfiber_system = M_equationSystems.get_system < ParameterSystem > ("xfibers");
    ParameterSystem &conductivity_system = M_equationSystems.get_system < ParameterSystem > ("conductivity");

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const libMesh::DofMap &dof_map_monodomain = monodomain_system.get_dof_map();
    const libMesh::DofMap &dof_map_fibers = fiber_system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    libMesh::FEType fe_type_qp1 = dof_map_monodomain.variable_type(0);
    libMesh::FEType fe_type_qp2 = dof_map_monodomain.variable_type(0);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr < libMesh::FEBase > fe_qp1(libMesh::FEBase::build(dim, fe_type_qp1));
    std::unique_ptr < libMesh::FEBase > fe_qp2(libMesh::FEBase::build(dim, fe_type_qp2));

    // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule_stiffness(dim, libMesh::SECOND);
    // A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule_mass(dim, libMesh::THIRD);

    // Tell the finite element object to use our quadrature rule.
    fe_qp1->attach_quadrature_rule(&qrule_stiffness);
    fe_qp2->attach_quadrature_rule(&qrule_mass);
    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<libMesh::Real> &JxW_qp1 = fe_qp1->get_JxW();
    const std::vector<libMesh::Real> &JxW_qp2 = fe_qp2->get_JxW();

    // The physical XY locations of the quadrature points on the element.
    // These might be useful for evaluating spatially varying material
    // properties at the quadrature points.
    const std::vector<libMesh::Point> &q_point_qp1 = fe_qp1->get_xyz();
    const std::vector<libMesh::Point> &q_point_qp2 = fe_qp2->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<libMesh::Real> > &phi_qp1 = fe_qp1->get_phi();
    const std::vector<std::vector<libMesh::Real> > &phi_qp2 = fe_qp2->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_qp1 = fe_qp1->get_dphi();
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_qp2 = fe_qp2->get_dphi();

    const std::vector<std::vector<libMesh::Real> > &dphidx_qp1 = fe_qp1->get_dphidx();
    const std::vector<std::vector<libMesh::Real> > &dphidy_qp1 = fe_qp1->get_dphidy();
    const std::vector<std::vector<libMesh::Real> > &dphidz_qp1 = fe_qp1->get_dphidz();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseMatrix<libMesh::Number> Mel;
    libMesh::DenseVector<libMesh::Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector < libMesh::dof_id_type > dof_indices;
    std::vector < libMesh::dof_id_type > dof_indices_fibers;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    //
    // Element iterators are a nice way to iterate through all the
    // elements, or all the elements that have some property.  The
    // iterator el will iterate from the first to// the last element on
    // the local processor.  The iterator end_el tells us when to stop.
    // It is smart to make this one const so that we don't accidentally
    // mess it up!  In case users later modify this program to include
    // refinement, we will be safe and will only consider the active
    // elements; hence we use a variant of the active_elem_iterator.
    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // Loop over the elements.  Note that  ++el is preferred to
    // el++ since the latter requires an unnecessary temporary
    // object.
    double f0[3];
    double s0[3];
    double n0[3];
    libMesh::RealGradient DgradV;
    libMesh::TensorValue<double> D0;

    /* for coduction block */
    // Declare a special finite element object for
    // boundary integration.
    std::unique_ptr < libMesh::FEBase > fe_face(libMesh::FEBase::build(dim, fe_type_qp2));
    libMesh::QGauss qface(dim - 1, libMesh::FOURTH);
    fe_face->attach_quadrature_rule(&qface);

    std::unique_ptr < libMesh::FEBase > fe_neighbor_face(libMesh::FEBase::build(dim, fe_type_qp1));
    fe_neighbor_face->attach_quadrature_rule(&qface);
    const std::vector<std::vector<libMesh::Real> > &phi_neighbor_face = fe_neighbor_face->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_neighbor_face = fe_neighbor_face->get_dphi();

    std::srand(std::time(0));
    std::cout << "* \t looping over elements ...  " << std::endl;

    for (; el != end_el; ++el)
    {
        const libMesh::Elem *elem = *el;
        const unsigned int elem_id = elem->id();
        dof_map_monodomain.dof_indices(elem, dof_indices);
        dof_map_fibers.dof_indices(elem, dof_indices_fibers);

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_qp1->reinit(elem);
        fe_qp2->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        auto n_dofs = dof_indices.size();
        Ke.resize(n_dofs, n_dofs);
        Me.resize(n_dofs, n_dofs);
        Mel.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);

        //        std::cout << "Fibers" << std::endl;
        // fiber direction
        try
        {
            dof_indices_fibers.at(2);
        }
        catch (std::exception &e)
        {
            std::cout << "dof indices fibers 2 not found: " << e.what() << std::endl;
        }


        switch (M_anisotropy)
        {
            case Anisotropy::Isotropic:
            {
                break;
            }
            case Anisotropy::TransverselyIsotropic:
            {
                std::cout << "f0, " << std::flush;
                f0[0] = (*fiber_system.solution)(dof_indices_fibers[0]);
                std::cout << "f1, " << std::flush;
                f0[1] = (*fiber_system.solution)(dof_indices_fibers[1]);
                std::cout << "f2, " << std::flush;
                f0[2] = (*fiber_system.solution)(dof_indices_fibers[2]);
                std::cout << "done.  " << std::flush;
                break;
            }
            case Anisotropy::Orthotropic:
            default:
            {
                f0[0] = (*fiber_system.solution)(dof_indices_fibers[0]);
                f0[1] = (*fiber_system.solution)(dof_indices_fibers[1]);
                f0[2] = (*fiber_system.solution)(dof_indices_fibers[2]);
                // sheet direction
                s0[0] = (*sheets_system.solution)(dof_indices_fibers[0]);
                s0[1] = (*sheets_system.solution)(dof_indices_fibers[1]);
                s0[2] = (*sheets_system.solution)(dof_indices_fibers[2]);
                // crossfiber direction
                n0[0] = (*xfiber_system.solution)(dof_indices_fibers[0]);
                n0[1] = (*xfiber_system.solution)(dof_indices_fibers[1]);
                n0[2] = (*xfiber_system.solution)(dof_indices_fibers[2]);
                break;
            }
        }

        // Conductivity tensor
        double Dff = (*conductivity_system.solution)(dof_indices_fibers[0]);
        double Dss = (*conductivity_system.solution)(dof_indices_fibers[1]);
        double Dnn = (*conductivity_system.solution)(dof_indices_fibers[2]);

        setup_local_conductivity(D0, Dff, Dss, Dnn, f0, s0, n0);
        D0 /= Chi;

//        int random_el = std::rand() % 100 + 1;
//        if(random_el <= 60)
//        {
//            D0 *= 0.0;
//        }

        // Assemble Mass terms
        for (unsigned int qp = 0; qp < qrule_mass.n_points(); qp++)
        {
            //  Matrix
            for (unsigned int i = 0; i < phi_qp2.size(); i++)
            {
                for (unsigned int j = 0; j < phi_qp2.size(); j++)
                {
                    // Mass term
                    Me(i, j) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                    Mel(i, i) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                    Fe(i) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                }
            }
        }
        monodomain_system.get_matrix("mass").add_matrix(Me, dof_indices);
        monodomain_system.get_matrix("lumped_mass").add_matrix(Mel, dof_indices);
        monodomain_system.get_vector("lumped_mass_vector").add_vector(Fe, dof_indices);
        for (unsigned int qp = 0; qp < qrule_stiffness.n_points(); qp++)
        {
            for (unsigned int i = 0; i < dphi_qp1.size(); i++)
            {
                DgradV = D0 * dphi_qp1[i][qp];

                for (unsigned int j = 0; j < dphi_qp1.size(); j++)
                {
                    // stiffness term
                    Ke(i, j) += JxW_qp1[qp] * DgradV * dphi_qp1[j][qp];
                }
            }
        }

//        int random_el = rand() % 100 + 1;
//        if(random_el <= 5)
//        {
//        // TEST CONDUCTION BLOCKS
//        for (unsigned int side = 0; side < elem->n_sides(); side++)
//        {
//            int random = rand() % 100 + 1;
//            std::cout << "random: " << random << std::endl;
//            if (random <= 2)
//            {
//                const std::vector<libMesh::Real> & JxW_face = fe_face->get_JxW();
//                const std::vector<std::vector<libMesh::Real> > & phi_face = fe_face->get_phi();
//                const std::vector<libMesh::Point> & qface_point = fe_face->get_xyz();
//                const std::vector<std::vector<libMesh::RealGradient> > & dphi_face = fe_face->get_dphi();
//                const std::vector<libMesh::Point>&  normals = fe_face->get_normals();
//                fe_face->reinit(elem, side);
//                // The location on the boundary of the current
//                // face quadrature point.
//
//                for (unsigned int qp=0; qp<qface.n_points(); qp++)
//                {
//
//                    const double xf = qface_point[qp](0);
//                    const double yf = qface_point[qp](1);
//                    const double zf = qface_point[qp](2);
//                   // if(yf > 0.3 && xf > 0.5)
//                {
//
//                double beta = 1e8;
//
//                    for (unsigned int i=0; i<phi_face.size(); i++)
//                    {
//                        for (unsigned int j=0; j<phi_face.size(); j++)
//                        {
//                            Ke(i,j) += JxW_face[qp] * beta * phi_face[i][qp] * phi_face[j][qp];
//                        }
//                    }
//                }
//                }
//            }
//        }
//            if (random <= 2) break;
//
//        }

        monodomain_system.get_matrix("stiffness").add_matrix(Ke, dof_indices);
    }
    std::cout << "* \t looping over elements done. Closing ...  " << std::endl;

    // closing matrices and vectors
    monodomain_system.get_matrix("mass").close();
    monodomain_system.get_matrix("lumped_mass").close();
    monodomain_system.get_matrix("stiffness").close();
    monodomain_system.get_vector("lumped_mass_vector").close();
    monodomain_system.get_matrix("high_order_mass").close();
    monodomain_system.get_matrix("high_order_mass").add(0.5, monodomain_system.get_matrix("mass"));
    monodomain_system.get_matrix("high_order_mass").add(0.5, monodomain_system.get_matrix("lumped_mass"));
    monodomain_system.get_matrix("high_order_mass").close();

    form_system_matrix(dt, false, "lumped_mass");
}

void MonodomainExplicit::setup_local_conductivity(libMesh::TensorValue<double> &D0, double Dff, double Dss, double Dnn, double *f0, double *s0, double *n0)
{
    const unsigned int max_dim = 3;
    switch (M_anisotropy)
    {
    case Anisotropy::Isotropic:
    {
        for (int idim = 0; idim < max_dim; ++idim)
        {
            D0(idim, idim) = Dff;
        }

        break;
    }

    case Anisotropy::TransverselyIsotropic:
    {
        for (int idim = 0; idim < max_dim; ++idim)
        {
            for (int jdim = 0; jdim < max_dim; ++jdim)
            {

                D0(idim, jdim) = (Dff - Dss) * f0[idim] * f0[jdim];
            }
            D0(idim, idim) += Dss;
        }
        break;
    }

    case Anisotropy::Orthotropic:
    default:
    {
        for (int idim = 0; idim < max_dim; ++idim)
        {
            for (int jdim = 0; jdim < max_dim; ++jdim)
            {
                D0(idim, jdim) = Dff * f0[idim] * f0[jdim] + Dss * s0[idim] * s0[jdim] + Dnn * n0[idim] * n0[jdim];

            }
        }
        break;
    }
    }
}

void MonodomainExplicit::assemble_dg_matrices(double dt)
{
    std::cout << "* MONODOMAIN: Assembling DG matrices ... " << std::endl;
    using std::unique_ptr;

    const libMesh::MeshBase &mesh = M_equationSystems.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const unsigned int max_dim = 3;
    const libMesh::Real Chi = M_equationSystems.parameters.get < libMesh::Real > ("Chi");

// Get a reference to the LinearImplicitSystem we are solving
    ElectroSystem &monodomain_system = M_equationSystems.get_system < ElectroSystem > (M_model);
    IonicModelSystem &ionic_model_system = M_equationSystems.get_system < IonicModelSystem > ("ionic_model");

    monodomain_system.get_matrix("mass").zero();
    monodomain_system.get_matrix("lumped_mass").zero();
    monodomain_system.get_matrix("high_order_mass").zero();
    monodomain_system.get_matrix("stiffness").zero();
    monodomain_system.get_vector("lumped_mass_vector").zero();

//     MatSetOption( (dynamic_cast<libMesh::PetscMatrix<libMesh::Number> >(monodomain_system.get_matrix("stiffness"))).mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//     MatSetOption( (dynamic_cast<libMesh::PetscMatrix<libMesh::Number> * >(monodomain_system.matrix))->mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    ParameterSystem &fiber_system = M_equationSystems.get_system < ParameterSystem > ("fibers");
    ParameterSystem &sheets_system = M_equationSystems.get_system < ParameterSystem > ("sheets");
    ParameterSystem &xfiber_system = M_equationSystems.get_system < ParameterSystem > ("xfibers");
    ParameterSystem &conductivity_system = M_equationSystems.get_system < ParameterSystem > ("conductivity");

// A reference to the  DofMap object for this system.  The  DofMap
// object handles the index translation from node and element numbers
// to degree of freedom numbers.  We will talk more about the  DofMap
// in future examples.
    const libMesh::DofMap &dof_map_monodomain = monodomain_system.get_dof_map();
    const libMesh::DofMap &dof_map_fibers = fiber_system.get_dof_map();

// Get a constant reference to the Finite Element type
// for the first (and only) variable in the system.
    libMesh::FEType fe_type_qp1 = dof_map_monodomain.variable_type(0);
    libMesh::FEType fe_type_qp2 = dof_map_monodomain.variable_type(0);

// Build a Finite Element object of the specified type.  Since the
// FEBase::build() member dynamically creates memory we will
// store the object as a std::unique_ptr<FEBase>.  This can be thought
// of as a pointer that will clean up after itself.  Introduction Example 4
// describes some advantages of  std::unique_ptr's in the context of
// quadrature rules.
    std::unique_ptr < libMesh::FEBase > fe_qp1(libMesh::FEBase::build(dim, fe_type_qp1));
    std::unique_ptr < libMesh::FEBase > fe_qp2(libMesh::FEBase::build(dim, fe_type_qp2));

// A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule_stiffness(dim, libMesh::SECOND);
// A 5th order Gauss quadrature rule for numerical integration.
    libMesh::QGauss qrule_mass(dim, libMesh::THIRD);

// Tell the finite element object to use our quadrature rule.
    fe_qp1->attach_quadrature_rule(&qrule_stiffness);
    fe_qp2->attach_quadrature_rule(&qrule_mass);
// Here we define some references to cell-specific data that
// will be used to assemble the linear system.
//
// The element Jacobian * quadrature weight at each integration point.
    const std::vector<libMesh::Real> &JxW_qp1 = fe_qp1->get_JxW();
    const std::vector<libMesh::Real> &JxW_qp2 = fe_qp2->get_JxW();

// The physical XY locations of the quadrature points on the element.
// These might be useful for evaluating spatially varying material
// properties at the quadrature points.
    const std::vector<libMesh::Point> &q_point_qp1 = fe_qp1->get_xyz();
    const std::vector<libMesh::Point> &q_point_qp2 = fe_qp2->get_xyz();

// The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<libMesh::Real> > &phi_qp1 = fe_qp1->get_phi();
    const std::vector<std::vector<libMesh::Real> > &phi_qp2 = fe_qp2->get_phi();

// The element shape function gradients evaluated at the quadrature
// points.
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_qp1 = fe_qp1->get_dphi();
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_qp2 = fe_qp2->get_dphi();

    const std::vector<std::vector<libMesh::Real> > &dphidx_qp1 = fe_qp1->get_dphidx();
    const std::vector<std::vector<libMesh::Real> > &dphidy_qp1 = fe_qp1->get_dphidy();
    const std::vector<std::vector<libMesh::Real> > &dphidz_qp1 = fe_qp1->get_dphidz();

// Define data structures to contain the element matrix
// and right-hand-side vector contribution.  Following
// basic finite element terminology we will denote these
// "Ke" and "Fe".  These datatypes are templated on
//  Number, which allows the same code to work for real
// or complex numbers.
    libMesh::DenseMatrix<libMesh::Number> Ke;
    libMesh::DenseMatrix<libMesh::Number> Me;
    libMesh::DenseMatrix<libMesh::Number> Mel;
    libMesh::DenseVector<libMesh::Number> Fe;

// for interior penalty
    std::unique_ptr < libMesh::FEBase > fe_elem_face(libMesh::FEBase::build(dim, fe_type_qp1));
    std::unique_ptr < libMesh::FEBase > fe_neighbor_face(libMesh::FEBase::build(dim, fe_type_qp1));
// Tell the finite element object to use our quadrature rule.
    libMesh::QGauss qface(dim - 1, fe_type_qp1.default_quadrature_order());

    fe_elem_face->attach_quadrature_rule(&qface);
    fe_neighbor_face->attach_quadrature_rule(&qface);
// Data for surface integrals on the element boundary
    const std::vector<std::vector<libMesh::Real> > &phi_face = fe_elem_face->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_face = fe_elem_face->get_dphi();
    const std::vector<libMesh::Real> &JxW_face = fe_elem_face->get_JxW();
    const std::vector<libMesh::Point> &qface_normals = fe_elem_face->get_normals();
    const std::vector<libMesh::Point> &qface_points = fe_elem_face->get_xyz();
// Data for surface integrals on the neighbor boundary
    const std::vector<std::vector<libMesh::Real> > &phi_neighbor_face = fe_neighbor_face->get_phi();
    const std::vector<std::vector<libMesh::RealGradient> > &dphi_neighbor_face = fe_neighbor_face->get_dphi();
// Data structures to contain the element and neighbor boundary matrix
// contribution. This matrices will do the coupling beetwen the dofs of
// the element and those of his neighbors.
// Ken: matrix coupling elem and neighbor dofs
    libMesh::DenseMatrix<libMesh::Number> Kne;
    libMesh::DenseMatrix<libMesh::Number> Ken;
    libMesh::DenseMatrix<libMesh::Number> Kee;
    libMesh::DenseMatrix<libMesh::Number> Knn;

    double deltaKn = 0.0;
    double deltaKn_neighobor = 0.0;

// This vector will hold the degree of freedom indices for
// the element.  These define where in the global system
// the element degrees of freedom get mapped.
    std::vector < libMesh::dof_id_type > dof_indices;
    std::vector < libMesh::dof_id_type > dof_indices_fibers;

// Now we will loop over all the elements in the mesh.
// We will compute the element matrix and right-hand-side
// contribution.
//
// Element iterators are a nice way to iterate through all the
// elements, or all the elements that have some property.  The
// iterator el will iterate from the first to// the last element on
// the local processor.  The iterator end_el tells us when to stop.
// It is smart to make this one const so that we don't accidentally
// mess it up!  In case users later modify this program to include
// refinement, we will be safe and will only consider the active
// elements; hence we use a variant of the active_elem_iterator.
    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

// Loop over the elements.  Note that  ++el is preferred to
// el++ since the latter requires an unnecessary temporary
// object.
    double f0[3];
    double s0[3];
    double n0[3];
    libMesh::RealGradient DgradV;
    libMesh::TensorValue<double> D0;
    libMesh::TensorValue<double> D0_neighbor;

    for (; el != end_el; ++el)
    {
        const libMesh::Elem *elem = *el;
        const unsigned int elem_id = elem->id();
        dof_map_monodomain.dof_indices(elem, dof_indices);
        dof_map_fibers.dof_indices(elem, dof_indices_fibers);

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe_qp1->reinit(elem);
        fe_qp2->reinit(elem);

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        auto n_dofs = dof_indices.size();
        Ke.resize(n_dofs, n_dofs);
        Me.resize(n_dofs, n_dofs);
        Mel.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);

        //        std::cout << "Fibers" << std::endl;
        // fiber direction
        f0[0] = (*fiber_system.solution)(dof_indices_fibers[0]);
        f0[1] = (*fiber_system.solution)(dof_indices_fibers[1]);
        f0[2] = (*fiber_system.solution)(dof_indices_fibers[2]);
        // sheet direction
        s0[0] = (*sheets_system.solution)(dof_indices_fibers[0]);
        s0[1] = (*sheets_system.solution)(dof_indices_fibers[1]);
        s0[2] = (*sheets_system.solution)(dof_indices_fibers[2]);
        // crossfiber direction
        n0[0] = (*xfiber_system.solution)(dof_indices_fibers[0]);
        n0[1] = (*xfiber_system.solution)(dof_indices_fibers[1]);
        n0[2] = (*xfiber_system.solution)(dof_indices_fibers[2]);
        // Conductivity tensor
        double Dff = (*conductivity_system.solution)(dof_indices_fibers[0]);
        double Dss = (*conductivity_system.solution)(dof_indices_fibers[1]);
        double Dnn = (*conductivity_system.solution)(dof_indices_fibers[2]);

        setup_local_conductivity(D0, Dff, Dss, Dnn, f0, s0, n0);
        D0 /= Chi;

        // Assemble Mass terms
        for (unsigned int qp = 0; qp < qrule_mass.n_points(); qp++)
        {
            //  Matrix
            for (unsigned int i = 0; i < phi_qp2.size(); i++)
            {
                for (unsigned int j = 0; j < phi_qp2.size(); j++)
                {
                    // Mass term
                    Me(i, j) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                    Mel(i, i) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                    Fe(i) += JxW_qp2[qp] * (phi_qp2[i][qp] * phi_qp2[j][qp]);
                }
            }
        }
        monodomain_system.get_matrix("mass").add_matrix(Me, dof_indices);
        monodomain_system.get_matrix("lumped_mass").add_matrix(Mel, dof_indices);
        monodomain_system.get_vector("lumped_mass_vector").add_vector(Fe, dof_indices);
        // Assemble stiffness matrix
        for (unsigned int qp = 0; qp < qrule_stiffness.n_points(); qp++)
        {
            for (unsigned int i = 0; i < dphi_qp1.size(); i++)
            {
                DgradV = D0 * dphi_qp1[i][qp];

                for (unsigned int j = 0; j < dphi_qp1.size(); j++)
                {
                    // stiffness term
                    Ke(i, j) += JxW_qp1[qp] * DgradV * dphi_qp1[j][qp];
                }
            }
        }
        monodomain_system.get_matrix("stiffness").add_matrix(Ke, dof_indices);

        // If the element has no neighbor on a side then that
        // side MUST live on a boundary of the domain.
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            double gamma = 10.0;
            if (elem->neighbor_ptr(side))
            {
                // Store a pointer to the neighbor we are currently
                // working on.
                const libMesh::Elem *neighbor = elem->neighbor_ptr(side);
                // Get the global id of the element and the neighbor
                const unsigned int neighbor_id = neighbor->id();

                // WARNING!!!! NOTE!!!!
                // Here I should use some check for amr:
                // check libmesh test: miscellaneous_ex5
                int random = rand() % 100 + 1;
                if (random >= 10)
                {
                    // Pointer to the element side
                    std::unique_ptr<const libMesh::Elem> elem_side(elem->build_side_ptr(side));

                    // h dimension to compute the interior penalty penalty parameter
                    const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
                    const unsigned int neighbor_b_order = static_cast<unsigned int>(fe_neighbor_face->get_order());
                    const double side_order = (elem_b_order + neighbor_b_order) / 2.;
                    // check against Ern, Stephansen, Zunino
                    //const double penalty = 0.5 * gamma * elem_side->volume() * elem_side->volume();
                    //const double h_elem = (elem->volume() / elem_side->volume()) / std::pow(side_order, 2.);
                    const double h_elem = elem->hmax();

                    // The quadrature point locations on the neighbor side
                    std::vector < libMesh::Point > qface_neighbor_point;

                    // The quadrature point locations on the element side
                    std::vector < libMesh::Point > qface_point;

                    // Reinitialize shape functions on the element side
                    fe_elem_face->reinit(elem, side);

                    // Get the physical locations of the element quadrature points
                    qface_point = fe_elem_face->get_xyz();

                    // Find their locations on the neighbor
                    unsigned int side_neighbor = neighbor->which_neighbor_am_i(elem);
                    libMesh::FEInterface::inverse_map(elem->dim(), fe_qp1->get_fe_type(), neighbor, qface_point, qface_neighbor_point);
                    // Calculate the neighbor element shape functions at those locations
                    fe_neighbor_face->reinit(neighbor, &qface_neighbor_point);
                    // Get the degree of freedom indices for the
                    // neighbor.  These define where in the global
                    // matrix this neighbor will contribute to.
                    std::vector < libMesh::dof_id_type > neighbor_dof_indices;
                    std::vector < libMesh::dof_id_type > neighbor_fiber_dof_indices;
                    dof_map_monodomain.dof_indices(neighbor, neighbor_dof_indices);
                    dof_map_fibers.dof_indices(neighbor, neighbor_fiber_dof_indices);

                    const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
                    // Zero the element and neighbor side matrix before
                    // summing them.  We use the resize member here because
                    // the number of degrees of freedom might have changed from
                    // the last element or neighbor.
                    // Note that Kne and Ken are not square matrices if neighbor
                    // and element have a different p level
                    Kne.resize(n_neighbor_dofs, n_dofs);
                    Ken.resize(n_dofs, n_neighbor_dofs);
                    Kee.resize(n_dofs, n_dofs);
                    Knn.resize(n_neighbor_dofs, n_neighbor_dofs);

                    //        std::cout << "Fibers" << std::endl;
                    // fiber direction
                    f0[0] = (*fiber_system.solution)(neighbor_fiber_dof_indices[0]);
                    f0[1] = (*fiber_system.solution)(neighbor_fiber_dof_indices[1]);
                    f0[2] = (*fiber_system.solution)(neighbor_fiber_dof_indices[2]);
                    // sheet direction
                    s0[0] = (*sheets_system.solution)(neighbor_fiber_dof_indices[0]);
                    s0[1] = (*sheets_system.solution)(neighbor_fiber_dof_indices[1]);
                    s0[2] = (*sheets_system.solution)(neighbor_fiber_dof_indices[2]);
                    // crossfiber direction
                    n0[0] = (*xfiber_system.solution)(neighbor_fiber_dof_indices[0]);
                    n0[1] = (*xfiber_system.solution)(neighbor_fiber_dof_indices[1]);
                    n0[2] = (*xfiber_system.solution)(neighbor_fiber_dof_indices[2]);
                    // Conductivity tensor
                    double Dff = (*conductivity_system.solution)(neighbor_fiber_dof_indices[0]);
                    double Dss = (*conductivity_system.solution)(neighbor_fiber_dof_indices[1]);
                    double Dnn = (*conductivity_system.solution)(neighbor_fiber_dof_indices[2]);                    //

                    setup_local_conductivity(D0_neighbor, Dff, Dss, Dnn, f0, s0, n0);
                    D0_neighbor /= Chi;

                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        double deltaKn = qface_normals[qp] * (D0 * qface_normals[qp]);
                        double deltaKn_neighbor = qface_normals[qp] * (D0_neighbor * qface_normals[qp]);

                        // e stands for elem
                        // n for neighbor
                        double we = deltaKn_neighbor / (deltaKn_neighbor + deltaKn);
                        double wn = deltaKn / (deltaKn_neighbor + deltaKn);
                        double gamma_K = we * deltaKn;
                        double alpha = 1.0;
                        double penalty = alpha * gamma_K / h_elem;
                        // Kee Matrix. Integrate the element test function i
                        // against the element test function j
                        for (unsigned int i = 0; i < n_dofs; i++)
                        {
                            for (unsigned int j = 0; j < n_dofs; j++)
                            {
                                // consistency
                                Kee(i, j) -= JxW_face[qp] * we * (qface_normals[qp] * (D0 * dphi_face[i][qp])) * phi_face[j][qp];
                                Kee(i, j) -= JxW_face[qp] * we * (qface_normals[qp] * (D0 * dphi_face[j][qp])) * phi_face[i][qp];
                                // stability
                                Kee(i, j) += JxW_face[qp] * penalty * phi_face[j][qp] * phi_face[i][qp];
                            }
                        }

                        // Knn Matrix. Integrate the neighbor test function i
                        // against the neighbor test function j
                        for (unsigned int i = 0; i < n_neighbor_dofs; i++)
                        {
                            for (unsigned int j = 0; j < n_neighbor_dofs; j++)
                            {
                                // consistency
                                //Kee(i, j) += penalty * JxW_face[qp] * (qface_normals[qp] * dphi_neighbor_face[i][qp])
                                //       * (qface_normals[qp] * dphi_neighbor_face[j][qp]);

//                                Knn(i,j) +=
//                                        penalty * JxW_face[qp] *
//                                  (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) +
//                                   phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

                                Knn(i, j) -= JxW_face[qp] * wn * (qface_normals[qp] * (D0_neighbor * dphi_neighbor_face[i][qp])) * phi_neighbor_face[j][qp];
                                Knn(i, j) -= JxW_face[qp] * wn * (qface_normals[qp] * (D0_neighbor * dphi_neighbor_face[j][qp])) * phi_neighbor_face[i][qp];
                                // stability
                                Knn(i, j) += JxW_face[qp] * penalty * phi_neighbor_face[j][qp] * phi_neighbor_face[i][qp];
                            }
                        }

                        // Kne Matrix. Integrate the neighbor test function i
                        // against the element test function j
                        for (unsigned int i = 0; i < n_neighbor_dofs; i++)
                        {
                            for (unsigned int j = 0; j < n_dofs; j++)
                            {
                                // consistency
                                Kne(i, j) -= JxW_face[qp] * wn * (qface_normals[qp] * (D0 * dphi_neighbor_face[i][qp])) * phi_face[j][qp];
                                Kne(i, j) -= JxW_face[qp] * we * (qface_normals[qp] * (D0 * dphi_face[j][qp])) * phi_neighbor_face[i][qp];
                                // consistency
//                                 Kne(i,j) +=
//                                         penalty * JxW_face[qp] *
//                                  (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) -
//                                   phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]));

                                // stability
                                Kne(i, j) -= JxW_face[qp] * penalty * phi_face[j][qp] * phi_neighbor_face[i][qp];
                            }
                        }

                        // Ken Matrix. Integrate the element test function i
                        // against the neighbor test function j
                        for (unsigned int i = 0; i < n_dofs; i++)
                        {
                            for (unsigned int j = 0; j < n_neighbor_dofs; j++)
                            {
                                // consistency
                                Ken(i, j) -= JxW_face[qp] * we * (qface_normals[qp] * (D0 * dphi_face[i][qp])) * phi_neighbor_face[j][qp];
                                Ken(i, j) -= JxW_face[qp] * wn * (qface_normals[qp] * (D0 * dphi_neighbor_face[j][qp])) * phi_face[i][qp];
//                                Ken(i,j) +=
//                                        penalty * JxW_face[qp] *
//                                  (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) -
//                                   phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));

                                // stability
                                Ken(i, j) -= JxW_face[qp] * penalty * phi_face[i][qp] * phi_neighbor_face[j][qp];
                            }
                        }
                    }

                    monodomain_system.get_matrix("stiffness").add_matrix(Kee, dof_indices);
                    monodomain_system.get_matrix("stiffness").add_matrix(Knn, neighbor_dof_indices);
                    monodomain_system.get_matrix("stiffness").add_matrix(Kne, neighbor_dof_indices, dof_indices);
                    monodomain_system.get_matrix("stiffness").add_matrix(Ken, dof_indices, neighbor_dof_indices);

                }

            }
            else // we are on the boundary
            {
                // Pointer to the element side
                std::unique_ptr<const libMesh::Elem> elem_side(elem->build_side_ptr(side));

                // h dimension to compute the interior penalty penalty parameter
                const unsigned int elem_b_order = static_cast<unsigned int>(fe_elem_face->get_order());
                const double side_order = elem_b_order;
                // check against Ern, Stephansen, Zunino
                //const double penalty = 0.5 * gamma * elem_side->volume() * elem_side->volume();
                const double h_elem = elem->hmax();
                //(elem->volume() / elem_side->volume()) / std::pow(side_order, 2.);

                // The quadrature point locations on the element side
                std::vector < libMesh::Point > qface_point;

                // Reinitialize shape functions on the element side
                fe_elem_face->reinit(elem, side);

                // Get the physical locations of the element quadrature points
                qface_point = fe_elem_face->get_xyz();

                // Zero the element and neighbor side matrix before
                // summing them.  We use the resize member here because
                // the number of degrees of freedom might have changed from
                // the last element or neighbor.
                // Note that Kne and Ken are not square matrices if neighbor
                // and element have a different p level
                Kee.resize(n_dofs, n_dofs);

                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                {
                    double deltaKn = qface_normals[qp] * (D0 * qface_normals[qp]);
                    double gamma_K = deltaKn;
                    double we = 1.0;
                    double alpha = 1.0;
                    double penalty = alpha * gamma_K / h_elem;
                    // Kee Matrix. Integrate the element test function i
                    // against the element test function j
                    for (unsigned int i = 0; i < n_dofs; i++)
                    {
                        for (unsigned int j = 0; j < n_dofs; j++)
                        {
                            // consistency
                            Kee(i, j) -= JxW_face[qp] * we * (qface_normals[qp] * (D0 * dphi_face[i][qp])) * phi_face[j][qp];
                            Kee(i, j) -= JxW_face[qp] * we * (qface_normals[qp] * (D0 * dphi_face[j][qp])) * phi_face[i][qp];
                            // stability
                            //Kee(i, j) += JxW_face[qp] * penalty * phi_face[j][qp] * phi_face[i][qp];
                        }
                    }
                }

                //monodomain_system.get_matrix("stiffness").add_matrix(Kee, dof_indices);
            }
        }

    }
// closing matrices and vectors
    monodomain_system.get_matrix("mass").close();
    monodomain_system.get_matrix("lumped_mass").close();
    monodomain_system.get_matrix("stiffness").close();
    monodomain_system.get_vector("lumped_mass_vector").close();
    monodomain_system.get_matrix("high_order_mass").close();
    monodomain_system.get_matrix("high_order_mass").add(0.5, monodomain_system.get_matrix("mass"));
    monodomain_system.get_matrix("high_order_mass").add(0.5, monodomain_system.get_matrix("lumped_mass"));

    form_system_matrix(dt, false, "lumped_mass");
}

void MonodomainExplicit::form_system_matrix(double dt, bool /*useMidpoint */, const std::string &mass)
{
    ElectroSystem &monodomain_system = M_equationSystems.get_system < ElectroSystem > (M_model);
// WAVE
    std::cout << "* MONODOMAIN EXP: forming system matrix using the " << mass << " matrix" << std::endl;
    M_systemMass = mass;
    double Cm = 1.0; // M_ionicModelPtr->membraneCapacitance();

    const libMesh::Real tau = M_equationSystems.parameters.get < libMesh::Real > ("tau");

    monodomain_system.matrix->zero();
    monodomain_system.matrix->close();

// Coefficient for matrix
// SBDF1
    double cdt = dt;
// SBDF2
    if (M_timestep_counter > 0 && M_timeIntegrator == TimeIntegrator::SecondOrderIMEX)
    {
        cdt = 2.0 / 3.0 * dt;
    }
// Matrix part
//    if(tau>0)
    {
        // TEST --------------------------------
        // S = Cm M
        monodomain_system.matrix->add(Cm , monodomain_system.get_matrix(mass));

        // //---------------------------------------
        // // S = Cm M Q^n+1 + tau / (c * dt) * Cm * M dQ + c dt K Q^n+1
        // monodomain_system.matrix->add(Cm * (1.0 + tau / (cdt)), monodomain_system.get_matrix(mass));
        // monodomain_system.matrix->add(cdt, monodomain_system.get_matrix("stiffness"));
    }
//    else
//    {
//        monodomain_system.matrix->add(Cm / cdt, monodomain_system.get_matrix(mass));
//        //monodomain_system.matrix->add(1.0, monodomain_system.get_matrix("stiffness"));
//    }
}

void MonodomainExplicit::form_system_rhs(double dt, bool useMidpoint, const std::string &mass)
{
    MonodomainSystem &monodomain_system = M_equationSystems.get_system < MonodomainSystem > (M_model);
// WAVE
    ElectroSystem &wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
    IonicModelSystem &iion_system = M_equationSystems.get_system < IonicModelSystem > ("iion");
    IonicModelSystem &istim_system = M_equationSystems.get_system < IonicModelSystem > ("istim");

    double Cm = 1.0; //M_ionicModelPtr->membraneCapacitance();
    const libMesh::Real tau = M_equationSystems.parameters.get < libMesh::Real > ("tau"); // time constant

    monodomain_system.rhs->zero();
    monodomain_system.rhs->close();

    auto &aux1 = monodomain_system.get_vector("aux1");
    aux1.zero();
    aux1.close();
    auto &aux2 = monodomain_system.get_vector("aux2");
    aux2.zero();
    aux2.close();

    auto &total_current = iion_system.get_vector("total_current");
    total_current.zero();
    total_current.close();
// RHS part:
// Evaluate

// commmon part is: -dt * Chi * M * I^n
//iion_system.solution->scale(dt*Chi);
    if (M_timestep_counter > 0 && TimeIntegrator::SecondOrderIMEX == M_timeIntegrator)
    {
        // SBDF2
        double cdt = 2.0 / 3.0 * dt;
        // RHS WAVE = Cm * tau / cdt * M * ( 4/3 Q^n - 1/3 Q^n-1 )
        //          - K * Z^n                            // Z^n = 4/3 V^n - 1/3V^n-1
        //          - M * ( 2*I^n-I^n-1) - tau * M * ( 2*dI^n-dI^n-1) - M * Istim

        // First compute -( 2*I^n - I^n-1+ 2*tau * dI^n - tau * dI^n-1 + Istim )
        // 2 * I^n
        total_current.add(-2.0, *iion_system.solution);
        // - I^n-1
        total_current.add(1.0, *iion_system.old_local_solution);
        // 2 * tau * dI^n
        total_current.add(-2.0 * tau, iion_system.get_vector("diion"));
        // - tau * dI^n-1
        total_current.add(tau, iion_system.get_vector("diion_old"));
        // Istim
        total_current.add(-1.0, *istim_system.solution);

        // Then compute M * ( I^n + tau * dI^n + Istim ) = M * total_current
        // adding it to the system RHS
        monodomain_system.get_matrix(mass).vector_mult_add(*monodomain_system.rhs, total_current);

        {
            // Second we compute the time derivative term
            // Cm * tau / cdt * M * ( 4/3 Q^n - 1/3 Q^n-1 )
            //
            // We store Cm * tau / cdt * M * ( 4/3 Q^n - 1/3 Q^n-1 ) in aux1
            aux1.add(4.0 / 3.0, *monodomain_system.old_local_solution);
            aux1.add(-1.0 / 3.0, *monodomain_system.older_local_solution);
            aux1.scale(Cm * tau / cdt);
            // Then we added to the RHS
            // M * aux1
            monodomain_system.get_matrix(M_systemMass).vector_mult_add(*monodomain_system.rhs, aux1);

            // Third compute - K * Z^n
            // We store -Z^n in aux2,  -Z^n = -4/3 V^n + 1/3V^n-1
            aux2.add(-4.0 / 3.0, *wave_system.old_local_solution);
            aux2.add(1.0 / 3.0, *wave_system.older_local_solution);
            // Then we added to the RHS
            // K * aux2
            monodomain_system.get_matrix("stiffness").vector_mult_add(*monodomain_system.rhs, aux2);
//            }
//            else
//            {
//                // We store Cm * tau / cdt * M * ( 4/3 Q^n - 1/3 Q^n-1 ) in aux1
//                aux1.add( 4.0/3.0, *monodomain_system.old_local_solution);
//                aux1.add(-1.0/3.0, *monodomain_system.older_local_solution);
//                aux1.scale(Cm/cdt);
//                // Then we added to the RHS
//                // M * aux1
//                monodomain_system.get_matrix(M_systemMass).vector_mult_add(*monodomain_system.rhs, aux1);
        }

    }
    else
    {
        // TEST----------------------------------------------
        // RHS = (1+dt*K)*V^n - dt*M*I^n - M*Istim
        double cdt = dt;
        // -I^n
        total_current.add(-1.0, *iion_system.solution);
        // -Istim
        total_current.add(-1.0, *istim_system.solution);
        // Then compute M * ( - I^n  + Istim ) = M * total_current
        monodomain_system.get_matrix(mass).vector_mult_add(*monodomain_system.rhs, total_current);
        //dt * V^n
        aux1.add(cdt,*monodomain_system.old_local_solution);
        // K * dt * V^n
        monodomain_system.get_matrix("stiffness").vector_mult_add(*monodomain_system.rhs, aux1);
        //+ V^n
        aux1.add(1.0,*monodomain_system.old_local_solution);
        //Sum it all to RHS
        *monodomain_system.rhs+= aux1;
    }
}

void MonodomainExplicit::solve_diffusion_step(double dt, double time, bool useMidpoint, const std::string &mass, bool reassemble)
{
// FORM RHS
    ElectroSystem &monodomain_system = M_equationSystems.get_system < ElectroSystem > (M_model);
//std::cout << "form_system_rhs" << std::endl;
    form_system_rhs(dt, useMidpoint, mass);
//std::cout << "form_system_rhs done" << std::endl;
    const libMesh::Real tau = M_equationSystems.parameters.get < libMesh::Real > ("tau"); // time constant

// If we are using SBDF2, we need to compute the system matrix again
// In fact we do a first step with Forward-Backward Euler
// and then we proceed with SBDF2
// std::cout << "call from system matrix? " << M_timestep_counter << std::endl;
    if (M_timestep_counter == 1 && TimeIntegrator::SecondOrderIMEX == M_timeIntegrator)
    {
        std::cout << "Yes! call from system matrix: " << M_timestep_counter << std::endl;
        form_system_matrix(dt, false, M_systemMass);
    }
//++M_timestep_counter;

    double tol = 1e-12;
    double max_iter = 2000;

    std::pair<unsigned int, double> rval = std::make_pair(0, 0.0);

//std::cout << "Solving" << std::endl;
//
//monodomain_system.matrix->print();
//monodomain_system.rhs->print();
    rval = M_linearSolver->solve(*monodomain_system.matrix, *monodomain_system.solution, *monodomain_system.rhs, tol, max_iter);

// std::cout << "solve done" << std::endl;
// WAVE
    ElectroSystem &wave_system = M_equationSystems.get_system < ElectroSystem > ("wave");
    {
        if (0 < M_timestep_counter && TimeIntegrator::SecondOrderIMEX == M_timeIntegrator)
        {
            // Use BDF2
            // V^n+1 = 4/3 V^n - 1/3 V^n + 2/3 dt * Q^n+1
            *wave_system.solution = *monodomain_system.solution;
            wave_system.solution->scale(2.0 / 3.0 * dt);
            wave_system.solution->add(4.0 / 3.0, *wave_system.old_local_solution);
            wave_system.solution->add(-1.0 / 3.0, *wave_system.older_local_solution);
        }

        else
        {
            // *wave_system.solution = *monodomain_system.solution;
            // wave_system.solution->scale(dt);
            // *wave_system.solution += *wave_system.old_local_solution;
            
            *wave_system.solution = *monodomain_system.solution;
            wave_system.solution->scale(dt);
            *wave_system.solution += *wave_system.old_local_solution;

        }
    }

    M_timestep_counter++;
}

} /* namespace BeatIt */
