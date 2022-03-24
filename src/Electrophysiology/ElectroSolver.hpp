/*
 * ElectroSolver.hpp
 *
 *  Created on: Oct 19, 2017
 *      Author: srossi
 */

#ifndef SRC_ELECTROPHYSIOLOGY_ELECTROSOLVER_HPP_
#define SRC_ELECTROPHYSIOLOGY_ELECTROSOLVER_HPP_


#include "Util/TimeData.hpp"
// Include file that defines (possibly multiple) systems of equations.
#include "libmesh/equation_systems.h"
#include "libmesh/linear_solver.h"
#include "libmesh/petsc_linear_solver.h"

#include <memory>
#include "Util/SpiritFunction.hpp"
#include "Util/Enums.hpp"
#include "Util/Factory.hpp"
#include "Util/Timer.hpp"
#include "libmesh/id_types.h"
#include "BoundaryConditions/BCHandler.hpp"

// Forward Definition
namespace libMesh
{
class Mesh;
class ExodusII_IO;
class Nemesis_IO;
class VTKIO;
class GMVIO;
class MeshRefinement;
class TimeData;
class ExplicitSystem;
class LinearImplicitSystem;
class PacingProtocol;
class ErrorVector;
class  MeshRefinement;
class BoundaryMesh;
}

namespace BeatIt
{

/*!
 *  forward declarations
 */
class BCHandler;
class IonicModel;
class PacingProtocol;

enum class Anisotropy { Isotropic,
                        TransverselyIsotropic,
                        Orthotropic,
                        UserDefined };

enum class EquationType { ReactionDiffusion,
                          Wave,
                          ParabolicEllipticBidomain,
                          ParabolicEllipticHyperbolic,
                          ParabolicParabolicHyperbolic   };
enum class ModelType { Monodomain, Bidomain, BidomainWithBath };
enum class TimeIntegrator { FirstOrderIMEX,     // FORWARD-BACKWARD EULER
                            SecondOrderIMEX  }; // SBDF2

enum class Ground { Nullspace,
                    GroundNode,
                    Dirichlet };


class ElectroSolver
{
public:
    typedef FactoryArg<ElectroSolver, std::string, libMesh::EquationSystems >     ElectroFactory;
    typedef libMesh::VTKIO Exporter;
    typedef libMesh::Nemesis_IO NemesisIO;
    // Another alternative when not using AMR
    typedef libMesh::ExodusII_IO EXOExporter;

    ElectroSolver( libMesh::EquationSystems& es, std::string model = "bidomain" );
    virtual ~ElectroSolver();


    void setup(GetPot& data, std::string section);
    virtual void setup_systems(GetPot& data, std::string section) = 0;
    void setup_ODE_systems(GetPot& data, std::string section);

    void restart( EXOExporter& importer, int step = 0, bool restart = true );
    void read_fibers( EXOExporter& importer, int step = 1);

    void init(double time);
    void init_systems(double time);
    void save(int step);
    void save_exo_timestep(int step, double time);
    void save_ve_timestep(int step, double time);
    void init_exo_output();
    void save_potential(int step, double time = 0.0);
    void save_potential_nemesis(int step, double time = 0.0);
    void save_parameters();
    void save_activation_times(int step = 1);
    void save_activation_times_nemesis(int step = 1);
    void save_conduction_velocity(int step = 1);

    virtual void amr( libMesh:: MeshRefinement& mesh_refinement, const std::string& type = "kelly" ) {}
    void reinit_linear_solver();
    //void update_pacing(double time);
    void update_activation_time(double time, double threshold = 0.8);
    void evaluate_conduction_velocity();



    //void cut(double time, std::string f);
    virtual void assemble_matrices(double dt = 1.0)  = 0;
    virtual void form_system_matrix(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass") = 0;
    virtual void form_system_rhs(double dt, bool useMidpoint = true, const std::string& mass = "lumped_mass") = 0;
    void advance();
    virtual void solve_reaction_step( double dt,
                              double time,
                              int step = 0,
                              bool useMidpoint = true,
                              const std::string& mass = "mass",
                              libMesh::NumericVector<libMesh::Number>* I4f_ptr = nullptr);

    virtual void solve_reaction_step_cg( double dt,
                              double time,
                              int step = 0,
                              bool useMidpoint = true,
                              const std::string& mass = "mass",
                              libMesh::NumericVector<libMesh::Number>* I4f_ptr = nullptr);

    virtual void solve_reaction_step_dg( double dt,
                              double time,
                              int step = 0,
                              bool useMidpoint = true,
                              const std::string& mass = "mass",
                              libMesh::NumericVector<libMesh::Number>* I4f_ptr = nullptr);

    virtual void solve_diffusion_step(double dt, double time,  bool useMidpoint = true, const std::string& mass = "lumped_mass", bool reassemble = true) = 0;
    virtual void generate_fibers(   const GetPot& data,
                            const std::string& section = "rule_based_fibers" ) {}
   double last_activation_time();
   double potential_norm();

   void set_potential_on_boundary(unsigned int boundID, double value = 1.0, int subdomain = -1);
   void setup_ic(libMesh::FunctionBase<libMesh::Number>& ic, double time=0.0); // setup initial conditions

   std::string  get_ionic_model_name(unsigned int key) const;


   const std::unique_ptr<libMesh::NumericVector<libMesh::Number> >&
   get_fibers();
   const std::unique_ptr<libMesh::NumericVector<libMesh::Number> >&
   get_sheets();
   const std::unique_ptr<libMesh::NumericVector<libMesh::Number> >&
   get_xfibers();
    //protected:


    /// input file
    GetPot                     M_datafile;
    /// Store pointer to the ionic model
//    std::unique_ptr<IonicModel> M_ionicModelPtr;
    std::map<unsigned int, std::shared_ptr<IonicModel> > M_ionicModelPtrMap;
    std::map<unsigned int, std::string > M_ionicModelNameMap;
    /// Equation Systems: One for the potential and one for the other variables
    /*!
     *  Use separate systems to avoid saving in all the variables
     */
    libMesh::EquationSystems&  M_equationSystems;

    std::unique_ptr<Exporter> M_exporter;
    std::unique_ptr<NemesisIO> M_nemesis_exporter;
    std::set<std::string> M_exporterNames;
    std::unique_ptr<Exporter> M_ionicModelExporter;
    std::set<std::string> M_ionicModelExporterNames;
    std::unique_ptr<EXOExporter> M_parametersExporter;
    std::set<std::string> M_parametersExporterNames;
    std::unique_ptr<EXOExporter> M_EXOExporter;
    std::unique_ptr<EXOExporter> M_potentialEXOExporter;

    std::string  M_outputFolder;
    bool M_assembleMatrix;
    bool M_useAMR;
    std::string  M_systemMass;

    std::unique_ptr<PacingProtocol> M_pacing; // intracellular
    std::unique_ptr<PacingProtocol> M_pacing_i; // intracellular
    std::unique_ptr<PacingProtocol> M_pacing_e; // extracellular
    std::unique_ptr<PacingProtocol> M_surf_pacing_i; // intracellular
    std::unique_ptr<PacingProtocol> M_surf_pacing_e; // extracellular


    std::unique_ptr<libMesh::PetscLinearSolver<libMesh::Number> > M_linearSolver;
    std::vector<double>  M_intraConductivity;
    std::vector<double>  M_extraConductivity;

    //int M_tissueBlockID;
    std::set<short unsigned int> M_tissueBlockIDs;
    int secret_blockID_key;
    std::vector < std::string > M_ionic_models_systems_name_vec;
    std::vector < std::string > M_ionic_models_vec;
    std::vector < unsigned int > M_ionic_models_IDs_vec;

    Anisotropy M_anisotropy;
    std::vector<double>  M_conductivity;
    EquationType M_equationType;
    DynamicTimeIntegratorType M_timeIntegratorType;
    TimeIntegrator M_timeIntegrator;
protected:
    long int M_timestep_counter;
    bool M_symmetricOperator;

public:

    long int& timestep_counter() { return M_timestep_counter; };
    std::string M_model;
    std::string model() const {return M_model;}

    double M_meshSize;
    ModelType M_modelType;
    std::string M_section;

    std::set<unsigned int> M_transmembranePotentialActiveSubdomains;
    int M_constraint_dof_id;
    Ground M_ground_ve;
    libMesh::Order M_order;
    libMesh::FEFamily M_FEFamily;
    BCHandler M_bch;

    struct EndocardialVe
    {
        std::unique_ptr<libMesh::BoundaryMesh> M_endocardium;
        std::unique_ptr<libMesh::EquationSystems>  M_boundary_es;
        std::map<libMesh::dof_id_type, libMesh::dof_id_type> M_reverse_node_id_map;
        std::map<libMesh::dof_id_type, libMesh::dof_id_type> M_node_id_map;

        std::unique_ptr<EXOExporter> M_EXOExporter;
    };
    EndocardialVe M_boundary_ve;
    void init_endocardial_ve(std::set<libMesh::boundary_id_type>& IDs, std::set<unsigned short>& subdomainIDs);

    Timer::duration_Type M_elapsed_time;
    unsigned int M_num_linear_iters;


};

} /* namespace BeatIt */

#endif /* SRC_ELECTROPHYSIOLOGY_ELECTROSOLVER_HPP_ */
