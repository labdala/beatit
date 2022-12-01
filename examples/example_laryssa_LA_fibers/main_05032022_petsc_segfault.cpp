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
#include "Electrophysiology/Bidomain/BidomainWithBath.hpp"
#include "Electrophysiology/Monodomain/MonodomainUtil.hpp"

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"

#include "Util/IO/io.hpp"
#include "libmesh/mesh.h"
#include "libmesh/elem.h"
//#include "timpi/communicator.h"
#include "libmesh/communicator.h"


#include <random>


// Packages used to create the fibers

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// Use input file
#include "libmesh/getpot.h"

// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/vtk_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"

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
#include "libmesh/enum_solver_package.h"

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Include headerfile to detect boundary id
#include "libmesh/boundary_info.h"

// To read the IDs from the input file
#include "Util/IO/io.hpp"

// Bring in everything from the libMesh namespace
using namespace libMesh;


// Subroutines to create the fibers

class AnatomicalParameters;

class Poisson{
	public:
		// Member variable Dirichlet BC list - can use inside poisson_assemble
		// each list stores subregion IDs to where that boundary condition should be imposed
		std::vector<int> dirichlet_side_set_list;
		std::vector<double> dirichlet_bc_list;


		// Stores a list of the names of all input files for each problem we are solving
		std::vector<std::string> input_list;

		// For each input file i we need:
        std::string dirichlet_side_set_list_aux;
        std::string dirichlet_bc_list_aux;
        std::string dirichlet_x_list;
        std::string dirichlet_y_list;
        std::string dirichlet_z_list;
        std::string dirichlet_r_list;


		// Pointer to a vector of points x,y,z
		std::vector<libMesh::Point> dirichlet_id_points;

		// Radius of sphere (neighborhood) where we impose Dirichlet BC
		std::vector<double> dirichlet_id_radius;

		//Reinitializes all the member variables
		void reinit(std::string, std::string,std::string,
				 std::string, std::string, std::string);

		//reads input secondary input file i
		void read_input_i(double, EquationSystems &);

		// assemble_poisson
		void assemble_poisson(EquationSystems&, const std::string&, int i, AnatomicalParameters &);

};

void Poisson::reinit( std::string dirichlet_x_list, std::string dirichlet_y_list, std::string dirichlet_z_list,
					 std::string dirichlet_r_list, std::string dirichlet_side_set_list_aux, std::string dirichlet_bc_list_aux)
{

    //erase what was stored in member variables
    dirichlet_side_set_list.clear();
    dirichlet_bc_list.clear();

    dirichlet_id_points.clear();

    dirichlet_id_radius.clear();

    //Read the BC list for the current problem
    BeatIt::readList(dirichlet_side_set_list_aux, dirichlet_side_set_list);
    BeatIt::readList(dirichlet_bc_list_aux, dirichlet_bc_list);

    //for general
    std::vector<double>  x, y, z, r;
    BeatIt::readList(dirichlet_x_list ,     x);
    BeatIt::readList(dirichlet_y_list ,     y);
    BeatIt::readList(dirichlet_z_list ,     z);
    BeatIt::readList(dirichlet_r_list ,     r);

    for (int jj = 0; jj < r.size(); jj++) {
        dirichlet_id_points.emplace_back(x[jj], y[jj], z[jj]);
        dirichlet_id_radius.push_back(r[jj]);
    }

    //check for error
    if (x.size()!=y.size() || x.size()!=z.size()|| x.size()!=r.size())
    {
       	std::cout << "the sizes of x, y, z, and r lists should be the same "<< std::endl;
       	std::cout<< "x.size() = " << x.size()  <<std::endl;
       	std::cout<< "y.size() = " << y.size()  <<std::endl;
       	std::cout<< "z.size() = " << z.size()  <<std::endl;
       	std::cout<< "r.size() = " << r.size()  <<std::endl;
       	throw std::runtime_error("the sizes of x, y, z, and r lists should be the same");
    }

    if(dirichlet_bc_list.size() != dirichlet_side_set_list.size() || dirichlet_bc_list.size() != dirichlet_id_radius.size()){
    	std::cout << "Size of Idirichlet_bc_list and Idirichlet_side_set_list have to be equal"<<std::endl;
    	std::cout << "dirichlet_bc_list.size() = " << dirichlet_bc_list.size() << ", dirichlet_side_set_list.size() = " <<
    			dirichlet_side_set_list.size() << " , dirichlet_id_radius.size() = " << dirichlet_id_radius.size()<<std::endl;
    	throw std::runtime_error("Size of Idirichlet_bc_list and Idirichlet_side_set_list has to be equal");
    }

    //Check
    if(dirichlet_id_points.size()!= dirichlet_id_radius.size()){
            	std::cout << " The number of points does not match the number of radii" << std::endl;
            	std::cout << "the size of vector points is "<< dirichlet_id_points.size()<<std::endl;
            	std::cout << "the size of vector radius is "<< dirichlet_id_radius.size()<<std::endl;
            	throw std::runtime_error(" The number of points does not match the number of radii"); // throw an exception
    }
	   //Print info on screen
     std::cout << "Setting Dirichlet BCs on " << dirichlet_id_points.size() << " points" << std::endl;

}

void Poisson::read_input_i( double i, EquationSystems & equation_systems)
{
  	std::cout << "----------------------------------------\n";
	std::cout << "----------------------------------------\n";

    //name of the problem
    GetPot data_i(input_list[i]);

    dirichlet_side_set_list_aux.clear();
    dirichlet_bc_list_aux.clear();
    dirichlet_x_list.clear();
    dirichlet_y_list.clear();
    dirichlet_z_list.clear();
    dirichlet_r_list.clear();

    // Read Dirichlet BC ID list from input file
    dirichlet_side_set_list_aux = data_i("Idirichlet_side_set_list", " ");
    dirichlet_bc_list_aux       = data_i("Idirichlet_bc_list"      , " ");

    // Read Dirichlet BC ID list for septum and laa points
    dirichlet_x_list = data_i("Ix", " ");
    dirichlet_y_list = data_i("Iy", " ");
    dirichlet_z_list = data_i("Iz", " ");
    dirichlet_r_list = data_i("Ir", " ");

   std::cout << "Ix, Iy, Iz and Ir are respectively: "<< dirichlet_x_list
		     << " | " << dirichlet_y_list << " | "<< dirichlet_z_list <<
			    " | " << dirichlet_r_list << std::endl;

   std::cout << "Idirichlet_side_set_list = " <<dirichlet_side_set_list_aux <<std::endl;
   std::cout << "Idirichlet_bc_list = " <<dirichlet_bc_list_aux <<std::endl;

}

// Save the fibers
void save_fibers(std::string e_test, std::string e_output_file_name, std::string e_output_number,
							EquationSystems & equation_systems, ParallelMesh & mesh){
	    std::string output_path = e_test + e_output_file_name + e_output_number + ".pvtu";
	    VTKIO(mesh).write_equation_systems(output_path, equation_systems);

	    // Export fibers in nemesis format
	    ExodusII_IO nemesis_exporter(mesh); // creates the exporter
	    std::vector<std::string> output_variables(4); // creates a vector that stores the names of vectors
	    /*output_variables[0] = "fibersx";
	    output_variables[1] = "fibersy";
	    output_variables[2] = "fibersz";
	    output_variables[3] = "sheetsx";
	    output_variables[4] = "sheetsy";
	    output_variables[5] = "sheetsz";
	    output_variables[6] = "xfibersx";
	    output_variables[7] = "xfibersy";
	    output_variables[8] = "xfibersz";
	    output_variables[9]  = "u0";
	    output_variables[10] = "u1";
	    output_variables[11] = "u2";
	    output_variables[12] = "u3";*/
	    output_variables[0] = "sheetsx";
	    output_variables[1] = "sheetsy";
	    output_variables[2] = "sheetsz";
	    output_variables[3] = "u0";
	    mesh.get_boundary_info().sideset_name(4)="mvr";
        nemesis_exporter.set_output_variables(output_variables);
	    nemesis_exporter.write_equation_systems(e_test + e_output_file_name + e_output_number + ".e", equation_systems); // saves the variables at the nodes
	    nemesis_exporter.write_element_data(equation_systems); // this saves the variables at the centroid of the elements
}

// The AnatomicalParameters conntains all parameters that have to do with
// region selection and building the fiber
class AnatomicalParameters{
	public:

	// this info is read from the main input file (input.main)
	// i_region is used to define the anatomical regions
	// i_epic_region is used to define the fibers in the anatomical region on the epicardium
	// i_endo_region is used to define the fibers in the anatomical region on the endocardium

	double milano; // could be boolean

    // radius sphere LAA tip
    double radius_LAA_tip;
	libMesh::Gradient f0, s0, n0;

	//Reinitialize variables
	void reinit(std::string, std::string);

	// Define regions
	int define_regions (std::vector<double> &);
	// function to define fibers  | input: u, output: fiber field??
	void define_fibers(std::vector<double> & , std::vector<libMesh::Gradient> & , int, const auto&);
};

void AnatomicalParameters::reinit( std::string milano_aux, std::string radius_laa_tip_aux)
	{
	std::vector<double> milano_vec;
	std::vector<double> radius_LAA_tip_vec;

	BeatIt::readList(milano_aux, milano_vec);
	milano = milano_vec[0];

	BeatIt::readList(radius_laa_tip_aux, radius_LAA_tip_vec);
	radius_LAA_tip = radius_LAA_tip_vec[0];

}

int AnatomicalParameters::define_regions (std::vector<double> & u){
	int blockid;
	//MILANO
	if(milano==1){
		if ((u[3]>=0.7)){     //MV
			blockid = 1;
		}
		else if ((u[2]>=0.85) || (u[2]<=0.1)) { // PVs
			blockid = 2;
		}
		else
			blockid = 3;
	}
	//OUR MODEL
	else{
		throw std::runtime_error("NOT MILANO!");
	}
    return blockid;
}

void AnatomicalParameters::define_fibers(std::vector<double> & u, std::vector<libMesh::Gradient> & du, int block_id, const auto& elem){

    // Solve the transmural problem as the first one, such that
    // the next line is always true
	s0 = du[0].unit();
/*
    //	FIBERS FIRST WAY - TRADITIONAL - LIKE MILANO PAPER 2021
	if(milano==1){
		switch (block_id)
		{
		//
		case 1:
		{
				n0 = du[3].unit();
				f0 = s0.cross(n0);
				break;
		}
		//
		case 2:
		{
				n0 = du[2].unit();
				f0 = s0.cross(n0);
				break;
		}
		//
		case 3:
		{
				n0 = du[1].unit();
				f0 = s0.cross(n0);
				break;
		}
		default:
		{
			std::cout << "Case default, element = "<< elem <<"\n";
			f0(0) = 1.0; f0(1) = 0.0; f0(2) = 0.0;
			s0(0) = 0.0; s0(1) = 1.0; s0(2) = 0.0;
			n0(0) = 0.0; n0(1) = 0.0; n0(2) = 1.0;
			break;
		}
		}
    }
	else{
		throw std::runtime_error("NOT MILANO!");
		}

f0 = f0.unit(); */
s0 = s0.unit();
//n0 = n0.unit();
}


// We now define the matrix assembly function for the
// Poisson system.  We need to first compute element
// matrices and right-hand sides, and then take into
// account the boundary conditions, which will be handled
// via a penalty method.
void Poisson::assemble_poisson(EquationSystems& es,
    const std::string& system_name, int i_solution, AnatomicalParameters & anatomic_parameters)
{

    // It is a good idea to make sure we are assembling
    // the proper system.
    libmesh_assert_equal_to(system_name, "Poisson");

    // Get a constant reference to the mesh object.
    const MeshBase& mesh = es.get_mesh();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the LinearImplicitSystem we are solving
    LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("Poisson");

    //Zero out the solution
    system.solution->zero();

    //Zero out the rhs
    system.rhs->zero();

    //Zero out the solution
    system.matrix->zero();

    // A reference to the  DofMap object for this system.  The  DofMap
    // object handles the index translation from node and element numbers
    // to degree of freedom numbers.  We will talk more about the  DofMap
    // in future examples.
    const DofMap& dof_map = system.get_dof_map();

    // Get a constant reference to the Finite Element type
    // for the first (and only) variable in the system.
    FEType fe_type = dof_map.variable_type(0);

    // Build a Finite Element object of the specified type.  Since the
    // FEBase::build() member dynamically creates memory we will
    // store the object as a std::unique_ptr<FEBase>.  This can be thought
    // of as a pointer that will clean up after itself.  Introduction Example 4
    // describes some advantages of  std::unique_ptr's in the context of
    // quadrature rules.
    std::unique_ptr<FEBase> fe(FEBase::build(dim, fe_type));

    // A 5th order Gauss quadrature rule for numerical integration.
    QGauss qrule(dim, FIFTH);

    // Tell the finite element object to use our quadrature rule.
    fe->attach_quadrature_rule(&qrule);

    // Declare a special finite element object for
    // boundary integration.
    std::unique_ptr<FEBase> fe_face(FEBase::build(dim, fe_type));

    // Boundary integration requires one quadrature rule,
    // with dimensionality one less than the dimensionality
    // of the element.
    QGauss qface(dim - 1, FIFTH);

    // Tell the finite element object to use our
    // quadrature rule.
    fe_face->attach_quadrature_rule(&qface);

    // Here we define some references to cell-specific data that
    // will be used to assemble the linear system.
    //
    // The element Jacobian * quadrature weight at each integration point.
    const std::vector<Real>& JxW = fe->get_JxW();

    // The physical XY locations of the quadrature points on the element.
    // These might be useful for evaluating spatially varying material
    // properties at the quadrature points.
    const std::vector<Point>& q_point = fe->get_xyz();

    // The element shape functions evaluated at the quadrature points.
    const std::vector<std::vector<Real>>& phi = fe->get_phi();

    // The element shape function gradients evaluated at the quadrature
    // points.
    const std::vector<std::vector<RealGradient>>& dphi = fe->get_dphi();

    // Define data structures to contain the element matrix
    // and right-hand-side vector contribution.  Following
    // basic finite element terminology we will denote these
    // "Ke" and "Fe".  These datatypes are templated on
    //  Number, which allows the same code to work for real
    // or complex numbers.
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    // This vector will hold the degree of freedom indices for
    // the element.  These define where in the global system
    // the element degrees of freedom get mapped.
    std::vector<dof_id_type> dof_indices;

    // Now we will loop over all the elements in the mesh.
    // We will compute the element matrix and right-hand-side
    // contribution.
    //
    // Element iterators are a nice way to iterate through all the
    // elements, or all the elements that have some property.  The
    // iterator el will iterate from the first to the last element on
    // the local processor.  The iterator end_el tells us when to stop.
    // It is smart to make this one const so that we don't accidentally
    // mess it up!  In case users later modify this program to include
    // refinement, we will be safe and will only consider the active
    // elements; hence we use a variant of the active_elem_iterator.
    for (const auto& elem : mesh.active_local_element_ptr_range())
    {
        // Get the degree of freedom indices for the
        // current element.  These define where in the global
        // matrix and right-hand-side this element will
        // contribute to.
        dof_map.dof_indices(elem, dof_indices);

        // Cache the number of degrees of freedom on this element, for
        // use as a loop bound later.  We use cast_int to explicitly
        // convert from size() (which may be 64-bit) to unsigned int
        // (which may be 32-bit but which is definitely enough to count
        // *local* degrees of freedom.
        const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

        // Compute the element-specific data for the current
        // element.  This involves computing the location of the
        // quadrature points (q_point) and the shape functions
        // (phi, dphi) for the current element.
        fe->reinit(elem);

        // With one variable, we should have the same number of degrees
        // of freedom as shape functions.
        libmesh_assert_equal_to(n_dofs, phi.size());

        // Zero the element matrix and right-hand side before
        // summing them.  We use the resize member here because
        // the number of degrees of freedom might have changed from
        // the last element.  Note that this will be the case if the
        // element type is different (i.e. the last element was a
        // triangle, now we are on a quadrilateral).

        // The  DenseMatrix::resize() and the  DenseVector::resize()
        // members will automatically zero out the matrix  and vector.
        Ke.resize(n_dofs, n_dofs);

        Fe.resize(n_dofs);

        // Now loop over the quadrature points.  This handles
        // the numeric integration.
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {

            // Now we will build the element matrix.  This involves
            // a double loop to integrate the test functions (i) against
            // the trial functions (j).
            for (unsigned int i = 0; i != n_dofs; i++)
                for (unsigned int j = 0; j != n_dofs; j++)
                {
                    Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                }

            // This is the end of the matrix summation loop
            // Now we build the element right-hand-side contribution.
            // This involves a single loop in which we integrate the
            // "forcing function" in the PDE against the test functions.
            {
                const Real x = q_point[qp](0);
                const Real y = q_point[qp](1);
                const Real eps = 1.e-3;


                // "fxy" is the forcing function for the Poisson equation.
                // In this case we set fxy to be a finite difference
                // Laplacian approximation to the (known) exact solution.
                //
                // We will use the second-order accurate FD Laplacian
                // approximation, which in 2D is
                //
                // u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
                //                u(i-1,j) + u(i+1,j) +
                //                -4*u(i,j))/h^2
                //
                // Since the value of the forcing function depends only
                // on the location of the quadrature point (q_point[qp])
                // we will compute it here, outside of the i-loop
                const Real fxy = 0;

                for (unsigned int i = 0; i != n_dofs; i++)
                {
                    Fe(i) += JxW[qp] * fxy * phi[i][qp];
                }
            }
        }

        // We have now reached the end of the RHS summation,
        // and the end of quadrature point loop, so
        // the interior element integration has
        // been completed.  However, we have not yet addressed
        // boundary conditions.  For this example we will only
        // consider simple Dirichlet boundary conditions.
        //
        // There are several ways Dirichlet boundary conditions
        // can be imposed.  A simple approach, which works for
        // interpolary bases like the standard Lagrange polynomials,
        // is to assign function values to the
        // degrees of freedom living on the domain boundary. This
        // works well for interpolary bases, but is more difficult
        // when non-interpolary (e.g Legendre or Hierarchic) bases
        // are used.
        //
        // Dirichlet boundary conditions can also be imposed with a
        // "penalty" method.  In this case essentially the L2 projection
        // of the boundary values are added to the matrix. The
        // projection is multiplied by some large factor so that, in
        // floating point arithmetic, the existing (smaller) entries
        // in the matrix and right-hand-side are effectively ignored.
        //
        // This amounts to adding a term of the form (in latex notation)
        //
        // \frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i
        //
        // where
        //
        // \frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1

        {

            // The following loop is over the sides of the element.
            // If the element has no neighbor on a side then that
            // side MUST live on a boundary of the domain.
            for (auto side : elem->side_index_range()){
                if (elem->neighbor_ptr(side) == nullptr)
                {
                    // get boundary info
                    //auto boundaryid = mesh.get_boundary_info().boundary_id(elem, side);
                    unsigned int n_boundary_ids=mesh.boundary_info->n_boundary_ids(elem,side);
                    std::vector<short int> boundary_ids_vec(n_boundary_ids);
                    mesh.boundary_info->boundary_ids(elem,side, boundary_ids_vec);
                    auto boundaryid=boundary_ids_vec[0];
                    //std::cout << boundaryid << "\n";

                    //check if boundaryid of element elem is in the Dirichlet BC ID list

                    double elem_dirichlet_bc=0;
                    bool found_side_set_in_list =0;
                    double radius;

                    //Check if side set is on the Idirichlet_side_set_list from input
                    for(int it=0; it<dirichlet_side_set_list.size(); it++)
                    {
                    	if(dirichlet_side_set_list[it]==boundaryid){
                    	elem_dirichlet_bc=dirichlet_bc_list[it];
                    	found_side_set_in_list =1;
                    	radius = dirichlet_id_radius[it];
                    	}
                    }


                    // The value of the shape functions at the quadrature
                    // points.
                    const std::vector<std::vector<Real>>& phi_face = fe_face->get_phi();

                    // The Jacobian * Quadrature Weight at the quadrature
                    // points on the face.
                    const std::vector<Real>& JxW_face = fe_face->get_JxW();

                    // The XYZ locations (in physical space) of the
                    // quadrature points on the face.  This is where
                    // we will interpolate the boundary value function.
                    const std::vector<Point>& qface_point = fe_face->get_xyz();

                    // Compute the shape function values on the element
                    // face.
                    fe_face->reinit(elem, side);

                    // Some shape functions will be 0 on the face, but for
                    // ease of indexing and generality of code we loop over
                    // them anyway
                    libmesh_assert_equal_to(n_dofs, phi_face.size());

                    // Loop over the face quadrature points for integration.
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        // The location on the boundary of the current
                        // face quadrature point.
                        const Real xf = qface_point[qp](0);
                        const Real yf = qface_point[qp](1);
                        bool is_on_neighborhood = 0;
                        bool is_on_set_fossa =0;
                        bool is_on_set_laa =0;


                        // Verify if the quadrature points are within the neighborhood of one of the id list points
                        for (int jj = 0; jj < dirichlet_id_points.size(); jj++)
                        {
                            libMesh::Point p_id(qface_point[qp] - dirichlet_id_points[jj]);

                            if ((dirichlet_id_radius[jj]>0) && (dirichlet_id_radius[jj]<47) )
                            {
								is_on_neighborhood = (p_id.norm() <= dirichlet_id_radius[jj]) ? 1 : 0;
							    // imposing BC on sphere
								if (is_on_neighborhood == 1){
										elem_dirichlet_bc=dirichlet_bc_list[jj];
										radius = dirichlet_id_radius[jj];
										break;
								}
							}

                        }


                        // The penalty value.  \frac{1}{\epsilon}
                        // in the discussion above.
                        const Real penalty = 1.e10;

                        // The boundary value.
                       //const Real value = exact_solution(xf, yf);

                        bool cond1 = (found_side_set_in_list && (radius<0)) ? 1 : 0;
                        bool cond2 = (found_side_set_in_list && (radius>0) && (radius<47) && is_on_neighborhood)? 1 : 0;


                        // Matrix contribution of the L2 projection.
                        if(cond1 || cond2) //
                        {
                            for (unsigned int i = 0; i != n_dofs; i++)
                            {
                                for (unsigned int j = 0; j != n_dofs; j++)
                                    Ke(i, j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];
                            }
                        }
                        // Right-hand-side contribution of the L2
                        // projection.
                        for (unsigned int i = 0; i != n_dofs; i++)
                        {
                             if(cond1 || cond2)
                             {
                                Fe(i) += JxW_face[qp] * penalty * elem_dirichlet_bc * phi_face[i][qp]; //dirichlet
                             }
                            else
                                Fe(i) += 0;         //neumann
                        }
                    }
                }
            }
        }

        // We have now finished the quadrature point loop,
        // and have therefore applied all the boundary conditions.

        // If this assembly program were to be used on an adaptive mesh,
        // we would have to apply any hanging node constraint equations
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);

        // The element matrix and right-hand-side are now built
        // for this element.  Add them to the global matrix and
        // right-hand-side vector.  The  SparseMatrix::add_matrix()
        // and  NumericVector::add_vector() members do this for us.
        system.matrix->add_matrix(Ke, dof_indices);
        system.rhs->add_vector(Fe, dof_indices);

    }

	// implement one point bc condition

    double radius;
    libMesh::Point p;
    double x;
    double y;
    double z;
    double bc_value;
    for (int jj = 0; jj < dirichlet_id_points.size(); jj++) {
    		   if(dirichlet_id_radius[jj]==49){
    			   p = dirichlet_id_points[jj];
    			   radius = anatomic_parameters.radius_LAA_tip; // small enough to enclose only one node
    			   x = p(0);
    			   y = p(1);
    			   z = p(2);
    			   bc_value = dirichlet_bc_list[jj];
    			   std::cout << "x= " <<x << " y= " <<y<< " z= " <<z<< " r= " <<radius << " bc_value= " <<bc_value   << std::endl;
               }
    }
    // loop over all the nodes
	for(auto & node : mesh.local_node_ptr_range() ){
		libMesh::Point point = * node;
        libMesh::Point laa_point(x,y,z);
        laa_point -= point;
        if(laa_point.norm()< radius){
			unsigned int dn = node->dof_number(system.number(), 0, 0);
			double penalty = 1e8;
			system.rhs->add(dn, bc_value*penalty);
			system.matrix->add(dn, dn, penalty);
	}

}

}


// MAIN

int main(int argc, char ** argv)
{
    // Read input file
    GetPot data = BeatIt::readInputFile(argc, argv);

    // Initialize libMesh
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // Create empty mesh
    ParallelMesh mesh(init.comm());
    mesh.allow_renumbering(false);

    std::cout<< "Reading mesh" << std::endl;
    // Read mesh filename from input file
    std::string mesh_name = data("mesh", "NONE");

    // Modifying code to solve all problems at once
    int N= data("number_of_problems", 0);

    //Create an object of class Poisson
    Poisson poisson_solver;

    // Read name of the input files corresponding to each problem
    std::string input_list_aux = data("inputs", "");
    BeatIt::readList(input_list_aux, poisson_solver.input_list);

    // This example requires a linear solver package.
    libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
        "--enable-petsc, --enable-trilinos, or --enable-eigen");

    // Brief message to the user regarding the program name
    // and command line arguments.
    libMesh::out << "Running " << argv[0];

    // This fiber example requires a linear solver package.
    libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
        "--enable-petsc, --enable-trilinos, or --enable-eigen");


    // Brief message to the user regarding the program name
        // and command line arguments.
        libMesh::out << "Running " << argv[0];

        for (int i = 1; i < argc; i++)
            libMesh::out << " " << argv[i];

        libMesh::out << std::endl << std::endl;

        // Skip this 2D example if libMesh was compiled as 1D-only.
        libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

        // read from input
        // Get the ids where Dirichlet BC are imposed
        std::string ids_dirichlet_test = data("dirichletbc", "");  //> did not work
        std::cout << ids_dirichlet_test << " ";
        int size_dirichlet_bcs = data.vector_variable_size("dirichletbc");
        std::vector<int> ids_dirichlet_bcs;
        for (int i = 0; i < size_dirichlet_bcs; i++) {
            ids_dirichlet_bcs.push_back(data("dirichletbc", i, -1));
        }

        // BeatIt::readList(ids_dirichlet_test, ids_dirichlet_bcs);               // IMPORTANT LINE HERE
         //std::cout <<"ids dirichlet:"<< ids_dirichlet_bcs << "\n";
        for (int i = 0; i < size_dirichlet_bcs; i++) {
            std::cout << ids_dirichlet_bcs[i] << " ";
        }
        // Read problem parameters
        std::string radius_laa_tip_aux = data("radius_LAA_tip"    , "0.5");
        //Fiber type is Milano's paper?
        //https://doi.org/10.1016/j.cma.2020.113468
        std::string milano_aux = data("milano", " 0");

        AnatomicalParameters anatomic_parameters;

        anatomic_parameters.reinit(milano_aux, radius_laa_tip_aux);


    //--------------------
    //    Used from example_laryssa and moved here

    // If we want to import the fiber from the meshfile we need to create
    // I/O object we will use to import the data
    // Before reading the fibers we need to create the correct
    // FE spaces to store them
    // This will be done after initializing the Electrophysiology solver
    libMesh::ExodusII_IO importer(mesh);

    // Create mesh:
    // If we passed a specific filename for the mesh let'd read that
    if ("NONE" != mesh_name)
    {
        int n_refinements = data("refs", 0);
        std::cout << "n_refs: " << n_refinements << std::endl;
        BeatIt::serial_mesh_partition(init.comm(), mesh_name, &mesh, n_refinements);
        std::cout << "after serial mesh partition" << std::endl;
    }
    // If no mesh file has been specified
    // we run on a cube
    else
    {
        // number of elements in the x,y and z direction
        int nelx = data("nelx", 10);
        int nely = data("nely", 10);
        // if nelz = 0 we run in 2D
        int nelz = data("nelz", 10);
        // the cube dimensions are defined as [minx, maxx] x [miny, maxy] x [minz, maxz]
        double maxx = data("maxx", 1.0);
        double maxy = data("maxy", 1.0);
        double maxz = data("maxz", 1.0);

        double minx = data("minx", 0.0);
        double miny = data("miny", 0.0);
        double minz = data("minz", 0.0);

        // Create a tetrahedral mesh
        auto elType = TET4;
        // If we are in 2D, create a triangle mesh
        if (nelz == 0)
            elType = TRI3;

        std::cout << "Creating the cube [" << minx << ", " << maxx << "] x ["
                                           << miny << ", " << maxy << "] x ["
                                           << minx << ", " << maxx << "] " << std::endl;
        std::cout << "Using " << nelx << " x " << nely << " x " << nelz << " elements." << std::endl;
        if(TET4 == elType) std::cout << "Element type TET4" << std::endl;
        else if(TRI3 == elType) std::cout << "Element type TRI3" << std::endl;
        else std::cout << "NO ELEMENT TYPE!!!" << std::endl;

        // Create mesh
        MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, minx, maxx, miny, maxy, minz, maxz, elType);
        // Usually we run using cm as dimensions
        // use this part to scale the mesh if it exported in mm (or um)
    }

    double scale = data("scale", 1.0);
    MeshTools::Modification::scale(mesh, scale, scale, scale);

    // output the details about the mesh
    mesh.print_info();

    //-----------------------------------------

    // Create an equation systems object.
    EquationSystems equation_systems(mesh);

    // Define fiber systems:
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& fiber_system = equation_systems.add_system<FiberSystem>("fibers");
    fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.init();
    auto& f_v = fiber_system.solution;
    FiberSystem& sheet_system = equation_systems.add_system<FiberSystem>("sheets");
    sheet_system.add_variable("sheetsx", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheet_system.add_variable("sheetsy", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheet_system.add_variable("sheetsz", libMesh::CONSTANT, libMesh::MONOMIAL);
    sheet_system.init();
    auto& s_v = sheet_system.solution;
    FiberSystem& xfiber_system = equation_systems.add_system<FiberSystem>("xfibers");
    xfiber_system.add_variable("xfibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.add_variable("xfibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.add_variable("xfibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    xfiber_system.init();
    auto& n_v = xfiber_system.solution;

    // Declare the Poisson system and its variables.
    // The Poisson system is another example of a steady system.
    LinearImplicitSystem& system = equation_systems.add_system<LinearImplicitSystem>("Poisson");

    // Adds the variable "u" to "Poisson".  "u"
    // will be approximated using second-order approximation.
    if ("NONE" != mesh_name)
    { equation_systems.get_system("Poisson").add_variable("u", FIRST); //TET4 ELEMENTS
    }
    else {
    	throw std::runtime_error("NO MESH AS INPUT!");
    }

    // Initialize the data structures for the equation system.
    //equation_systems.init();
    system.init();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    // Reading the input for each of the problems
    std::cout << "Reading the input for each of the" << N << "  problems" <<std::endl;
    for (int i = 0; i < N; i++) {

        std::string ui = "u" + std::to_string(i);    //how cool! to use i in the name here!!!
        ExplicitSystem& si = equation_systems.add_system<ExplicitSystem>(ui);
        si.add_variable(ui, FIRST);
        si.init();
        
        std::cout << "read input i=" << i << std::endl;
 	   //READ INPUT I > BE CAREFUL BECAUSE HERE I WILL BE OVERWRITTING WHAT I HAD BEFORE
 	   poisson_solver.read_input_i(i, equation_systems);

        std::cout << " reinit " << std::endl;
       // reinit Poisson member variables using input file values
 	  poisson_solver.reinit(poisson_solver.dirichlet_x_list, poisson_solver.dirichlet_y_list,  poisson_solver.dirichlet_z_list,
 			  poisson_solver.dirichlet_r_list,  poisson_solver.dirichlet_side_set_list_aux, poisson_solver.dirichlet_bc_list_aux);

      std::cout << "assemble_poisson" << std::endl;
 	  poisson_solver.assemble_poisson(equation_systems, "Poisson", i, anatomic_parameters);

 /*       //  Libmesh zeros out matrix and RHS when i call solve
        // preventing that to happen:
        equation_systems.get_system<ImplicitSystem>("Poisson").zero_out_matrix_and_rhs=false;
        std::cout << "solve poisson system i= "<< i << std::endl;
        equation_systems.get_system("Poisson").solve();


        // the processor number
		int rank = si.comm().rank();
		// size of solution vector at each mesh node - global vector (including all the processors)
		int size_vector = si.solution->size();
		// size of solution vector at node in this processor - local vector
		int local_size = si.solution->local_size();
		int size_vector_poisson = equation_systems.get_system("Poisson").solution->size();
		int local_size_poisson = equation_systems.get_system("Poisson").solution->local_size();
		std::cout << "rank=" << rank << " , local_size=" <<local_size <<" , local_size_poisson=" <<local_size_poisson  << std::endl;
		std::cout << "rank=" << rank << " , size_vector=" << size_vector << " , size_vector_poisson=" << size_vector_poisson << std::endl;
		//si.comm().barrier(); //all the processors have to be synchronized by this point before it proceeds
		si.solution->close();

        // Copy the content of "u" (solution to Poisson problem) to "ui"
        *si.solution = *equation_systems.get_system("Poisson").solution; // "unique vector"
        si.update(); //distributes/copies the solution of the global on the shared interface > local processor vectors - "repeated vector"
*/
    }
    std::cout << "Read output folder name" << std::endl;
    // Read output folder name from input file
    std::string output_folder_name = data("output_name", "Output_default");
    std::string output_number = data("output_number", "");
    std::string test = output_folder_name + output_number;
    BeatIt::createOutputFolder(test);

    // Output file name
    std::string output_file_name = data("output_file", "/out");

    // Solve the system "Poisson".  Note that calling this
    // member will assemble the linear system and invoke
    // the default numerical solver.  With PETSc the solver can be
    // controlled from the command line.  For example,
    // you can invoke conjugate gradient with:
    //
    // ./executable -ksp_type cg
    //
    // You can also get a nice X-window that monitors the solver
    // convergence with:
    //
    // ./introduction-ex3 -ksp_xmonitor
    //
    // if you linked against the appropriate X libraries when you
    // built PETSc.
/*
    const DofMap& dofmap = system.get_dof_map(); // replicating assembly function

    std::cout << " Looping over elements to set up subdomain IDs"<< std::endl;

    for (const auto& elem: mesh.active_local_element_ptr_range()) {
       // elem->subdomain_id()=1;
        int blockid = elem->subdomain_id();
        libMesh::Point centroid = elem->vertex_average();
        std::vector<double> u(N);
        std::vector<libMesh::Gradient> du(N);

        // loop over the number of laplace problems we are solving
        for (int i = 0; i < N; i++)
        {
            std::string ui = "u" + std::to_string(i);
            ExplicitSystem& si = equation_systems.get_system<ExplicitSystem>(ui);
            u[i] = si.point_value(si.variable_number(ui), centroid, *elem);
            du[i] = si.point_gradient(si.variable_number(ui), centroid, *elem);
        }

        //Define anatomical regions
        blockid = anatomic_parameters.define_regions(u);

        elem->Elem::subdomain_id()=blockid; // this line does not work


        // Setup fibers
        anatomic_parameters.define_fibers(u, du, blockid, elem);

        elem->Elem::subdomain_id()=1; //reset block ids so that no problem with EP

        std::vector<dof_id_type> fibers_dof_indices; //putting this back to the system so that we can export it

        fiber_system.get_dof_map().dof_indices(elem, fibers_dof_indices);
        for (int idim = 0; idim < 3; ++idim)
        {
            fiber_system.solution ->set(fibers_dof_indices[idim], anatomic_parameters.f0(idim));
            sheet_system.solution ->set(fibers_dof_indices[idim], anatomic_parameters.s0(idim));
            xfiber_system.solution->set(fibers_dof_indices[idim], anatomic_parameters.n0(idim));
        }

    }

    // After solving the system write the solution
    // to a VTK-formatted plot file.
    std::cout << "saving fibers" << std::endl;
    save_fibers(test, output_file_name, output_number, equation_systems, mesh);

*/

/*    // Information about the time, such as
    // current iteration, timestep, final time etc
    // can be conveniently read from the input file
    // and stored in a TimeData abject
    // In the input file specify
    //
    //    [time]
    //        # Timestep
    //        dt = 0.125      # Default: 1.0
    //        # Simulation initial time
    //        init_time = 0.0 # Default: 0.0
    //        # Simulation end time
    //        final_time = 10  # Default: 1.0
    //        # Maximum number of timesteps
    //        max_iter = 200000000   # Default: 99999999
    //        # Export the solution every save_iter iterations
    //        save_iter = 8   # Default: 1
    //    [../]
    //
    // Create the TimeData object
    BeatIt::TimeData datatime;
    // Set it up using the input file
    datatime.setup(data, "");
    // Output on screen the stored variables
    datatime.print();

    // Create libMesh Equations systems
    // This will hold the mesh and create the corresponding
    // finite element spaces
    //libMesh::EquationSystems es(mesh);


    // Create the Electrophyiosiology model based on the input file
    //
    //   model = monowave
    //   # The input parameters of the model are written
    //   # in a section with the same name as 'model'
    //   [monowave]
    //       # PARAMETERS HERE
    //   [../]
    //
    // Read the model from the input file
    std::string model = data("model", "monowave");
    std::cout << "Create Electrophysiology model ..." << std::endl;
    // Create the model
    BeatIt::ElectroSolver* solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, equation_systems);
    // Setup the EP model using the input file
    std::cout << "Calling setup..." << std::endl;

    solver->setup(data, model);

    // Initialize systems
    std::cout << "Eqution_systems info ..." << std::endl;
    equation_systems.print_info();

    // Set up initial conditions at time
    std::cout << "Calling init ..." << std::endl;
    solver->init(datatime.M_startTime);
    // Now read the fibers if wanted

    bool read_fibers = data("read_fibers", false);
    if(read_fibers)
    {
    	std::cout << "~~~~ ~~~~~~~~~~Reading fibers from input file~~~~ ~~~~~~~~~~" << std::endl;
        // First show the elemental variables that can be imported
        auto elemental_variables = importer.get_elem_var_names();
        for (auto && var : elemental_variables) std::cout << var << std::endl;
        solver->read_fibers(importer,1); //solver->fiber_system?
    }

    if (solver == NULL){
    	std::cout << " !!!!!!!!!!!!         SOLVER IS A NULL PTR !!!!!!!!" << std::endl;
    }

    // Export simulation parameters
    // This will also export the fiber field
    //    solver->save_parameters();
    // Assemble matrices
    std::cout << "Assembling matrices" << std::endl;

    solver->assemble_matrices(datatime.M_dt);                           // PROBLEM HERE WHEN RUN IN PARALLEL
    // output file counter

    int save_iter = 0;
    // Export initial condition at time
//    solver->save_exo_timestep(save_iter, datatime.M_time);
    solver->save_potential(save_iter, datatime.M_startTime);
//    solver->save_parameters();

    // Parameters to save the activation times
    // A node is activated if the transmembrane potential > threshold
    double threshold = data("threshold", -10.0);
    // We export the activation times every at_save_iter iterations
    int at_save_iter = data("at_save_iter", 25);

    // Old parameters that define the method, not to be changed
    // TODO: clean this part
    std::string system_mass = data(model + "/diffusion_mass", "mass");
    std::string iion_mass = data(model + "/reaction_mass", "lumped_mass");
    bool useMidpointMethod = false;
    int step0 = 0;
    int step1 = 1;

    //solver->save_parameters();


    // Start loop in time
    std::cout << "Time loop starts:" << std::endl;
    // Control the time loop using the TimeData object
    for (; datatime.M_iter < datatime.M_maxIter && datatime.M_time < datatime.M_endTime;)
    {
        // We are doing a new iteration
        // let's update first the information in the TimeData object
        datatime.advance();
        // Ouput to screen the current time and the current iteration
        std::cout << "Time:" << datatime.M_time << ", Iter: " << datatime.M_iter << std::endl;
        // Advance the solution in time: u_n <- u_n+1
        solver->advance();
        // Solve ionic model and evaluate ionic currents
        solver->solve_reaction_step(datatime.M_dt, datatime.M_time, step0, useMidpointMethod, iion_mass);
        // Solve monodomain model
        solver->solve_diffusion_step(datatime.M_dt, datatime.M_time, useMidpointMethod, iion_mass);
        // Update the activation times
        solver->update_activation_time(datatime.M_time, threshold);
        //Export the solution if at the right timestep
        if (0 == datatime.M_iter % datatime.M_saveIter)
        {
            // update output file counter
            save_iter++;
            // export current solution
            solver->save_potential(save_iter, datatime.M_time);
//            solver->save_exo_timestep(save_iter, datatime.M_time);
        }
        // export the activation times if at the corresponding timestep
        if (0 == datatime.M_iter % (at_save_iter * datatime.M_saveIter))
        {
            // export activation times
            solver->save_activation_times(save_iter);
        }

    }

    // // export activation times again
    solver->save_activation_times(save_iter);

    // delete solver before ending the simulation
    // avoiding memory leaks
    delete solver;

*/
    // The end
    std::cout << "Good luck with your simulation :P" << std::endl;
    return 0;
}
