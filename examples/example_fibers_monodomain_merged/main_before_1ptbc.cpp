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
	    std::vector<std::string> output_variables(13); // creates a vector that stores the names of vectors
	    output_variables[0] = "fibersx";
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
	    output_variables[12] = "u3";
	    nemesis_exporter.set_output_variables(output_variables);
	    nemesis_exporter.write_equation_systems(e_test + e_output_file_name + e_output_number + ".e", equation_systems); // saves the variables at the nodes
	    nemesis_exporter.write_element_data(equation_systems); // this saves the variables at the centroid of the elements
}

// The AnatomicalParameters constains all parameters that have to do with
// region selection and building the fiber
class AnatomicalParameters{
	public:

	// this info is read from the main input file (input.main)
	// i_region is used to define the anatomical regions
	// i_epic_region is used to define the fibers in the anatomical region on the epicardium
	// i_endo_region is used to define the fibers in the anatomical region on the endocardium

	double milano; // could be boolean
	int i_PV1, i_epic_PV1, i_endo_PV1;
	int i_PV2, i_epic_PV2, i_endo_PV2;
	int i_PV3, i_epic_PV3, i_endo_PV3;
	int i_PV4, i_epic_PV4, i_endo_PV4;
	int i_antra1, i_epic_antra1_rspv, i_endo_antra1_rspv;
	int           i_epic_antra1_ripv, i_endo_antra1_ripv;
	int i_antra2, i_epic_antra2, i_endo_antra2;
	int i_carina1, i_epic_carina1, i_endo_carina1;
	int i_carina2, i_epic_carina2, i_endo_carina2;
	int i_floor, i_epic_floor, i_endo_floor;
	int i_LAA, i_epic_LAA, i_endo_LAA;
	int i_fossa, i_epic_fossa, i_endo_fossa;
	int i_strip, i_epic_strip, i_endo_strip;
	int i_ripv, i_epic_ripv, i_endo_ripv;
	int i_septum, i_epic_septum, i_endo_septum;
	int i_anterior, i_epic_anterior, i_endo_anterior;
	int i_posterior, i_epic_posterior, i_endo_posterior;
	int i_lateral, i_epic_lateral, i_endo_lateral;
	int i_epic_posterior_bottom, i_endo_posterior_bottom;
	int i_epic_anterior_bottom, i_endo_anterior_bottom;
	int i_epic_septum_bottom_rspv, i_endo_septum_bottom_rspv;
	int i_epic_septum_bottom_ripv, i_endo_septum_bottom_ripv;
	double thresholdA_PV1, thresholdB_PV1;
	double thresholdA_PV2, thresholdB_PV2;
	double thresholdA_PV3, thresholdB_PV3;
	double thresholdA_PV4, thresholdB_PV4;
	double thresholdA_antra1_rspv, thresholdB_antra1_rspv;
	double thresholdA_antra1_ripv, thresholdB_antra1_ripv;
	double thresholdA_antra2, thresholdB_antra2;
	double thresholdA_carina1, thresholdB_carina1;
	double thresholdA_carina2, thresholdB_carina2;
	double thresholdA_floor, thresholdB_floor;
	double thresholdA_LAA, thresholdB_LAA;
	double thresholdA_septum, thresholdB_septum;
	double thresholdA_anterior, thresholdB_anterior;
	double thresholdA_posterior, thresholdB_posterior;
	double thresholdA_lateral, thresholdB_lateral;
	double thresholdA_strip, thresholdB_strip;
	double thresholdA_ripv, thresholdB_ripv;
	double thresholdA_fossa, thresholdB_fossa;
	libMesh::Gradient f0, s0, n0;

	//Reinitialize variables
	void reinit(std::string, std::string, std::string, std::string,
				std::string, std::string, std::string, std::string,std::string,
				std::string, std::string, std::string, std::string,std::string,
				std::string, std::string, std::string, std::string,
				std::string, std::string, std::string, std::string, std::string);

	// Define regions
	int define_regions (std::vector<double> &);
	// function to define fibers  | input: u, output: fiber field??
	void define_fibers(std::vector<double> & , std::vector<libMesh::Gradient> & , int, const auto&);
};

void AnatomicalParameters::reinit(std::string Ithreshold_pv1, std::string threshold_pv2, std::string threshold_pv3, std::string threshold_pv4,
									std::string threshold_floor, std::string threshold_laa, std::string threshold_antra1_rspv, std::string threshold_antra1_ripv, std::string threshold_antra2,
									std::string threshold_septum, std::string threshold_anterior, std::string threshold_posterior,
									std::string threshold_carina1, std::string threshold_carina2, std::string threshold_lateral, std::string threshold_strip, std::string threshold_ripv,
									std::string threshold_septum_bottom_rspv, std::string threshold_septum_bottom_ripv, std::string threshold_anterior_bottom,
									std::string threshold_posterior_bottom, std::string threshold_fossa, std::string milano_aux)
	{
	std::vector<double> tpv1, tpv2, tpv3, tpv4, tfloor, tlaa, tantra1_rspv, tantra1_ripv, tantra2, tseptum, tanterior, tposterior, tcarina1, tcarina2, tlateral, tstrip;
	std::vector<double> tanterior_bottom, tposterior_bottom, tseptum_bottom_rspv, tseptum_bottom_ripv, tfossa, tripv;
	std::vector<double> milano_vec;

	BeatIt::readList(milano_aux, milano_vec);
	milano = milano_vec[0];

	BeatIt::readList(Ithreshold_pv1, tpv1);
    i_PV1          = tpv1[0];
    thresholdA_PV1 = tpv1[1];
    thresholdB_PV1 = tpv1[2];
    i_epic_PV1     = tpv1[3];
    i_endo_PV1     = tpv1[4];

    BeatIt::readList(threshold_pv2, tpv2);
    i_PV2          = tpv2[0];
    thresholdA_PV2 = tpv2[1];
    thresholdB_PV2 = tpv2[2];
    i_epic_PV2     = tpv2[3];
    i_endo_PV2     = tpv2[4];

    BeatIt::readList(threshold_pv3, tpv3);
    i_PV3          = tpv3[0];
    thresholdA_PV3 = tpv3[1];
    thresholdB_PV3 = tpv3[2];
    i_epic_PV3     = tpv3[3];
    i_endo_PV3     = tpv3[4];

    BeatIt::readList(threshold_pv4, tpv4);
    i_PV4          = tpv4[0];
    thresholdA_PV4 = tpv4[1];
    thresholdB_PV4 = tpv4[2];
    i_epic_PV4     = tpv4[3];
    i_endo_PV4     = tpv4[4];

    BeatIt::readList(threshold_floor,tfloor );
    i_floor = tfloor[0];
    thresholdA_floor = tfloor[1];
    thresholdB_floor = tfloor[2];
    i_epic_floor     = tfloor[3];
    i_endo_floor     = tfloor[4];

    BeatIt::readList(threshold_laa,tlaa );
    i_LAA          = tlaa[0];
    thresholdA_LAA = tlaa[1];
    thresholdB_LAA = tlaa[2];
    i_epic_LAA     = tlaa[3];
    i_endo_LAA     = tlaa[4];

    BeatIt::readList(threshold_fossa,tfossa );
    i_fossa          = tfossa[0];
    thresholdA_fossa = tfossa[1];
    thresholdB_fossa = tfossa[2];
    i_epic_fossa     = tfossa[3];
    i_endo_fossa     = tfossa[4];

    BeatIt::readList(threshold_strip,tstrip );
    i_strip          = tstrip[0];
    thresholdA_strip = tstrip[1];
    thresholdB_strip = tstrip[2];
    i_epic_strip     = tstrip[3];
    i_endo_strip     = tstrip[4];

    BeatIt::readList(threshold_ripv,tripv );
    i_ripv          = tripv[0];
    thresholdA_ripv = tripv[1];
    thresholdB_ripv = tripv[2];
    i_epic_ripv     = tripv[3];
    i_endo_ripv     = tripv[4];

    BeatIt::readList(threshold_antra1_rspv, tantra1_rspv);
    i_antra1          = tantra1_rspv[0];
    thresholdA_antra1_rspv = tantra1_rspv[1];
    thresholdB_antra1_rspv = tantra1_rspv[2];
    i_epic_antra1_rspv = tantra1_rspv[3];
    i_endo_antra1_rspv = tantra1_rspv[4];

    BeatIt::readList(threshold_antra1_ripv, tantra1_ripv);
    thresholdA_antra1_ripv = tantra1_ripv[1];
    thresholdB_antra1_ripv = tantra1_ripv[2];
    i_epic_antra1_ripv = tantra1_ripv[3];
    i_endo_antra1_ripv = tantra1_ripv[4];
    std::cout << "tantra1_ripv[3]=" << tantra1_ripv[3] << ", tantra1_ripv[4]" << tantra1_ripv[4] <<std::endl;

    BeatIt::readList(threshold_antra2, tantra2);
    i_antra2          = tantra2[0];
    thresholdA_antra2 = tantra2[1];
    thresholdB_antra2 = tantra2[2];
    i_epic_antra2     = tantra2[3];
    i_endo_antra2     = tantra2[4];

    BeatIt::readList(threshold_septum, tseptum);
    i_septum          = tseptum[0];
    thresholdA_septum = tseptum[1];
    thresholdB_septum = tseptum[2];
    i_epic_septum     = tseptum[3];
    i_endo_septum     = tseptum[4];

    BeatIt::readList(threshold_anterior, tanterior);
    i_anterior          = tanterior[0];
    thresholdA_anterior = tanterior[1];
    thresholdB_anterior = tanterior[2];
    i_epic_anterior     = tanterior[3];
    i_endo_anterior     = tanterior[4];

    BeatIt::readList(threshold_posterior, tposterior);
    i_posterior          = tposterior[0];
    thresholdA_posterior = tposterior[1];
    thresholdB_posterior = tposterior[2];
    i_epic_posterior     = tposterior[3];
    i_endo_posterior     = tposterior[4];

    BeatIt::readList(threshold_carina1, tcarina1);
    i_carina1          = tcarina1[0];
    thresholdA_carina1 = tcarina1[1];
    thresholdB_carina1 = tcarina1[2];
    i_epic_carina1     = tcarina1[3];
    i_endo_carina1     = tcarina1[4];

    BeatIt::readList(threshold_carina2, tcarina2);
    i_carina2          = tcarina2[0];
    thresholdA_carina2 = tcarina2[1];
    thresholdB_carina2 = tcarina2[2];
    i_epic_carina2     = tcarina2[3];
    i_endo_carina2     = tcarina2[4];

    BeatIt::readList(threshold_lateral, tlateral);
    i_lateral          = tlateral[0];
    thresholdA_lateral = tlateral[1];
    thresholdB_lateral = tlateral[2];
    i_epic_lateral     = tlateral[3];
    i_endo_lateral     = tlateral[4];

    BeatIt::readList(threshold_septum_bottom_rspv, tseptum_bottom_rspv);
    i_epic_septum_bottom_rspv     = tseptum_bottom_rspv[3];
    i_endo_septum_bottom_rspv     = tseptum_bottom_rspv[4];

    BeatIt::readList(threshold_septum_bottom_ripv, tseptum_bottom_ripv);
    i_epic_septum_bottom_ripv     = tseptum_bottom_ripv[3];
    i_endo_septum_bottom_ripv     = tseptum_bottom_ripv[4];

    BeatIt::readList(threshold_anterior_bottom, tanterior_bottom);
    i_epic_anterior_bottom   = tanterior_bottom[3];
    i_endo_anterior_bottom   = tanterior_bottom[4];

    BeatIt::readList(threshold_posterior_bottom, tposterior_bottom);
    i_epic_posterior_bottom  = tposterior_bottom[3];
    i_endo_posterior_bottom  = tposterior_bottom[4];

}

int AnatomicalParameters::define_regions (std::vector<double> & u){
	int blockid;
	//MILANO
	if(milano==1){
		if ((u[3]>=0.7)){     //MV
			blockid = 1;
		}
		else if ((u[2]>=0.8) || (u[2]<=0.35)) { // PVs
			blockid = 2;
		}
		else
			blockid = 3;
	}
	//OUR MODEL
	else{
		/*if ((u[i_PV1]>thresholdA_PV1)*(u[i_PV1]<thresholdB_PV1)){
			blockid = 12;
		}
		else if ((u[i_PV2]>thresholdA_PV2)*(u[i_PV2]<thresholdB_PV2)) { 						//pv
			blockid = 1;
		}
		else if ((u[i_PV3]>thresholdA_PV3)*(u[i_PV3]<thresholdB_PV3)) { 						//pv
			blockid = 2;
		}
		else if ((u[i_PV4]>thresholdA_PV4)*(u[i_PV4]<thresholdB_PV4)) { 						//pv
			blockid = 3;
		}*/
		if ((u[i_floor]>thresholdA_floor)*(u[i_floor]<thresholdB_floor)) { 				//floor
			blockid = 4;
		}
		else{
			if ((u[i_LAA]>thresholdA_LAA)*(u[i_LAA]<thresholdB_LAA)) {							       //laa
				blockid = 5;
			}
			else if ((u[i_fossa]>thresholdA_fossa)*(u[i_fossa]<thresholdB_fossa)) {							//fossa
				blockid = 18;
			}
			else if ((u[i_antra1]>thresholdA_antra1_rspv)*(u[i_antra1]<thresholdB_antra1_rspv)&&
					((u[i_ripv]>thresholdB_ripv)||(u[i_ripv]<thresholdA_ripv))) { 	                      //antra1 - rspv
				blockid = 6;
			}
			else if ((u[i_antra1]>thresholdA_antra1_ripv)*(u[i_antra1]<thresholdB_antra1_ripv)&&
					(u[i_ripv]<thresholdB_ripv)*(u[i_ripv]>thresholdA_ripv)) { 			                  //antra1 - ripv
				blockid = 16;
			}
			else if ((u[i_lateral]>thresholdA_lateral)*(u[i_lateral]<thresholdB_lateral)&& u[i_strip]>thresholdB_strip) { //lateral
				blockid = 7;
			}
			else if ((u[i_antra2]>thresholdA_antra2)*(u[i_antra2]<thresholdB_antra2)) { 			//antra2
				blockid = 8;
			}
			else if (//(u[i_strip]>thresholdA_strip)*(u[i_strip]<thresholdB_strip) &&
					(u[i_septum]>thresholdA_septum)*(u[i_septum]<thresholdB_septum) &&
					((u[i_ripv]>thresholdB_ripv)||(u[i_ripv]<thresholdA_ripv))) {                    //septum bottom rspv
				blockid = 13;
			}
			else if (//(u[i_strip]>thresholdA_strip)*(u[i_strip]<thresholdB_strip) &&
					(u[i_septum]>thresholdA_septum)*(u[i_septum]<thresholdB_septum) &&
					(u[i_ripv]<thresholdB_ripv)*(u[i_ripv]>thresholdA_ripv)) {                     //septum bottom ripv
				blockid = 17;
			}
			/*else if ((u[i_septum]>thresholdA_septum)*(u[i_septum]<thresholdB_septum)) { 			//septum
				blockid = 9;
			}*/
			else if ((u[i_strip]>thresholdB_strip) &&
					(u[i_anterior]>thresholdA_anterior)*(u[i_anterior]<thresholdB_anterior)&&
					u[i_ripv]>thresholdB_ripv) {	                                              //anterior_bottom
						blockid = 14;
					}
			else if ((u[i_anterior]>thresholdA_anterior)*(u[i_anterior]<thresholdB_anterior) &&
					u[i_ripv]>thresholdB_ripv) {                                                    //anterior
				blockid = 10;
			}
			else if ((u[i_strip]>thresholdA_strip)*(u[i_strip]<thresholdB_strip)) //&&
					//(u[i_posterior]>thresholdA_posterior)*(u[i_posterior]<thresholdB_posterior)) { //posterior bottom
			{
					blockid = 15;
			}
			else {               																	  //posterior
				blockid = 11;
			}
		}
	}
    return blockid;
}

void AnatomicalParameters::define_fibers(std::vector<double> & u, std::vector<libMesh::Gradient> & du, int block_id, const auto& elem){

    // Solve the transmural problem as the first one, such that
    // the next line is always true
	s0 = du[0].unit();

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
	// OUR MODEL-----------------------
	else{
		// Fibers created so we have epi and endo layers in some regions
		switch (block_id)
		{
		//PV left front
		case 12:
		{
				n0 = du[i_epic_PV1].unit();
				f0 = s0.cross(n0);
				if(i_epic_PV1==0 && i_endo_PV1==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				else if(i_endo_PV1==90){
					f0 = du[i_epic_PV1].unit();
					n0 = s0.cross(f0);
				}
				break;
		}
		//PV left back
		case 1:
		{
				n0 = du[i_epic_PV2].unit();
				f0 = s0.cross(n0);
				if(i_epic_PV2==0 && i_endo_PV2==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				else if(i_endo_PV2==90){
					f0 = du[i_epic_PV2].unit();
					n0 = s0.cross(f0);
				}
				break;
		}
		//PV right front
		case 2:
		{
				n0 = du[i_epic_PV3].unit();
				f0 = s0.cross(n0);
				if(i_epic_PV3==0 && i_endo_PV3==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				else if(i_endo_PV3==90){
					f0 = du[i_epic_PV3].unit();
					n0 = s0.cross(f0);
				}
				break;
		}
		//PV right front
		case 3:
		{
				n0 = du[i_epic_PV4].unit();
				f0 = s0.cross(n0);
				if(i_epic_PV4==0 && i_endo_PV4==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				else if(i_endo_PV4==90){
					f0 = du[i_epic_PV4].unit();
					n0 = s0.cross(f0);
				}
				break;
		}
		// Floor
		case 4:
		{
			//if(u[0]>0.5){
				n0 = du[i_epic_floor].unit();
				f0 = s0.cross(n0);
				if(i_epic_floor==0 && i_endo_floor==0){
					f0(0) = 1.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
			   else if(i_endo_floor==90){
				   f0 = du[i_epic_floor].unit();
				   n0 = s0.cross(f0);
				}
			break;
		}
		// LAA
		case 5:
		{
			//if(u[0]>0.5){
			if 	(i_endo_LAA!=90 && i_epic_LAA!=90 ){
				n0 = du[i_epic_LAA].unit();
				f0 = s0.cross(n0);
				break;
					}                          // endocardium 90 degrees wrt to epicardium
//        	else if(u[0]<=0.5 && i_endo_LAA!=90){
			else if(i_epic_LAA==0 && i_endo_LAA==0){
				f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
				s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
				n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				break;
			}
			else{
				f0 = du[i_epic_LAA].unit();
				n0 = s0.cross(f0);
				break;
			}
			/*else{
				f0 = du[i_endo_LAA].unit();
				n0 = s0.cross(f0);
				break;
			}*/
		}
		// fossa
		case 18:
		{
			if 	(i_endo_fossa!=90 && i_epic_fossa!=90 ){
				n0 = du[i_epic_fossa].unit();
				f0 = s0.cross(n0);
				break;
					}
			else if(i_epic_fossa==0 && i_endo_fossa==0){
				f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
				s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
				n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				break;
			}
			else{
				f0 = du[i_epic_fossa].unit();
				n0 = s0.cross(f0);
				break;
			}
		}
		//Antra between pv2 and pv3
		case 6:
		{
				n0 = du[i_epic_antra1_rspv].unit();
				f0 = s0.cross(n0);
				if(i_epic_antra1_rspv==0 && i_endo_antra1_rspv==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				   else if(i_endo_antra1_rspv==90){
					   f0 = du[i_epic_antra1_rspv].unit();
					   n0 = s0.cross(f0);
					}
				break;
		}
		case 16:
		{
				n0 = du[i_epic_antra1_ripv].unit();
				f0 = s0.cross(n0);
				if(i_epic_antra1_ripv==0 && i_endo_antra1_ripv==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				else if(i_endo_antra1_ripv==90){
					   f0 = du[i_epic_antra1_ripv].unit();
					   n0 = s0.cross(f0);
					}
				break;
		}
		//Antra between pv0 and pv1
		case 8:
		{
				n0 = du[i_epic_antra2].unit();
				f0 = s0.cross(n0);
				if(i_epic_antra2==0 && i_endo_antra2==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				break;
		}
		//Lateral
		case 7:
		{
			//if(u[0]>0.5){
				n0 = du[i_epic_lateral].unit();
				f0 = -s0.cross(n0);
				if(i_epic_lateral==0 && i_endo_lateral==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				break;
			/*}
			else if(u[0]<=0.5 && i_endo_lateral!=90){
				n0 = du[i_epic_lateral].unit();
				f0 = s0.cross(n0);
				break;
			}
			else{
				f0 = du[i_endo_lateral].unit();
				n0 = s0.cross(f0);
				break;
			}*/
		}
		//Septum
		case 9:
		{
			/*if(u[0]>0.5){
				n0 = du[i_epic_septum].unit();
				f0 = s0.cross(n0);
				break;
			}
			else if(u[0]<=0.5 && i_endo_septum!=90){*/
				n0 = du[i_endo_septum].unit();
				f0 = s0.cross(n0);
			/*}
			else{
				f0 = du[i_endo_septum].unit();
				n0 = s0.cross(f0);
				break;
			}*/
				if(i_epic_septum==0 && i_endo_septum==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				break;
		}
		//Anterior
		case 10:
		{
			if(u[2]>0.5){
				n0 = du[i_epic_anterior].unit();
				f0 = s0.cross(n0);

				if(i_epic_anterior==0 && i_endo_anterior==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
			}
			else{
				n0 = du[5].unit();
				f0 = s0.cross(n0);
			}
				break;
		}
		//Posterior
		case 11:
		{
				n0 = du[i_epic_posterior].unit();
				f0 = s0.cross(n0);
				if(i_epic_posterior==0 && i_endo_posterior==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				break;
		}
		//septum_bottom_rspv
		case 13:
		{
			//if(u[0]>0.5){
				n0 = du[i_epic_septum_bottom_rspv].unit();
				f0 = -s0.cross(n0);
				if(i_epic_septum_bottom_rspv==0 && i_endo_septum_bottom_rspv==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				break;
			/*}
			else if(u[0]<=0.5 && i_endo_septum_bottom_rspv!=90){
				n0 = du[i_epic_septum_bottom_rspv].unit();
				f0 = s0.cross(n0);
				break;
			}
			else{
				f0 = du[i_endo_septum_bottom_rspv].unit();
				n0 = s0.cross(f0);
				break;
			}*/
		}
		case 17:
		{
				n0 = du[i_epic_septum_bottom_ripv].unit();
				f0 = -s0.cross(n0);
				if(i_epic_septum_bottom_ripv==0 && i_endo_septum_bottom_ripv==0){
					f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
					s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
					n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
				}
				else if(i_endo_septum_bottom_ripv==90){
					   f0 = -du[i_epic_septum_bottom_ripv].unit();
					   n0 = s0.cross(f0);
					}
				break;

		}
			//Anterior_bottom
		case 14:
			{
				//if(u[0]>0.5){
					n0 = du[i_epic_anterior_bottom].unit();
					f0 = s0.cross(n0);
					if(i_epic_anterior_bottom==0 && i_endo_anterior_bottom==0){
						f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
						s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
						n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
					}
					break;
				/*}
				else if(u[0]<=0.5 && i_endo_anterior_bottom!=90){
					n0 = du[i_epic_anterior_bottom].unit();
					f0 = s0.cross(n0);
					break;
				}
				else{
					f0 = du[i_endo_anterior_bottom].unit();
					n0 = s0.cross(f0);
					break;
				}*/
			}
			//Posterior_bottom
			case 15:
			{
				//if(u[0]>0.5){
					n0 = du[i_epic_posterior_bottom].unit();
					f0 = s0.cross(n0);
					if(i_epic_posterior_bottom==0 && i_endo_posterior_bottom==0){
						f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
						s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
						n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
					}
					break;
				/*}
				else if(u[0]<=0.5 && i_endo_posterior_bottom!=90){
					n0 = du[i_epic_posterior_bottom].unit();
					f0 = s0.cross(n0);
					break;
				}
				else{
					f0 = du[i_endo_posterior_bottom].unit();
					n0 = s0.cross(f0);
					break;
				}*/
		   }
		default:
		{
			std::cout << "Case default, element = "<< elem <<"\n";
			/*f0(0) = 1.0; f0(1) = 0.0; f0(2) = 0.0;
			s0(0) = 0.0; s0(1) = 1.0; s0(2) = 0.0;
			n0(0) = 0.0; n0(1) = 0.0; n0(2) = 1.0;*/
			f0(0) = 0.0; f0(1) = 0.0; f0(2) = 0.0;
			s0(0) = 0.0; s0(1) = 0.0; s0(2) = 0.0;
			n0(0) = 0.0; n0(1) = 0.0; n0(2) = 0.0;
			break;
		}
		}
}

f0 = f0.unit();
s0 = s0.unit();
n0 = n0.unit();

/*
if(block_id ==7){
std::cout << "block_id=7 and f0, n0 and s0 are equal to\n";
f0.print();
std::cout<<" ";
n0.print();
std::cout<<" ";
s0.print();
std::cout<<"\n";
}*/
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
            for (auto side : elem->side_index_range())
                if (elem->neighbor_ptr(side) == nullptr)
                {
                    // get boundary info
                    auto boundaryid = mesh.get_boundary_info().boundary_id(elem, side);
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
                        for (int jj = 0; jj < dirichlet_id_points.size(); jj++) {
                            libMesh::Point p_id(qface_point[qp] - dirichlet_id_points[jj]);

                            if ((dirichlet_id_radius[jj]>0) && (dirichlet_id_radius[jj]<47) ){
								is_on_neighborhood = (p_id.norm() <= dirichlet_id_radius[jj]) ? 1 : 0;
							    // imposing BC on sphere
								if (is_on_neighborhood == 1){
										elem_dirichlet_bc=dirichlet_bc_list[jj];
										radius = dirichlet_id_radius[jj];
										break;
								}
							}
                            if(dirichlet_id_radius[jj]==47){
								// imposing BC on set IDs
								 if(i_solution > 1){
									 ExplicitSystem& s1 = es.get_system<ExplicitSystem>("u1");
									 double u1 = s1.point_value(s1.variable_number("u1"), qface_point[qp], *elem);
									 //be careful because elem above is not a surface element, it is a volumetric element
									 is_on_set_fossa = ((u1>anatomic_parameters.thresholdA_fossa)*(u1<anatomic_parameters.thresholdB_fossa))? 1:0;
									 if (is_on_set_fossa){
			                        	 elem_dirichlet_bc=dirichlet_bc_list[jj];
										 radius = dirichlet_id_radius[jj];
										 break;
									 }
								 }
                            }
                            if(dirichlet_id_radius[jj]==48){
								// imposing BC on set IDs
								 if(i_solution > 1){
									 ExplicitSystem& s1 = es.get_system<ExplicitSystem>("u1");
									 double u1 = s1.point_value(s1.variable_number("u1"), qface_point[qp], *elem);
									 //be careful because elem above is not a surface element, it is a volumetric element
									 is_on_set_laa = ((u1>anatomic_parameters.thresholdA_LAA)*(u1<anatomic_parameters.thresholdB_LAA))? 1:0;
									 if (is_on_set_laa){
										 elem_dirichlet_bc=dirichlet_bc_list[jj];
			                        	 radius = dirichlet_id_radius[jj];
										 break;
									 }
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
                        bool cond3 = (found_side_set_in_list && radius==47 && is_on_set_fossa)? 1 : 0;
                        bool cond4 = (found_side_set_in_list && radius==48 && is_on_set_laa)? 1 : 0;


                        // Matrix contribution of the L2 projection.
                        if(cond1 || cond2 || cond3 || cond4) //
                        {
                            for (unsigned int i = 0; i != n_dofs; i++)
                            {
                                for (unsigned int j = 0; j != n_dofs; j++)
                                    Ke(i, j) += JxW_face[qp] * penalty * phi_face[i][qp] * phi_face[j][qp];
                            }
                        }
                        // Right-hand-side contribution of the L2
                        // projection.
                        for (unsigned int i = 0; i != n_dofs; i++) {
                             if(cond1 || cond2 || cond3 || cond4){
                                Fe(i) += JxW_face[qp] * penalty * elem_dirichlet_bc * phi_face[i][qp]; //dirichlet
                             }
                            else
                                Fe(i) += 0;         //neumann
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
        // Read problem solution index and thresholds
        std::string threshold_pv1       = data("pv1"      , " 1, 1, 0, 0, 0");
        std::string threshold_pv2       = data("pv2"      , " 1, 1, 0, 0, 0");
        std::string threshold_pv3       = data("pv3"      , " 1, 1, 0, 0, 0");
        std::string threshold_pv4       = data("pv4"      , " 1, 1, 0, 0, 0");
        std::string threshold_floor     = data("floor"    , " 1, 1, 0, 0, 0");
        std::string threshold_laa       = data("laa"      , " 1, 1, 0, 0, 0");
        std::string threshold_antra1_rspv= data("antra1_rspv"   , " 1, 1, 0, 0, 0");
        std::string threshold_antra1_ripv= data("antra1_ripv"   , " 1, 1, 0, 0, 0");
        std::string threshold_antra2    = data("antra2"   , " 1, 1, 0, 0, 0");
        std::string threshold_septum    = data("septum"   , " 1, 1, 0, 0, 0");
        std::string threshold_septum_bottom_rspv    = data("septum_bottom_rspv"     , " 1, 1, 0, 0, 0");
        std::string threshold_septum_bottom_ripv    = data("septum_bottom_ripv"     , " 1, 1, 0, 0, 0");
        std::string threshold_anterior_bottom  = data("anterior_bottom"   , " 1, 1, 0, 0, 0");
        std::string threshold_posterior_bottom = data("posterior_bottom"  , " 1, 1, 0, 0, 0");
        std::string threshold_anterior  = data("anterior" , " 1, 1, 0, 0, 0");
        std::string threshold_posterior = data("posterior", " 1, 1, 0, 0, 0");
        std::string threshold_carina1   = data("carina1"  , " 1, 1, 0, 0, 0");
        std::string threshold_carina2   = data("carina2"  , " 1, 1, 0, 0, 0");
        std::string threshold_lateral   = data("lateral"  , " 1, 1, 0, 0, 0");
        std::string threshold_strip     = data("strip"    , " 1, 1, 0, 0, 0");
        std::string threshold_ripv     = data("ripv"      , " 0,0 , 0, 0, 0");
        std::string threshold_fossa     = data("fossa"    , " 1, 1, 0, 0, 0");

        //Fiber type is Milano's paper?
        //https://doi.org/10.1016/j.cma.2020.113468
        std::string milano_aux = data("milano", " 0");

        AnatomicalParameters anatomic_parameters;

        anatomic_parameters.reinit(threshold_pv1, threshold_pv2, threshold_pv3, threshold_pv4,
        							threshold_floor, threshold_laa, threshold_antra1_rspv,threshold_antra1_ripv, threshold_antra2,
    								threshold_septum, threshold_anterior, threshold_posterior,
    								threshold_carina1, threshold_carina2, threshold_lateral, threshold_strip, threshold_ripv,
    								threshold_septum_bottom_rspv, threshold_septum_bottom_ripv, threshold_anterior_bottom,
    								threshold_posterior_bottom, threshold_fossa, milano_aux);


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
        // Use the Importer to read
        //importer.read(mesh_name);
        // If we did not use the importer we would have done
        //mesh.read(mesh_name);

        // We may want to refine this mesh n times
        // if the original mesh is too coarse
        // read the number of times we will refine the mesh, 0 by default
        // 0 does not refine
        // refs =3 in the input file will make it refine 3 times
        /*int num_refs = data("refs", 0);
        std::cout << "Refining the mesh " << num_refs << " times. " << std::endl;
        // Refine the mesh
        MeshRefinement(mesh).uniformly_refine(num_refs);
        // Finalize the mesh to run simulations
        mesh.prepare_for_use();*/
        int n_refinements = data("refs", 0);
        std::cout << "n_refs: " << n_refinements << std::endl;
        BeatIt::serial_mesh_partition(init.comm(), mesh_name, &mesh, n_refinements);
        
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
        equation_systems.get_system("Poisson").add_variable("u", SECOND); // QUAD9 ELEMENTS FOR SQUARE
    }

    // Initialize the data structures for the equation system.
    equation_systems.init();

    // Prints information about the system to the screen.
    equation_systems.print_info();

    // Reading the input for each of the problems
    for (int i = 0; i < N; i++) {

        std::string ui = "u" + std::to_string(i);    //how cool! to use i in the name here!!!
        ExplicitSystem& si = equation_systems.add_system<ExplicitSystem>(ui);
        si.add_variable(ui, FIRST);
        si.init();

 	   //READ INPUT I > BE CAREFUL BECAUSE HERE I WILL BE OVERWRITTING WHAT I HAD BEFORE
 	   poisson_solver.read_input_i(i, equation_systems);

       // reinit Poisson member variables using input file values
 	  poisson_solver.reinit(poisson_solver.dirichlet_x_list, poisson_solver.dirichlet_y_list,  poisson_solver.dirichlet_z_list,
 			  poisson_solver.dirichlet_r_list,  poisson_solver.dirichlet_side_set_list_aux, poisson_solver.dirichlet_bc_list_aux);

 	  poisson_solver.assemble_poisson(equation_systems, "Poisson", i, anatomic_parameters);

        //  Libmesh zeros out matrix and RHS when i call solve
        // preventing that to happen:
        equation_systems.get_system<ImplicitSystem>("Poisson").zero_out_matrix_and_rhs=false;

        equation_systems.get_system("Poisson").solve();

        // Copy the content of "u" (solution to Poisson problem) to "ui"
        *si.solution = *equation_systems.get_system("Poisson").solution; // "unique vector"
        si.update(); //distributes/copies the solution of the global on the shared interface > local processor vectors - "repeated vector"
    }

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

    const DofMap& dofmap = system.get_dof_map(); // replicating assembly function

    std::cout << " Looping over elements to set up subdomain IDs"<< std::endl;

    for (const auto& elem: mesh.active_local_element_ptr_range()) {
        int blockid = elem->subdomain_id();
        libMesh::Point centroid = elem->centroid();
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
    save_fibers(test, output_file_name, output_number, equation_systems, mesh);



    // Information about the time, such as
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
    libMesh::EquationSystems es(mesh);


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
    BeatIt::ElectroSolver* solver = BeatIt::ElectroSolver::ElectroFactory::Create(model, es);
    // Setup the EP model using the input file
    std::cout << "Calling setup..." << std::endl;
    solver->setup(data, model);
    // Initialize systems
    std::cout << "Calling init ..." << std::endl;
    es.print_info();
    // Set up initial conditions at time
    solver->init(datatime.M_startTime);
    // Now read the fibers if wanted
    bool read_fibers = data("read_fibers", false);
    if(read_fibers)
    {
        // First show the elemental variables that can be imported
        auto elemental_variables = importer.get_elem_var_names();
        for (auto && var : elemental_variables) std::cout << var << std::endl;
        solver->read_fibers(importer,1); //solver->fiber_system?
    }
    // Export simulation parameters
    // This will also export the fiber field
//    solver->save_parameters();
    // Assemble matrices
    std::cout << "Assembling matrices" << std::endl;
    solver->assemble_matrices(datatime.M_dt);
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

    // The end
    std::cout << "Good luck with your simulation :P" << std::endl;
    return 0;
}
