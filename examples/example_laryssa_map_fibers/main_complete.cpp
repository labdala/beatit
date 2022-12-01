/*
 * main_laryssa_define_subregions
 * based and inspired in main_remap_regions.cpp by srossi
 *
 *  Created on: Oct 25, 2021
 *      Author: labdala
 */

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include <libmesh/point_locator_tree.h>
#include "libmesh/elem.h"
#include "libmesh/getpot.h"
#include <algorithm>
#include <tuple>
#include "Util/IO/io.hpp"
#include <libmesh/boundary_mesh.h>
#include <libmesh/parallel_mesh.h>
#include <libmesh/replicated_mesh.h>
#include "libmesh/point_locator_nanoflann.h"
#include "libmesh/mesh_base.h"
#include "libmesh/boundary_info.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"
#include <libmesh/diff_context.h>
#include <libmesh/numeric_vector.h>

int main(int argc, char ** argv)
{

    // Initialize libmesh
    using namespace libMesh;
    // Import input file
    GetPot data = BeatIt::readInputFile(argc, argv);

    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // Mesh that has fiber data, multiple blocks and duplicate nodes
    libMesh::Mesh mesh_fiber(init.comm());

    // Mesh that has sidesets assigned and one blockID=1
    libMesh::Mesh mesh_sideset(init.comm());
    
    // Use exodus importer so we can read the fibers
    libMesh::ExodusII_IO importer(mesh_fiber);

    // Get info from input file
    std::string fiber_file = data("fiber_file", "NONE");
    std::string sideset_file = data("sideset_file", "NONE");


    // Read meshes
    if ("NONE" != fiber_file && "NONE" != sideset_file)
    {
        std::cout << "Reading mesh with fiber =  " << fiber_file << " and mapping fibers to = " << sideset_file << std::endl;
        mesh_sideset.read(sideset_file);
        importer.read(fiber_file);
        mesh_fiber.read(fiber_file);
    }
    else
    {
        std::cout << "Set input_mesh in input file!" << std::endl;
        throw std::runtime_error("No mesh given");
    }

    //mesh_fiber.prepare_for_use();
    mesh_fiber.print_info();


    // Move mesh_fiber to the a vector pointer
    //std::cout << " Creating locator ptr to point to mesh - HERE WE USE STD::MOVE" << fiber_file << std::endl; 
	//std::vector<std::unique_ptr<libMesh::Mesh>> fiber_mesh_ptr(1);
	//fiber_mesh_ptr[0].reset(new libMesh::Mesh(init.comm()));
	//fiber_mesh_ptr[0] -> read(fiber_file);
	//*fiber_mesh_ptr[0] =std::move(mesh_fiber);

    // Set tolerance for distance for point locator
	double bc_tolerance = data("bc_tolerance", 1e-4);
	std::cout<< "The tolerance for finding centroid is " << bc_tolerance <<std::endl;

    // Create point locator
	std::vector < std::unique_ptr<libMesh::PointLocatorBase> > fiber_locator_ptr(1);
//	fiber_mesh_ptr[0]->print_info();
	std::cout<< " reset point locator tree" <<std::endl;
	fiber_locator_ptr[0].reset(new libMesh::PointLocatorNanoflann(mesh_fiber));
	fiber_locator_ptr[0]->set_close_to_point_tol(bc_tolerance);
	fiber_locator_ptr[0]->_verbose=true;
	

    // Create fiber system for the mesh where we have the fiber field data
    libMesh::EquationSystems  es_fiber(mesh_fiber);
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& fiber_system = es_fiber.add_system<FiberSystem>("fibers");
    fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.init();
    auto& f_v = fiber_system.solution;

    std::cout << " +++++++++++++ EquationSystems setup in the mesh where we have fiber info \n";
    es_fiber.print_info();

    unsigned int step=1;
    std::cout << "Available variables in exodus file" << std::endl;
    auto& elem_var_names = importer.get_elem_var_names();
    for (auto && var : elem_var_names) std::cout << var << std::endl;

    
    std::vector<std::string> system_variables;
    system_variables.push_back("fibersx");
    system_variables.push_back("fibersy");
    system_variables.push_back("fibersz");

    //Read fiber data into mesh_fiber
    importer.copy_elemental_solution(fiber_system, "fibersx", "fibersx", step);
    importer.copy_elemental_solution(fiber_system, "fibersy", "fibersy", step);
    importer.copy_elemental_solution(fiber_system, "fibersz", "fibersz", step);
    std::cout << " +++++++++++++ EquationSystems in the mesh where we have fiber info after adding fiber data \n";
    fiber_system.print_info();
	
   // mesh_fiber.prepare_for_use();
    std::cout << "norm of fiber solution = " << fiber_system.solution->l2_norm() << std::endl;
    std::cout << "norm of fiber current local solution = " << fiber_system.current_local_solution->l2_norm() << std::endl;

    // Set up fFiber system in the mesh where we DO NOT have fiber field data
    libMesh::EquationSystems  es_sideset(mesh_sideset);
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& fiber_system_sideset = es_sideset.add_system<FiberSystem>("fibers");
    fiber_system_sideset.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system_sideset.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system_sideset.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system_sideset.init();

    std::cout << "----------------- EquationSystems setup in the mesh where we DO NOT have fiber info \n";
    es_sideset.print_info();
    std::cout << "----------------- \n";

	
    //mesh_sideset.prepare_for_use();
    
    // variable to store number of elements we cannot find
    int bc_missing_elements = 0;
    std::set<int> bc_missing_elem;

	// Loop over elements of mesh_sideset
	//libMesh::MeshBase::const_element_iterator el = mesh_sideset.active_local_elements_begin();
	//const libMesh::MeshBase::const_element_iterator end_el = mesh_sideset.active_local_elements_end();
	//libMesh::MeshBase::const_element_iterator el = mesh_fiber.active_local_elements_begin();
	//const libMesh::MeshBase::const_element_iterator end_el = mesh_fiber.active_local_elements_end();

	libMesh::MeshBase::const_element_iterator el = mesh_fiber.elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh_fiber.elements_end();

    std::cout << "Looping over elements of the fiber mesh: " << std::endl;
    //std::cout << "Looping over elements of the sideset mesh: " << std::endl;
/*	int c = 0;
	for (; el != end_el; ++el)
	{
		Elem * elem = *el;
        elem->subdomain_id()=1;    	// Calculate centroid of the element in the mesh_sideset
    	// Calculate centroid of the element in the mesh_sideset
		bool found = false;

	    // Check wether centroid is part of the mesh_fiber
    	const Elem* fiber_elem = (*fiber_locator_ptr[0])(centroid);
		if (fiber_elem)
		{
			found = true;
		    std::cout << "found element id=" << elem->id() << " and its correspondent in the fiber mesh " << fiber_elem->id() << std::endl;
            std::cout << "-----------------Elem sideset -----------" << std::endl;
            elem->print_info();
            std::cout << "-----------------Elem fibers -----------" << std::endl;
            fiber_elem->print_info();
            
            std::vector<double> ff(3);

            // Copy fiber info into ff vector
           // ff[0] = fiber_system.point_value(fiber_system.variable_number("fibersx"), centroid, *fiber_elem);
           // ff[1] = fiber_system.point_value(fiber_system.variable_number("fibersy"), centroid, *fiber_elem);
           // ff[2] = fiber_system.point_value(fiber_system.variable_number("fibersz"), centroid, *fiber_elem);
           // std::cout << "ff = (" << ff[0] << ", " << ff[1] << ", " << ff[2] <<") ; norm=" << std::sqrt(ff[0]*ff[0]+ff[1]*ff[1]+ff[2]*ff[2])  << std::endl;
			 break;
		} 
		if (!found)
		{
				int elID = elem->id();
				bc_missing_elements++;
				bc_missing_elem.insert(elID);
		}
        break;
	}*/
	
//	std::cout << "Missed " << bc_missing_elem.size() << " elements." << std::endl;

/*	el = mesh_fiber.active_local_elements_begin();
	for (; el != end_el; ++el)
	{
		Elem * elem = *el;
        if(elem->subdomain_id()!=1) std::cout << "elem block id = " << elem->subdomain_id() << std::endl;
    }
	mesh_fiber.prepare_for_use();
*/
    std::cout << "Export to new mesh" << std::endl;
    std::string output_file = data("output", "output_centroid.e");
    //ExodusII_IO exporter(mesh_sideset);
    ExodusII_IO exporter(mesh_fiber);
    exporter.write(output_file);
    exporter.write_equation_systems(output_file, es_fiber);
    exporter.write_element_data(es_fiber);

    std::cout << "Good luck!" << std::endl;

    //test
    //myfile.close();

    return 0;
}
