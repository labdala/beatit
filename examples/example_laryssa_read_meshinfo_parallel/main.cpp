/*
 * main_laryssa_read_meshinfo_parallel
 *
 *  Created on: Nov 26, 2021
 *      Author: labdala
 *
 *      Goal is to read mesh and fibers from an input file in parallel
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
#include "libmesh/equation_systems.h"
#include "libmesh/explicit_system.h"
#include "libmesh/communicator.h"
#include "libmesh/vtk_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include <fstream>
#include <sstream>
#include <string>


void read_fibers(libMesh::EquationSystems & equation_systems, libMesh::ReplicatedMesh * mesh_ptr, libMesh::Parallel::Communicator & rank_split_comm, libMesh::ExodusII_IO * importer_ptr)
{
	//std::cout << "Inside read_fibers function" << std::endl;
    libMesh::ExplicitSystem & system = equation_systems.add_system<libMesh::ExplicitSystem>("fibers");
    system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    equation_systems.init();
    equation_systems.print_info();

    auto& elem_var_names = importer_ptr->get_elem_var_names();
    std::cout << "Names in exodus file: " << std::endl;
    for (auto && v:elem_var_names){
    	std::cout << v << std::endl;
    }

    // copy_elemental_solution(system, system_variable_name, exodus_variable_name)
	importer_ptr->copy_elemental_solution(system, "fibersx", "fibersX");
	importer_ptr->copy_elemental_solution(system, "fibersy", "fibersY");
	importer_ptr->copy_elemental_solution(system, "fibersz", "fibersZ");
	double norm = system.solution->linfty_norm();
	std::cout << "Fiber norm = " << norm << std::endl;
}


int main(int argc, char ** argv)
{
	using namespace libMesh;

	// Initialization and processors
	LibMeshInit init(argc, argv, MPI_COMM_WORLD);
    int np = init.comm().size();
    std::cout << "reading mesh to run with " << np << " processors" << std::endl;
	const unsigned int rank = init.comm().rank();
	std::cout << "This is processor number: " << rank << std::endl;

	// Create serial mesh
	GetPot data = BeatIt::readInputFile(argc, argv);
	std::string mesh_name = data("mesh", " ");
	std::cout << "Mesh name = " << mesh_name << std::endl;
    std::cout << "Create replicated mesh ptr" << std::endl;
	ReplicatedMesh * original_mesh_ptr;
	ReplicatedMesh * mesh_to_split;

	// Use this communicator if I am doing something in a single processor
	std::cout << "Splitting communicator" << std::endl;
    libMesh::Parallel::Communicator rank_split_comm;
    init.comm().split(rank, rank, rank_split_comm);

	std::cout << "Create importer object" << std::endl;
    libMesh::ExodusII_IO * importer_ptr;

	// checkpoint output folder name
	std::string outputname = "03meshLA.cpr";
	// output in binary?
	bool binary = data("binary", true);

	if(rank ==0) {
		// need to have another mesh here otherwise get MPI_Abort
		mesh_to_split = new ReplicatedMesh(rank_split_comm);
		mesh_to_split->allow_renumbering(false);
		mesh_to_split->read(mesh_name);
		processor_id_type n_procs = np;
		MetisPartitioner partitioner;
		partitioner.partition(* mesh_to_split, n_procs);

		std::cout << "Splitting mesh" << std::endl;
		CheckpointIO cpr(* mesh_to_split);
		cpr.current_processor_ids().clear();
		for (unsigned int i = 0; i < n_procs; i++)
			cpr.current_processor_ids().push_back(i);
		cpr.current_n_processors() = n_procs;
		cpr.parallel() = true;
		cpr.binary() = binary;
		std::cout << outputname << std::endl;
		cpr.write(outputname);
	}
	init.comm().barrier();
	std::cout << "Done reading mesh, processor=" << rank << std::endl;

	std::cout << "create parallel mesh , processor = " << rank << std::endl;
	ParallelMesh parallel_mesh(init.comm());

	// Each processor reads its own piece
	parallel_mesh.allow_renumbering(false);
	parallel_mesh.read(outputname);
	parallel_mesh.print_info();

	std::cout << "Create equation systems for parallel mesh" << std::endl;
	EquationSystems parallel_es(parallel_mesh);

	ExplicitSystem & parallel_system = parallel_es.add_system<ExplicitSystem>("parallel_fs");
	std::cout << "Add variables" << std::endl;
	parallel_system.add_variable("fibersx", CONSTANT, MONOMIAL);
	parallel_system.add_variable("fibersy", CONSTANT, MONOMIAL);
	parallel_system.add_variable("fibersz", CONSTANT, MONOMIAL);
	std::cout << "Init System" << std::endl;
	parallel_es.init();
	// Prints information about the system to the screen.
	parallel_es.print_info();


	/*if(rank ==0) {

		original_mesh_ptr = new Mesh(rank_split_comm); // or replicatedmesh here?

		std::cout << "Create importer ptr" << std::endl;
		importer_ptr = new ExodusII_IO(* original_mesh_ptr);
		std::cout << "Read mesh using importer ptr" << std::endl;
		importer_ptr->read(mesh_name);
		original_mesh_ptr->print_info();


		std::cout << "Reading fibers on processor " << rank << std::endl;
		EquationSystems replicated_es(* original_mesh_ptr);
		read_fibers(replicated_es, original_mesh_ptr, rank_split_comm, importer_ptr);

		// read_fibers_inline
		{
			        libMesh::ExplicitSystem & system = replicated_es.add_system<libMesh::ExplicitSystem>("fibers");
			        system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
			        system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
			        system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
			        replicated_es.init();
			        replicated_es.print_info();

			        //auto& elem_var_names = importer_ptr->get_elem_var_names();
			        //std::cout << "Names in exodus file: " << std::endl;
			        //for (auto && v:elem_var_names){
			        //	std::cout << v << std::endl;
			        //}

					importer_ptr->copy_elemental_solution(system, "fibersx", "fibersx");
					importer_ptr->copy_elemental_solution(system, "fibersy", "fibersx", 1);
					importer_ptr->copy_elemental_solution(system, "fibersz", "fibersz", 1);
					double norm = system.solution->linfty_norm();
					std::cout << "Fiber norm = " << norm << std::endl;
		}
	}
*/


	bool create_mesh_info = data("create_mesh_info", false);
	int nel = 0;
	if (rank==0){

		original_mesh_ptr = new Mesh(rank_split_comm); // or replicatedmesh here?
		original_mesh_ptr->allow_renumbering(false);

		std::cout << "Create importer ptr" << std::endl;
		importer_ptr = new ExodusII_IO(* original_mesh_ptr);
		std::cout << "Read mesh using importer ptr" << std::endl;
		importer_ptr->read(mesh_name);
		std::cout << "This is inside if condition. Info about original_mesh" << std::endl;
		original_mesh_ptr->print_info();

		EquationSystems replicated_es(* original_mesh_ptr);
		read_fibers(replicated_es, original_mesh_ptr, rank_split_comm, importer_ptr);
		libMesh::ExplicitSystem & fiber_system=replicated_es.get_system<ExplicitSystem>("fibers");
	    auto & f_v= fiber_system.solution;

	    if(create_mesh_info){
	    	std::ofstream mesh_info;
	    	mesh_info.open("mesh_info.txt");
			nel = original_mesh_ptr->n_elem();
			std::cout << "Looping over all the " << nel << " elements of mesh" << std::endl;
			for (int c = 0; c < nel; ++c)
			{
					bool found = false;
					Elem * elem = original_mesh_ptr->elem_ptr(c);
					const unsigned int elem_id = elem->id();
					Point centroid;
					double fx, fy, fz;
					std::vector<libMesh::dof_id_type> dof_indices;


					//if(rank==0){

					centroid = elem->vertex_average();
					auto & test = fiber_system.get_dof_map();
					test.dof_indices(elem, dof_indices);
					f_v->set(dof_indices[0], fx); // check get_dof_map
					f_v->set(dof_indices[1], fy);
					f_v->set(dof_indices[2], fx);
					Point fiber_el(fx, fy, fz);

					//}
					// can I use get_dof_map from serial to find indexes of elements in parallel?
					// serial mesh still partitions the mesh

					mesh_info << elem_id << " ";
					mesh_info << centroid(0) << " ";
					mesh_info << centroid(1) << " ";
					mesh_info << centroid(2) << " ";
					mesh_info << fiber_el(0) << " ";
					mesh_info << fiber_el(1) << " ";
					mesh_info << fiber_el(2) << " ";
					mesh_info << "\n";

					if(false){
					// Broadcast centroid    - below I actually wanted communicator from processor 0
					//double centroidx = centroid(0);
					//double centroidy = centroid(1);
					//double centroidz = centroid(2);
					//original_mesh_ptr->comm().broadcast<double>(centroidx, rank);
					//original_mesh_ptr->comm().broadcast<double>(centroidy, rank);
					//original_mesh_ptr->comm().broadcast<double>(centroidz, rank);


					//original_mesh_ptr->comm().broadcast<double>(fiber_el(0), rank);
					//original_mesh_ptr->comm().broadcast<double>(fiber_el(1), rank);
					//original_mesh_ptr->broadcast<double>(fiber_el(2), rank);


					// if (found == false) throw error
					//for (i = 0; i < n_procs; i++){
					// check if element is in one of the processors of parallel mesh;
					// >>>> const Elem* rank_element = (*SOMETHING[i])(centroid);    // this is using nanoflann

					//if(rank_element){
					//	assign fiber from elem in replicated_es to elem in parallel_es;
					//}

					//}
					}
				}
			mesh_info.close();
	    }
	}
	 // close set <<<<<<<<<<<<<<<<<

	if(create_mesh_info==false)
	{
		// Looping over mesh_info.txt lines
		std::ifstream infile("mesh_info.txt");

		std::cout << "Looping over mesh_info.txt lines" << std::endl;
		std::string line;
		std::getline(infile, line);
		while (std::getline(infile, line))
		{
			std::istringstream iss(line);
			std::cout << "line = " << line << std::endl;
			int a;
			double ca, cb, cc, fa, fb, fc;
			if (!(iss >> a >> ca >> cb >> cc >> fa >> fb >> fc)) { break; } // error

			Elem * el_in_processor = parallel_mesh.query_elem_ptr(a);
			if(el_in_processor){
				Point centroid_serial = el_in_processor->vertex_average();
				// check if centroid is the same in parallel and serial meshes
				std::cout << "a = " << a << ", centroid=(" << ca << ", "<< cb << ", " << cc << ")" << std::endl;
				std::cout << "centroid in serial=" <<centroid_serial(0) << std::endl;
				/*if(centroid_serial(0)!=ca || centroid_serial(1)!=cb || centroid_serial(2)!=cc );
					std::cout <<  "they are equal" << std::endl;*/
			}
		}
	}

	//VTKIO(parallel_mesh).write("03meshLA_parallel.pvtu");


	std::cout << "All done! :)" << std::endl;
    return 0;

}

