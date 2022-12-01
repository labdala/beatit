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

int main(int argc, char ** argv)
{
    std::cout << "HELLO!" <<std::endl;
    // Initialize libmesh
    using namespace libMesh;
    // Import input file
    GetPot data = BeatIt::readInputFile(argc, argv);

    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    // Create mesh:
    // one for the conforming mesh
    // one with the regions to be mapped
    libMesh::Mesh mesh(init.comm());
    //libMesh::DistributedMesh mesh(init.comm()); // does what the serial_mesh_partition is doing?
    //libMesh::ParallelMesh mesh(init.comm());
    //mesh.allow_renumbering(false);


	mesh.read("../example_laryssa_define_subregions/bin/03meshLA.e", nullptr, true);
	//BeatIt::serial_mesh_partition(init.comm(), "./bin/03meshLA.e", &mesh, 0);
	//init.comm().barrier();

	// set tolerance
	double bc_tolerance = data("bc_tolerance", 1e-1); // default = 1e-2
	std::cout<< "The tolerance for finding centroid is " << bc_tolerance <<std::endl;

    //init.comm().barrier();

	std::vector<std::unique_ptr<libMesh::Mesh>> subregions_mesh_ptr(1);
//	std::vector<std::unique_ptr<libMesh::ReplicatedMesh>> subregions_mesh_ptr(1);

	std::cout << "Reading mvr PTR mesh " << std::endl;
	//	subregions_mesh_ptr[0].reset(new libMesh::Mesh(init.comm()));

	init.comm().barrier();

	std::cout<< "THIS IS NEW" << std::endl;

	std::cout << "Reading endo ptr mesh " << std::endl;
	subregions_mesh_ptr[0].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[0] -> read("../example_laryssa_define_subregions/bin/LIPV.e");

	/*// Copying each subregion mesh to each processor
	subregions_mesh_ptr[1].reset(new libMesh::Mesh(init.comm()));
	//subregions_mesh_ptr[0].reset(new libMesh::ReplicatedMesh(init.comm()));
	subregions_mesh_ptr[1] -> read("./bin/LIPV.e");

	std::cout << "Reading mvr PTR mesh " << std::endl;
	subregions_mesh_ptr[2].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[2] -> read("./bin/LSPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[3].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[3] -> read("./bin/RSPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[4].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[4] -> read("./bin/RIPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[5].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[5] -> read("./bin/MVR.e");
*/
	// IDs for each region:
	// endo=2, epi=3, MV=4
	// LIPV=5, LSPV=6, RSPV=7, RIPV=8
	std::string subregions_ids = "5";
	std::vector < unsigned int > subregions_ids_vec;
	BeatIt::readList(subregions_ids, subregions_ids_vec);
	int default_id = 3;

	// subregions
	int num_bc_subregions = subregions_mesh_ptr.size();
	int num_ids = subregions_ids_vec.size();

	std::vector < std::unique_ptr<libMesh::PointLocatorTree> > subregions_locator_ptr(num_bc_subregions);
	for (int i = 0; i < num_bc_subregions; ++i)
	{
		subregions_mesh_ptr[i]->print_info();

		std::cout<< " reset point locator tree" <<std::endl;
		subregions_locator_ptr[i].reset(new libMesh::PointLocatorTree(*subregions_mesh_ptr[i]));
		std::cout<< " I DO NOT KNOW WHAT THIS IS FOR: set point close" <<std::endl;
		subregions_locator_ptr[i]->set_close_to_point_tol(bc_tolerance);

		// this gives us more information of what is happening - maybe it is doing a linear search through the mesh (not n*log(n) anymore, but rather n^2)
		subregions_locator_ptr[i]->_verbose=true;
	}

	libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	el = mesh.active_local_elements_begin();
	int bc_missing_elements = 0;
	std::set<int> bc_missing_elem;     // why do we use a set and not a vector here??

	std::cout << "Looping over elements of the original mesh: " << std::endl;


	int c = 0;
	for (; el != end_el; ++el)
	{
		Elem * elem = *el;

		for (int s = 0; s < elem->n_sides(); ++s)
		{
			// check if element is on the surface of the mesh
			if (nullptr == elem->neighbor_ptr(s))
			{
				// calculate centroid of the surface element (opposed to volumetric element)
				auto centroid = elem->side_ptr(s)->centroid();
				bool found = false;

				// Let's loop over the number of subregion meshes and look for the centroid
				// As soon as we find it, we can pass to the next element
				for (int i = 0; i < num_bc_subregions; ++i)
				{

					// see if centroid is part of that region

					//CHECK THIS OUT 	 -- IS IT FALSE OR NULL IF IT DOES NOT FIND THE CENTROID WITH THE TOLERANCE WE GAVE?
					const Elem* subregion_element = (*subregions_locator_ptr[i])(centroid);
					if (subregion_element)
					{
						unsigned int sideset_id = i+1;
						if(num_ids == num_bc_subregions) sideset_id = subregions_ids_vec[i];
						mesh.boundary_info->add_side (elem, s, sideset_id);
						found = true;
						break;
					}
				}

				if (!found)
				{
					mesh.boundary_info->add_side (elem, s, default_id);
					int elID = elem->id();
					bc_missing_elements++;
					bc_missing_elem.insert(elID);
				}

			}
		}
	}
	std::cout << "Missed " << bc_missing_elem.size() << " boundary elements." << std::endl;

	mesh.prepare_for_use();

    std::cout << "Export to new mesh" << std::endl;
    std::string output_file = data("output", "remapped_mesh_rpv_lpv_mvr.e");
    ExodusII_IO exporter(mesh);
    exporter.write(output_file);

    std::cout << "Good luck!" << std::endl;
    return 0;
}

