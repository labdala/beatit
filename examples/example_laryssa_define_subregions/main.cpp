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

/*int main(int argc, char ** argv)
{
    using namespace libMesh;
    GetPot data = BeatIt::readInputFile(argc, argv);
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);
    libMesh::Mesh mesh(init.comm());
    mesh.read("F72_m0_pv4_ext.e");
    std::string name_subdomain = "pv4";
    mesh.subdomain_name(1)=name_subdomain;
    std::string output_file = data("output", "F72_m0_pv4_ext.e");
    ExodusII_IO exporter(mesh);
    exporter.write(output_file);
    return 0;
}*/
/*
 * main_laryssa_define_subregions
 * based and inspired in main_remap_regions.cpp by srossi
 *
 *  Created on: Oct 25, 2021
 *      Author: labdala
 */

int main(int argc, char ** argv)
{

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
    libMesh::ExodusII_IO importer(mesh);


    std::cout << "Reading mesh" << std::endl;  

    mesh.read("/nas/longleaf/home/laryssa/beatit/build-tmp-1.7.0/examples/example_laryssa_define_subregions/F72_fsi_fibers.e");
    mesh.print_info();
    mesh.subdomain_name(1) = "left_atrium";
//    mesh.get_boundary_info().sideset_name(4)="mvr";

/*	libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	el = mesh.active_local_elements_begin();*/
//	std::string name_subdomain = "mvr";
/*	for (; el != end_el; ++el)
	{
		Elem * elem = * el; */
//      Name block id 
//		elem->subdomain_id()=name_subdomain;

//      Name side set
/*        for (auto side : elem->side_index_range()){
                 if (elem->neighbor_ptr(side) == nullptr)
                 {
                   // get boundary info
                     unsigned int n_boundary_ids=mesh.boundary_info->n_boundary_ids(elem,side);
            	     //std::cout <<* elem << std::endl;
                     //if(n_boundary_ids !=0) std::cout << "n_b_ids=" << n_boundary_ids << std::endl;
                     std::vector<short int> boundary_ids_vec(n_boundary_ids);
                     mesh.boundary_info->boundary_ids(elem,side, boundary_ids_vec);
                     short int boundaryid=boundary_ids_vec[0];
                     std::cout << boundaryid << std::endl;  
                  // if(boundaryid==4)     std::cout << "hi" << "\n";
                   }
        }*/
//    }
    std::cout << "Export mesh" << std::endl;
    std::string output_file = data("output", "F72_fsi_fibers_named.e");
    ExodusII_IO exporter(mesh);
    exporter.write(output_file);

    return 0;

/*

	if(false){


	//mesh.read("./bin-1.7.0/03meshLA.e", nullptr, true);
	//mesh.read("./bin-1.7.0/sphere.e", nullptr, true);
	//BeatIt::serial_mesh_partition(init.comm(), "./bin-1.7.0/03meshLA.e", &mesh, 0);

	//test
	std::ofstream myfile;
	myfile.open("centroid.txt");


	// set tolerance
	double bc_tolerance = data("bc_tolerance", 1e-2); // default = 1e-2
	std::cout<< "The tolerance for finding centroid is " << bc_tolerance <<std::endl;


	std::vector<std::unique_ptr<libMesh::Mesh>> subregions_mesh_ptr(7);
//	std::vector<std::unique_ptr<libMesh::ReplicatedMesh>> subregions_mesh_ptr(1);


	// Copying each subregion mesh to each processor
	subregions_mesh_ptr[0].reset(new libMesh::Mesh(init.comm()));
	//subregions_mesh_ptr[0].reset(new libMesh::ReplicatedMesh(init.comm()));
	subregions_mesh_ptr[0] -> read("./bin-1.7.0/LIPV.e");
	//subregions_mesh_ptr[0] -> read("./bin-1.7.0/MVR_paraview.e");
	//subregions_mesh_ptr[0] -> read("./bin-1.7.0/sphere_part1.e");

	std::cout << "Reading mvr PTR mesh " << std::endl;
	subregions_mesh_ptr[1].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[1] -> read("./bin-1.7.0/LSPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[2].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[2] -> read("./bin-1.7.0/RSPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[3].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[3] -> read("./bin-1.7.0/RIPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[4].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[4] -> read("./bin-1.7.0/MVR.e");

	std::cout << "Reading endo ptr mesh " << std::endl;
	subregions_mesh_ptr[5].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[5] -> read("./bin-1.7.0/endo.e");

	std::cout << "Reading endo ptr mesh " << std::endl;
	subregions_mesh_ptr[6].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[6] -> read("./bin-1.7.0/epi.e");

	// IDs for each region:
	// endo=2, epi=3, MV=4
	// LIPV=5, LSPV=6, RSPV=7, RIPV=8
	std::string subregions_ids = " 5, 6, 7, 8, 4, 2, 3";
	//std::string subregions_ids = " 4";
	std::vector < unsigned int > subregions_ids_vec;
	BeatIt::readList(subregions_ids, subregions_ids_vec);
	int default_id = 9;

	// subregions
	int num_bc_subregions = subregions_mesh_ptr.size();
	int num_ids = subregions_ids_vec.size();

	std::vector < std::unique_ptr<libMesh::PointLocatorBase> > subregions_locator_ptr(num_bc_subregions);
	for (int i = 0; i < num_bc_subregions; ++i)
	{
		subregions_mesh_ptr[i]->print_info();

		std::cout<< " reset point locator tree" <<std::endl;
		subregions_locator_ptr[i].reset(new libMesh::PointLocatorTree(*subregions_mesh_ptr[i]));
		std::cout<< " I DO NOT KNOW WHAT THIS IS FOR: set point close" <<std::endl;
		subregions_locator_ptr[i]->set_close_to_point_tol(bc_tolerance);
		//subregions_locator_ptr[i]->enable_out_of_mesh_mode();
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
				auto centroid = elem->side_ptr(s)->true_centroid();
				bool found = false;

				//myfile << "centroid = " << centroid << " \n";

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
						//myfile << "found = " << found << " \n";
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
    std::string output_file = data("output", "output_centroid.e");
    ExodusII_IO exporter(mesh);
    exporter.write(output_file);

    std::cout << "Good luck!" << std::endl;

    //test
    //myfile.close();

    return 0;

	}

*/
	//NEW VERSION
	//------------NEW IMPLEMENTATION
/*
	// Initialize libmesh
	using namespace libMesh;
	// Import input file
	GetPot data = BeatIt::readInputFile(argc, argv);

	// Initialize MPI
	LibMeshInit init(argc, argv, MPI_COMM_WORLD);

	// Create mesh:
	std::cout << "Reading mesh" << std::endl;
	libMesh::Mesh mesh(init.comm());
	mesh.read("./bin-1.7.0/03meshLA.e", nullptr, true);
	//mesh.read("./bin-1.7.0/03meshLA_nodeset.e", nullptr, true);
	//mesh.read("./bin-1.7.0/mv.e", nullptr, true);


	std::ofstream nodes_found;
	nodes_found.open("nodes_found.txt");

	std::ofstream nodes_subregion;
	nodes_subregion.open("nodes_subregion.txt");

	std::ofstream nodes_mesh;
	nodes_mesh.open("nodes_mesh.txt");


	std::cout << "WARNING! This code assumes each mesh node has only one ID associated with it" << std::endl;
*/
//if(false){
	// set tolerance
/*	double bc_tolerance = data("bc_tolerance", 1e-2);


	std::vector<std::unique_ptr<libMesh::Mesh>> subregions_mesh_ptr(1);


	// Copying each subregion mesh to each processor
	subregions_mesh_ptr[0].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[0] -> read("./bin-1.7.0/mv_part1.e");
*/
	/*std::cout << "Reading mvr PTR mesh " << std::endl;
	subregions_mesh_ptr[1].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[1] -> read("./bin-1.7.0/LIPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[2].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[2] -> read("./bin-1.7.0/LSPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[3].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[3] -> read("./bin-1.7.0/RSPV.e");

	std::cout << "Reading rspv PTR mesh " << std::endl;
	subregions_mesh_ptr[4].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[4] -> read("./bin-1.7.0/RIPV.e");

	std::cout << "Reading endo ptr mesh " << std::endl;
	subregions_mesh_ptr[5].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[5] -> read("./bin-1.7.0/endo.e");

	std::cout << "Reading endo ptr mesh " << std::endl;
	subregions_mesh_ptr[6].reset(new libMesh::Mesh(init.comm()));
	subregions_mesh_ptr[6] -> read("./bin-1.7.0/epi.e");*/

	// IDs for each region:
	// endo=2, epi=3, MV=4
	// LIPV=5, LSPV=6, RSPV=7, RIPV=8
//	std::string subregions_ids = "4";
	//std::string subregions_ids = "4, 5, 6, 7, 8, 2,3";
/*	std::vector < unsigned int > subregions_ids_vec;
	BeatIt::readList(subregions_ids, subregions_ids_vec);
	int default_id = 9;
*/
	// subregions
/*	int num_bc_subregions = subregions_mesh_ptr.size();
	int num_ids = subregions_ids_vec.size();

	std::vector < std::unique_ptr<libMesh::PointLocatorBase> > subregions_locator_ptr(num_bc_subregions);
	for (int i = 0; i < num_bc_subregions; ++i)
	{
		subregions_mesh_ptr[i]->print_info();

		std::cout<< " reset point locator tree" <<std::endl;
		subregions_locator_ptr[i].reset(new libMesh::PointLocatorNanoflann(*subregions_mesh_ptr[i]));
		//subregions_locator_ptr[i]->set_close_to_point_tol(bc_tolerance);
		subregions_locator_ptr[i]->enable_out_of_mesh_mode();
		subregions_locator_ptr[i]->_verbose=true;
	}

	std::cout << "Looping over nodes of the original mesh: " << std::endl;
	// set variables;
	int bc_missing_nodes = 0;
	std::set<int> bc_missing_node;
*/

//	for(auto & node : mesh.local_node_ptr_range() ){
/*		bool found=false;
		Point node_coords((* node)(0), (* node)(1), (* node)(2));
		nodes_mesh << "x = " << (* node)(0) << ", y = " << (* node)(1) << ", z = " << (* node)(0) << " \n";

		for(int i=0; i< num_bc_subregions; ++i){


			const Node* subregion_node = (subregions_locator_ptr[i]->locate_node(node_coords, nullptr ,bc_tolerance));
			if (subregion_node)
			{

				unsigned int node_id = i+1;
				if(num_ids == num_bc_subregions) node_id = subregions_ids_vec[i];
				mesh.boundary_info->add_node(node, node_id);
				found = true;
				nodes_found << "x = " << (* node)(0) << ", y = " << (* node)(1) << ", z = " << (* node)(0) << " \n";
				nodes_subregion << "x = " << (* subregion_node)(0) << ", y = " << (* subregion_node)(1) << ", z = " << (* subregion_node)(0) << " \n";
				break;
			}

			if (!found)
			{
				mesh.boundary_info->add_node(node, default_id);
				int nodeID = node->id();
				bc_missing_nodes++;
				bc_missing_node.insert(nodeID);
			}
*///		}

//		}


//}
//	if(false){
	// loop over elements, check 4 nodes, if same id, then assign set_id

/*	libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	el = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
	std::cout << "Looping over elements of the original mesh: " << std::endl;

	int c = 0;
	int id_endo = 2;
	int id_epi = 3;
	for (; el != end_el; ++el)
	{
		Elem * elem = *el;
		for (int s = 0; s < elem->n_sides(); ++s)
		{
					// check if element is on the surface of the mesh
					if (nullptr == elem->neighbor_ptr(s))
					{
						//std::cout << "Getting surface element" << std::endl;
						// surface element
						std::unique_ptr<libMesh::Elem> surf_elem = elem->side_ptr(s);

						// number of nodes of the surface element
						unsigned int n_nodes=surf_elem->n_nodes();
						//std::cout<< "n_nodes=" << n_nodes <<std::endl;

						//vector of ids associated with each node - here we are assuming each node has only one nodeset ID
						std::vector<short int> node_id_vec(n_nodes);
						//
	                    std::vector<short int> boundary_ids_vec(n_nodes);

	                    // this is specific for triangles

	                    //std::cout << "node 0" <<std::endl;
	                    const unsigned int i=0; // this is an important step
	                    Node * node_0=surf_elem->node_ptr(i);
	                    mesh.boundary_info->boundary_ids(node_0, boundary_ids_vec);
	                    node_id_vec[0]=boundary_ids_vec[0];

	                    //std::cout << "node 1" <<std::endl;
	                    const unsigned int j=1;
	                    Node * node_1=surf_elem->node_ptr(j);
						mesh.boundary_info->boundary_ids(node_1, boundary_ids_vec);
						node_id_vec[1]=boundary_ids_vec[0];

						//std::cout << "node 2" <<std::endl;
						const unsigned int k=2;
	                    Node * node_2=surf_elem->node_ptr(k);
						mesh.boundary_info->boundary_ids(node_2, boundary_ids_vec);
						node_id_vec[2]=boundary_ids_vec[0];


	                    // if all nodes have the same ID
	                    if(node_id_vec[0] > 3 && node_id_vec[0] < 9 ) mesh.boundary_info->add_side (elem, s, node_id_vec[0]);
	                    else if (node_id_vec[1] > 3 && node_id_vec[1]< 9) mesh.boundary_info->add_side (elem, s, node_id_vec[1]);
	                    else if (node_id_vec[2] > 3 && node_id_vec[2]<9) mesh.boundary_info->add_side (elem, s, node_id_vec[2]);
	                    else if( node_id_vec[0] ==2 || node_id_vec[1]==2 || node_id_vec[2]==2) mesh.boundary_info->add_side (elem, s, id_endo);
	                    else if( node_id_vec[0] ==3 && node_id_vec[1]==3 && node_id_vec[2]==3) mesh.boundary_info->add_side (elem, s, id_epi);
	                    //else mesh.boundary_info->add_side (elem, s, default_id);
*//*	                 }
		}
	}
}*/

/*	mesh.prepare_for_use();

	std::cout << "Export to new mesh" << std::endl;
	std::string output_file = data("output", "03meshLA_mv_out.e");
//	std::string output_file = data("output", "03meshLA_sidesets.e");
	ExodusII_IO exporter(mesh);
	exporter.write(output_file);

	//std::cout << "Missed " << bc_missing_node.size() << " boundary nodes." << std::endl;

	std::cout << "Good luck!" << std::endl;
	nodes_found.close();
	nodes_subregion.close();
	nodes_mesh.close();

    return 0;
*/
}


