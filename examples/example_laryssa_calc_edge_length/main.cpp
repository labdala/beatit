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
#include "libmesh/mesh_refinement.h"
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

// Initialize libmesh
using namespace libMesh;


double dist(Point v1, Point v2)
{
    double d=0;
    d = (v1(0)-v2(0))*(v1(0)-v2(0)) + (v1(1)-v2(1))*(v1(1)-v2(1)) + (v1(2)-v2(2))*(v1(2)-v2(2));
    d = std::sqrt(d);
    return d;
}

int main(int argc, char ** argv)
{


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

    mesh.read("./f72_cubit_test.e");
    mesh.print_info();


	libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
	const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

	el = mesh.active_local_elements_begin();
	std::cout << "Looping over elements of the original mesh: " << std::endl;

    int nodes =0;
	double average_edge_length = 0;
	double average_edge_length_local = 0;


    int i1=0;
    int i2=1;
    int i3=2;
    int i4=3;
    double e_min=1.0;
    double e_max=0.0;
    std::vector<double> e(6);
    int nelem=0;
    for (; el != end_el; ++el)
    {
        average_edge_length_local=0;
        Elem * elem = *el;

        Node * node1 = elem->node_ptr(i1);
        Node * node2 = elem->node_ptr(i2);
        Node * node3 = elem->node_ptr(i3);
        Node * node4 = elem->node_ptr(i4);

        Point x1((* node1)(0), (* node1)(1), (* node1)(2));
        Point x2((* node2)(0), (* node2)(1), (* node2)(2));
        Point x3((* node3)(0), (* node3)(1), (* node3)(2));
        Point x4((* node4)(0), (* node4)(1), (* node4)(2));
        
        e[0] = dist(x1,x2);
        e[1] = dist(x1,x3);
        e[2] = dist(x1,x4);
        e[3] = dist(x2,x3);
        e[4] = dist(x2,x4);
        e[5] = dist(x3,x4);
        for(int i=0; i<6; i++)
        {
           // std::cout << "e[i]=" << e[i] << std::endl;
            average_edge_length_local += e[i];
            if(e[i]<e_min) e_min=e[i];
            if(e[i]>e_max) e_max=e[i];
        }
        average_edge_length += average_edge_length_local/6;
        nelem++; 
    }
    average_edge_length = average_edge_length/nelem;

    std::cout << "average_edge_length = " << average_edge_length << std::endl;
    std::cout << "e_min = " << e_min << std::endl;
    std::cout << "e_max = " << e_max << std::endl;
    

	std::cout << "Good luck!" << std::endl;
    return 0;

}


