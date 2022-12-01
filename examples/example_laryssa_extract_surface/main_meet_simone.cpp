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
#include "libmesh/equation_systems.h"
#include "libmesh/vtk_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/explicit_system.h"

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
    libMesh::Mesh mesh_surface(init.comm());
    libMesh::Mesh mesh3(init.comm());
    libMesh::ExodusII_IO importer(mesh);
    std::string name_mesh = data("mesh", "NONE");   
    bool create_sync_mesh=data("create_sync_mesh", FALSE);    
 
    std::cout << "setting dimensions" << std::endl;
    mesh_surface.set_mesh_dimension(2); //CHECK THIS        
    mesh_surface.set_spatial_dimension(3);

    std::cout << "Reading mesh" << std::endl;  
    importer.read(name_mesh);
    mesh.prepare_for_use();
    mesh.print_info();

    if(create_sync_mesh){
        std::cout << "extracting surface automatically using get_boundary_info().sync" << std::endl;
        mesh.get_boundary_info().sync(mesh3);
        std::cout << "getting boundary info for mesh3" << std::endl;
        std::string output_path3="RESILIENT_LA_surface3.e";
        ExodusII_IO (mesh3).write(output_path3); // creates the exporter
    }    

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    // set block to 1 and count number of elements on the boundary
    double tol = 1e-6;
    int nelem_bd=0;
    el = mesh.active_local_elements_begin();
    for (; el != end_el; ++el)
    {
        Elem * elem = * el;
        // Get elem on the surface of the mesh
        for (auto side : elem->side_index_range())
        {
           // elem->Elem::subdomain_id()=1;  
            if (elem->neighbor_ptr(side) == nullptr)
            {
                 nelem_bd++;
                 auto side_el = elem->build_side_ptr(side);
                 Node * n0 = side_el->node_ptr(0); 
                 for(libMesh::MeshBase::const_element_iterator it=mesh.active_local_elements_begin(); it != end_el; it++)
                 {
                    if(elem->subdomain_id() != (*it)->subdomain_id()){
                        for(auto side2: (*it)->side_index_range()){
                        if ((*it)->neighbor_ptr(side2) == nullptr)
                        {
                            auto side2_el = (*it)->build_side_ptr(side2);
                            Node * n20 = side2_el->node_ptr(0);
                            Node * n21 = side2_el->node_ptr(1);
                            Node * n22 = side2_el->node_ptr(2);
                            *n20 -= *n0;
                            *n21 -= *n0;
                            *n22 -= *n0;
                            if(n20->norm() <=  tol || n21->norm() <=  tol ||n22->norm() <=  tol   ){
                                std::cout << "found a node, n0id = " << n0->id() << ", n20 = " << n20->id() <<", n21 = " << n21->id() << ", n22 = " << n22->id()  << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}
    //maybe dont need below
    mesh_surface.reserve_elem(nelem_bd);
    unsigned int node_id = 0;
    unsigned int elem_id = 0;
    el = mesh.active_local_elements_begin();
    // assign mesh nodes
    std::cout << "looping over elements" << std::endl;
    for (; el != end_el; ++el)
    {
        Elem * elem = * el;
        // Get elem on the surface of the mesh
        for (auto side : elem->side_index_range())
        {
            if (elem->neighbor_ptr(side) == nullptr)
            {
                        auto side_elem = elem->build_side_ptr(side);
                        Node * n0 = side_elem->node_ptr(0);
                        Node * n1 = side_elem->node_ptr(1);
                        Node * n2 = side_elem->node_ptr(2);
                        mesh_surface.add_node(n0);
                        mesh_surface.add_node(n1);
                        mesh_surface.add_node(n2);
                        mesh_surface.add_elem(side_elem.get());
            }
        }
    }
    std::cout << "preparing mesh for use" << std::endl;
    mesh_surface.prepare_for_use();
    std::string output_path="RESILIENT_LA_surface.e";


    // Export fibers in nemesis format
    std::cout << "exporting mesh" << std::endl;
    ExodusII_IO nemesis_exporter(mesh_surface); // creates the exporter
    std::cout << "Export mesh" << std::endl;
    nemesis_exporter.write(output_path);

    return 0;

}
