/*
 * brute_force_assign_material.cpp
 *
 *  Created on: Oct 20, 2022
 *      Author: labdala
 */


/* Main mesh should be a volume mesh
 * the "boundary condition" meshes should be surface meshes
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
#include <libmesh/equation_systems.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/dof_map.h>

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

    std::string mesh_file = data("input_mesh", "NONE");
    std::cout << "-------------------" <<std::endl;
    if ("NONE" != mesh_file)
    {
        std::cout << "Reading mesh " << mesh_file << std::endl;
        mesh.read(mesh_file);
    }
    else
    {
        std::cout << "Set input_mesh in input file!" << std::endl;
        throw std::runtime_error("No mesh given");
    }
    mesh.print_info();

    double tolerance = data("tolerance", 1e-3);
    int default_id = data("default_id", 666);

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    libMesh::EquationSystems es(mesh);
    LinearImplicitSystem& material_sys = es.add_system<LinearImplicitSystem>("material_system");
    material_sys.add_variable("Material", libMesh::FIRST, libMesh::LAGRANGE);
    material_sys.init();
    auto& mm = material_sys.solution;

    std::vector<dof_id_type> material_dof_indices; //putting this back to the system so that we can export it


// WE WANT MATERIAL TO BE CELL CENTERED. WE DO NOT NEED THE NODES OTHER THAN FOR CHECKING IF ELEMENT IS THE SAME
    // dof_indices.print();

    // set bc flags
    //bool set_bc_flags = data("bc", false);
    if(true)
    {
        std::string bc_meshes = data("bcs", "NONE");
        std::string material_ids = data("material", "");
        std::vector < std::string > bc_files_vec;
        BeatIt::readList(bc_meshes, bc_files_vec);
        std::cout << "bc meshes" << std::endl;
        for(auto && m : bc_files_vec) std::cout << m << std::endl;
        std::vector < unsigned int > material_ids_vec;
        BeatIt::readList(material_ids, material_ids_vec);
        std::cout << "material values for each bc mesh" << std::endl;
        for(auto && m : material_ids_vec) std::cout << m << std::endl;
        int num_bc_subregions = bc_files_vec.size();
        int num_ids = material_ids_vec.size();
        std::cout << "Num bc meshes: " << num_bc_subregions << ", num ids: " << num_ids << std::endl;
        double bc_tolerance = data("bc_tolerance", 1e-2);

        std::cout << "Reading "<< num_bc_subregions <<" bc meshes and creating point locator ...  " << std::endl;
        std::vector < std::unique_ptr<libMesh::Mesh> > bc_region_meshes_ptr(num_bc_subregions);

        for (int i = 0; i < num_bc_subregions; ++i)
        {
            std::cout << "-------------------" << std::endl;
            std::cout << "Reading bc mesh " << bc_files_vec[i] << ", " << i << std::endl;
            bc_region_meshes_ptr[i].reset(new libMesh::Mesh(init.comm()));
            bc_region_meshes_ptr[i]->read(bc_files_vec[i]);
            bc_region_meshes_ptr[i]->print_info();
        }

        int bc_missing_elements = 0;
        std::set<int> bc_missing_elem;

        std::cout << "Looping over elements: " << std::endl;

        int c = 0;
        for (; el != end_el; ++el)
        {
            c++;
            Elem * elem = *el;
            if( c % 20000 == 0){
                std::cout << "-----------Elem: " << c << std::endl;
                //elem->print_info();
            } 

            material_sys.get_dof_map().dof_indices(elem, material_dof_indices);

            //loop on the sides of the element
            for (int s = 0; s < elem->n_sides(); ++s)
            {
                // check if the element is on the boundary
               if (nullptr == elem->neighbor_ptr(s))
               {
                    auto side_elem_ptr  = elem->build_side_ptr(s);
                    auto node0 = side_elem_ptr->node_ptr(0);
                    auto node1 = side_elem_ptr->node_ptr(1);
                    auto node2 = side_elem_ptr->node_ptr(2);

                    bool found_elem = false;
                    int region = -1;
                    // Let's loop over the number of subregion meshes and look for the node
                    // As soon as we find it, we can pass to the next element
                    for (int i = 0; i < num_bc_subregions; ++i)
                    {
                        bool found_node = false;
                        libMesh::MeshBase::const_node_iterator node =  bc_region_meshes_ptr[i]->active_nodes_begin();
                        const libMesh::MeshBase::const_node_iterator end_nodes =  bc_region_meshes_ptr[i]->active_nodes_end();
                        //looping over bc mesh and checking if node0 is in any of them3e
                        
                        //std::cout << "looping over nodes of " << bc_meshes[i] << "mesh" << std::endl;
                        for (; node != end_nodes; ++node)
                        {
                            libMesh::Node * bc_node = *node;
                            libMesh::Point v(*bc_node - *node0);
                            double distance = v.norm();
                            if(distance < tolerance)
                            {
                                found_node = true;
                                region = i;
                                break;
                            }
                        }
                        // // if the node is in at least one of them, go check what is the specific node
                        if(found_node)
                        {

                            libMesh::MeshBase::const_element_iterator bc_el = bc_region_meshes_ptr[i]->active_local_elements_begin();
                            const libMesh::MeshBase::const_element_iterator end_bc_el = bc_region_meshes_ptr[i]->active_local_elements_end();
                            //loop over bc mesh elements
                            for (; bc_el != end_bc_el; ++bc_el)
                            {
                                Elem * bc_elem = *bc_el;
                                auto bc_node0 = bc_elem->node_ptr(0);
                                auto bc_node1 = bc_elem->node_ptr(1);
                                auto bc_node2 = bc_elem->node_ptr(2);

                                // check if node0 is any of the  vertices of the triangle
                                // here assuming that the bc mesh is a surface one and has TRI elements
                                libMesh::Point v0(*bc_node0 - *node0);
                                double distance0 = v0.norm();
                                libMesh::Point v1(*bc_node1 - *node0);
                                double distance1 = v1.norm();
                                libMesh::Point v2(*bc_node2 - *node0);
                                double distance2 = v2.norm();
                                //std::cout << "A - distance0 = " << distance0<< ", distance1 = " << distance1<< ", distance2 = " << distance2 << std::endl;
                                //check if any of the distances is small
                                if(distance0 < tolerance || distance1 < tolerance || distance2 < tolerance)
                                {
                                    v0  = *bc_node0 - *node1;
                                    distance0 = v0.norm();
                                    v1 = *bc_node1 - *node1;
                                    distance1 = v1.norm();
                                    v2 = *bc_node2 - *node1;
                                    distance2 = v2.norm();
                                    //std::cout << "*node1" << *node1 <<std::endl;
                                    //std::cout << "A - distance0 = " << distance0<< ", distance1 = " << distance1<< ", distance2 = " << distance2 << std::endl;

                                    // narrow down. if an element of the bc mesh shares three nodes with the 
                                    // main mesh, then they must be the same element
                                    if(distance0 < tolerance || distance1 < tolerance || distance2 < tolerance)
                                    {
                                         v0  = *bc_node0 - *node2;
                                         distance0 = v0.norm();
                                         v1 = *bc_node1 - *node2;
                                         distance1 = v1.norm();
                                         v2 = *bc_node2 - *node2;
                                         distance2 = v2.norm();
                                        //  std::string tf= (node2==nullptr)? "true" : "false";
                                        // std::cout << "tf" << tf <<std::endl;
                                        // std::cout << "v0 = " << v0 << std::endl;// distance0<< ", distance1 = " << distance1<< ", distance2 = " << distance2 << std::endl;
                                         // check last node in elements 
                                         if(distance0 < tolerance || distance1 < tolerance || distance2 < tolerance)
                                         {
                                           found_elem = true;
                                           unsigned int material_value = i+1;
                                           if(num_ids == num_bc_subregions) material_value = material_ids_vec[i];
                                            for (int idim = 0; idim < 3; ++idim)
                                            {
                                                //mesh.boundary_info->add_side (elem, s, sideset_id);
                                                mm->set(material_dof_indices[idim], material_value);
                                            } 
                                           break;
                                         }
                                     }
                                 } // end distance check
                             } // end loop over elem bcs
                        } // end if found node
                        if(found_elem) break;
                    } // end loop over mesh subregions

                    if (!found_elem)
                    {
                        bc_missing_elements++;
                        //std::cout << "\n Missed elem: " << bc_missing_elements << std::endl;
                    }

                }
            }
        }
        std::cout << "Missed " << bc_missing_elem.size() << " boundary elements." << std::endl;

        mesh.prepare_for_use();
     }
    std::cout << "Export to new mesh" << std::endl;
    std::string output_file = data("output", "remapped_mesh.e");
    ExodusII_IO exporter(mesh);
 
    std::vector<std::string> output_variables(1);
    output_variables[0] = "Material";
    exporter.set_output_variables(output_variables);
    exporter.write_equation_systems(output_file, es); // saves the variables at the nodes
	exporter.write_element_data(es); // this saves the variables at the centroid of the elements


    std::cout << "Good luck!" << std::endl;
    return 0;
}

