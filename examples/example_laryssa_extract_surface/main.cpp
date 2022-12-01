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
#include "libmesh/node.h"
#include <unordered_set>
#include <libmesh/mesh_communication.h>
struct intersection_of_blockIDs
{
    libMesh::dof_id_type elID;
    unsigned short side;
    libMesh::dof_id_type n1ID, n2ID;
};


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
    libMesh::ReplicatedMesh mesh(init.comm());
    libMesh::Mesh mesh3(init.comm());
    libMesh::ExodusII_IO importer(mesh);
    std::string name_mesh_sync = data("mesh_sync", "NONE");   
    bool assign_blockid_intersection = data("assign_blockid_intersection", 0);   
    bool create_sync_mesh=data("create_sync_mesh", 0);    
    std::cout << "create_sync_mesh = " << create_sync_mesh << std::endl;

// Step 1 - extract surface     
if(create_sync_mesh){
        std::cout << "reading mesh " << name_mesh_sync << std::endl;
        mesh.read(name_mesh_sync);
        std::cout << "extracting surface automatically using get_boundary_info().sync" << std::endl;
        mesh.get_boundary_info().sync(mesh3);
        std::cout << "getting boundary info for mesh3" << std::endl;
        mesh3.print_info();
        std::string output_path3="RESILIENT_LA_surface3.e";
        ExodusII_IO (mesh3).write(output_path3); // creates the exporter
        return 0;
}    

    std::string name_mesh = data("mesh", "NONE");   
    std::cout << "Reading mesh " << name_mesh << std::endl;  
    importer.read(name_mesh);
    mesh.prepare_for_use();
    mesh.print_info();

    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    double tol = 1e-6;

// Step 2 - assign block id 8 to intersection
if(assign_blockid_intersection){
    std::cout << "ASSIGNING BLOCK ID 8 TO INTERSECTION - THIS MAY TAKE A WHILE" << std::endl;
    int elnum=0;
    // set block to 1 and count number of elements on the boundary

    int i1=0;
    int i2=1;
    int i3=2;
    int i4=3;
    el = mesh.active_local_elements_begin();
    for (; el != end_el; ++el)
    {
        if(elnum % 100000==0) std::cout << "element number " << elnum << std::endl;
        Elem * elem = * el;
        Node * elnode1 = elem->Elem::node_ptr(i1);
        Node * elnode2 = elem->Elem::node_ptr(i2);
        Node * elnode3 = elem->Elem::node_ptr(i3);
        if ( elem->subdomain_id()==10)
        {
                 auto centroid1 = elem->true_centroid();
                 for(libMesh::MeshBase::const_element_iterator it=mesh.active_local_elements_begin(); it != end_el; it++)
                 {
                        if ((*it)->subdomain_id() == 92)
                        {

                            auto centroid2 = (*it)->true_centroid();
                            auto centroid_dif = centroid1-centroid2;
                            if(centroid_dif.norm() <=  tol)
                            {
                                elem->subdomain_id()=8;
                                (*it)->subdomain_id()=8;
                            }
                        }
                }
        }
        elnum++;
     }
 }

/*
// Step 3 - go through elements with blockid8 and create struct
 
    // reinitialize el if need to use it again
    el = mesh.active_local_elements_begin();
    std::vector<intersection_of_blockIDs> endocardium;
    std::vector<intersection_of_blockIDs> epicardium;    
    for (; el != end_el; ++el)
    {
        Elem * elem = *el;
        if(elem->subdomain_id()==8)   
        {
            std::cout << "---------- ELEM " << elem->id() << " -----------------" <<  std::endl;

            for (auto side : elem->side_index_range())
            {
                Elem * neighboor_elem = elem->neighbor_ptr(side);

                // ENDOCARDIUM
                if (neighboor_elem->subdomain_id()==1)  
                {
                    intersection_of_blockIDs intersection_endo;
                    bool found_first_node_endo=0;
                    const unsigned int ss = side;
                    for (unsigned int i=0; i<3; i++)
                    {
                        const unsigned int ii=i;
                        Node * nn = elem-> node_ptr(ii);
                        for (unsigned int j=0; j<3; j++)
                        {
                            const unsigned int jj=j;
                            Node * mm = neighboor_elem-> node_ptr(jj);
                            if(mm == nn){
                        //        std::cout << "Elements " << elem->id() << " and " << neighboor_elem->id() << " share node " << nn->id() << std::endl; 
                                intersection_endo.elID = neighboor_elem->id();
                                intersection_endo.side = side;
                                if(!found_first_node_endo){
                                    intersection_endo.n1ID = nn->id(); 
                                    found_first_node_endo=1;
                                }
                                else{
                                    intersection_endo.n2ID = nn->id();
                                    endocardium.push_back(intersection_endo);
                                }
                            }   
                        }
                    }   
                }
                // EPICARDIUM
                else if (neighboor_elem->subdomain_id()==2)  
                {
                    intersection_of_blockIDs intersection_epi;
                    bool found_first_node_epi=0;
                    const unsigned int ss = side;
                    for (unsigned int i=0; i<3; i++)
                    {
                        const unsigned int ii=i;
                        Node * nn = elem-> node_ptr(ii);
                        for (unsigned int j=0; j<3; j++)
                        {
                            const unsigned int jj=j;
                            Node * mm = neighboor_elem-> node_ptr(jj);
                            if(mm == nn){
                                intersection_epi.elID = neighboor_elem->id();
                                intersection_epi.side = side;
                                if(!found_first_node_epi){
                                    intersection_epi.n1ID = nn->id(); 
                                    found_first_node_epi=1;
                                }
                                else{
                                    intersection_epi.n2ID = nn->id();
                                    epicardium.push_back(intersection_epi);
                                }
                            }   
                        }
                    }   
                }
            }
        }
    }

// ASSIGN NODE IDS TO THE NODES IN THE INTERSECTIONS
    std::cout << "Assigning nodesets" << std::endl;
    int nodeset_id_epi=10;
    int nodeset_id_endo=9;
    std::cout << " -----------ENDOCARDIUM-------------" << std::endl;
    for( auto endo:endocardium)
    {
        std::cout << "Element " << endo.elID << " shares nodes " << endo.n1ID << " and " << endo.n2ID << " with the epicardium on side" << endo.side << std::endl;       
        mesh.boundary_info->add_node(endo.n1ID, nodeset_id_endo);
        mesh.boundary_info->add_node(endo.n2ID, nodeset_id_endo);
    }
    std::cout << " -----------EPICARDIUM-------------" << std::endl;
    for(auto epi:epicardium)
    {
        mesh.boundary_info->add_node(epi.n1ID, nodeset_id_epi);
        mesh.boundary_info->add_node(epi.n2ID, nodeset_id_epi);
        std::cout << "Element " << epi.elID << " shares nodes " << epi.n1ID << " and " << epi.n2ID << " with the endocardium on side" << epi.side << std::endl;       
    }
*/
    std::cout << "deleting elements in the intersection" << std::endl;
    std::unordered_set<int> nodes_to_delete;
    el = mesh.active_local_elements_begin();
    for (; el != end_el; ++el)
    {
        Elem * elem = *el;
        if(elem->subdomain_id()==8)
        {
         // delete nodes -- automatically happens when element id deleted   
          /*  for (unsigned int i=0; i<3; i++)
            {

                const unsigned int ii=i;
                const Node * nn = elem-> node_ptr(ii);
                
                std::vector<boundary_id_type> boundary_ids;
                mesh.boundary_info->boundary_ids(nn, boundary_ids); 
                short int i0=0;
                if(boundary_ids.empty())
                {
                    const unsigned int nn_id = nn->id();
                    std::cout << "deleting node" << std::endl;
                    std:: cout << "local node node id = " << i << std::endl;
                    std:: cout << "global node node id = " << nn_id << std::endl;
                    // delete note if have not deleted yet
                    if(nodes_to_delete.find(nn_id)!=nodes_to_delete.end()){
                        Node * nodedel = mesh.node_ptr(nn_id);
                        std::cout << *nodedel << std::endl;
                        nodes_to_delete.insert(nn_id);
                        mesh.delete_node(nodedel);
                    }
                }
            }*/
            mesh.delete_elem(elem);
        }
        else
        {
            elem->subdomain_id()=1;
        }
    }

/*
    libMesh::MeshBase::const_node_iterator nd1 = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator end_nd = mesh.active_nodes_end();
    
    MeshCommunication::assign_global_indices(& mesh);
    Elem * elem1 = mesh.elem_ptr(9);
    Node * node1 = mesh.node_ptr(0);
    Node * node2 = mesh.node_ptr(4);
    libMesh::dof_id_type n1_id = node1->id();
    std::cout << n1_id << std::endl; 
    unsigned int n1_id_local = elem1->local_node(n1_id);*/ 
    //*(elem1->node_ptr(n1_id_local)) =* node2;
    //mesh.delete_node(node1);                           
    /*for (; nd1 != end_nd; ++nd1)
    {
        Node * node1 = *nd1;
        std::vector<boundary_id_type> boundary_ids1;
        mesh.boundary_info->boundary_ids(node1, boundary_ids1); 
        if(!boundary_ids1.empty())
        {   
            if(boundary_ids1[0]==9)
            {
                std::cout << "node " << node1->id() << " has boundary id 9" << std::endl;
                for (libMesh::MeshBase::const_node_iterator nd2=mesh.active_nodes_begin(); nd2 != end_nd; nd2++)
                {
                    Node * node2 = *nd2;
                    std::vector<boundary_id_type> boundary_ids2;
                    mesh.boundary_info->boundary_ids(node2, boundary_ids2); 
                    if(!boundary_ids2.empty())
                    {   
                        if(boundary_ids2[0]==10)
                        {
                            std::cout << " node2 " << node2->id() << " has boundary id 10" << std::endl;
                            Point p((*node1)(0)-(*node2)(0), (*node1)(1)-(*node2)(1),  (*node1)(2)-(*node2)(2));
                            if(p.norm() < tol && node1->id() != node2->id()){
                                std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;
                                std::cout << "nodes " << node1->id() << " and " << node2->id() << " are the same" << std::endl;
                                std:: cout << ".........node " << node1->id() << std::endl;
                                std::cout << *node1 << std::endl;
                                std:: cout << ".........node " << node2->id() << std::endl;
                                std::cout << *node2 << std::endl;
                                mesh.delete_node(node1);                           
                                Elem * elem1 = mesh.elem_ptr()
                                elem1->local_node(node1->id()) = node2;
                                std::cout << "---------------------------------------------- after node1=node2" << std::endl;
                                std::cout << "nodes " << node1->id() << " and " << node2->id() << " are the same" << std::endl;
                                std:: cout << ".........node " << node1->id() << std::endl;
                                std::cout << *node1 << std::endl;
                                std:: cout << ".........node " << node2->id() << std::endl;
                                std::cout << *node2 << std::endl;
                            }
                        }
                    }
                 }
             }       
          }
    }
    nd1 = mesh.active_nodes_begin();
    for (; nd1 != end_nd; ++nd1)
    {
        Node * np = *nd1;
        std::cout << " node = " << np->id() << std::endl;
     //   if(np->id()==3) std::cout << *np << std::endl;
    }
    
    el = mesh.active_local_elements_begin();
    for (; el != end_el; ++el)
    {
        Elem * elem = *el;
        std::cout << "element " << elem->id() << std::endl ;
        elem->print_info();
    }*/
    std::cout << "preparing mesh for use" << std::endl;
    mesh.prepare_for_use();
    //mesh.renumber_nodes_and_elements();    

    /*std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" <<  std::endl;
    nd1 = mesh.active_nodes_begin();
    for (; nd1 != end_nd; ++nd1){
        Node * nn =*nd1;
        std:: cout << ".......node " << nn->id() << std::endl;
        std::cout << (*nn) << std::endl;
    }*/
    



    std::string output_path="RESILIENT_LA_SURFACE.e";


    // Export fibers in nemesis format
    std::cout << "exporting mesh" << std::endl;
    ExodusII_IO nemesis_exporter(mesh); // creates the exporter
    std::cout << "Export mesh" << std::endl;
    nemesis_exporter.write(output_path);   
    return 0;

}
