/*
 * main_laryssa_define_subregions
 * based and inspired in main_remap_regions.cpp by srossi
 *
 *  Created on: Oct 25, 2021
 *      Author: labdala
 */
#include<iostream>
#include<fstream>
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
#include "libmesh/enum_elem_quality.h"
#include "libmesh/elem_quality.h" 
#include "libmesh/dense_matrix.h"
#include <iomanip>
#include <Eigen/Dense>

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
    libMesh::ExodusII_IO importer(mesh);
    std::string name_mesh=data("mesh", " ");

    std::cout << "Reading mesh" << std::endl;  

    mesh.read(name_mesh);
    mesh.print_info();
    libMesh::MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    libMesh::MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    const ElemQuality jac= JACOBIAN; 
    std::string description= libMesh::Quality::describe(jac);
    std::cout << description << std::endl;
    std::ofstream quality_file;
    std::ofstream contravariant_file;
    quality_file.open("quality.txt");
    contravariant_file.open("contravariant.txt");
    long double determinant;
    quality_file  <<std::fixed  <<  std::setprecision(30) << endl;
    contravariant_file  << std::setprecision(30) << endl;
    for(; el!=end_el; el++){
        libMesh::Elem * elem = * el; 
        /*elem->print_info();
        break;*/
        const libMesh::Point p0 = elem->point(0);
        const libMesh::Point p1 = elem->point(1);
        const libMesh::Point p2 = elem->point(2);
        const libMesh::Point p3 = elem->point(3);

        Eigen::Matrix<long double, 3, 3> contravariant;
        //libMesh::DenseMatrix<double> contravariant(3,3);
        contravariant(0, 0) = p1(0) - p0(0);
        contravariant(0, 1) = p2(0) - p0(0);
        contravariant(0, 2) = p3(0) - p0(0);
        contravariant(1, 0) = p1(1) - p0(1);
        contravariant(1, 1) = p2(1) - p0(1);
        contravariant(1, 2) = p3(1) - p0(1);
        contravariant(2, 0) = p1(2) - p0(2);
        contravariant(2, 1) = p2(2) - p0(2);
        contravariant(2, 2) = p3(2) - p0(2);
/*        contravariant_file << contravariant(0, 0) << " " << contravariant(0, 1) << " " << contravariant(0, 2) << " "
        << contravariant(1, 0) << " " << contravariant(1, 1) << " " << contravariant(1, 2) << " " 
        << contravariant(2, 0) << " " << contravariant(2, 1) << " " << contravariant(2, 2) << "\n ";
*/
        if(contravariant(2,0)!=0 && contravariant(2,1)!=0 && contravariant(2,2)!=0)
        {
         contravariant_file << contravariant(2, 0) << " " << contravariant(2, 1) << " " << contravariant(2, 2) << "\n ";
        }
        long double x = ((contravariant(1,1) * contravariant(2,2)) - (contravariant(2,1) * contravariant(1,2)));
        long double y = ((contravariant(1,0) * contravariant(2,2)) - (contravariant(2,0) * contravariant(1,2)));
        long double z = ((contravariant(1,0) * contravariant(2,1)) - (contravariant(2,0) * contravariant(1,1)));
    
        determinant = contravariant.determinant();// ((contravariant(0,0) * x) - (contravariant(0,1) * y) + (contravariant(0,2) * z));
      // if(determinant < 0){
       //     std::cout << "hi" << std::endl;
            quality_file << "elem_id = " << elem->id() << ", jac = " << determinant << ", elem_volume = " << elem->volume() << std::endl; 
      // }
    }

    quality_file.close();
    contravariant_file.close();
    return 0;

}
