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
    libMesh::ExodusII_IO importer(mesh);


    std::cout << "Reading mesh" << std::endl;  


   // mesh.read("./F72_fsi_fibers.e");
    importer.read("./fastl_fibers.e");
    mesh.prepare_for_use();
    mesh.print_info();
   /* mesh.subdomain_name(1) = "left_atrium";
    mesh.get_boundary_info().sideset_name(2)="endo";
    mesh.get_boundary_info().sideset_name(3)="epi";
    mesh.get_boundary_info().sideset_name(4)="mvr";
    mesh.get_boundary_info().sideset_name(5)="pv3";
    mesh.get_boundary_info().sideset_name(6)="pv4";
    mesh.get_boundary_info().sideset_name(7)="pv2";
    mesh.get_boundary_info().sideset_name(8)="pv1a";
    mesh.get_boundary_info().sideset_name(9)="pv1b";
 */   libMesh::EquationSystems  equation_systems(mesh);    
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& fiber_system = equation_systems.add_system<FiberSystem>("fibers");
    fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.init();
    /*auto& f_v = fiber_system.solution;
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
*/
    equation_systems.print_info();

    unsigned int step=1;
    auto& elem_var_names = importer.get_elem_var_names();
    for (auto && var : elem_var_names) std::cout << var << std::endl;

    const int num_fiber_systems =1;// 3;
    std::vector < std::string > fibers(num_fiber_systems);
    fibers[0] = "fibers";
   /* fibers[1] = "sheets";
    fibers[2] = "xfibers";
*/

    for (auto && var : elem_var_names) std::cout << var << std::endl;
    auto elem_first = elem_var_names.begin();
    auto elem_end = elem_var_names.end();
    for (int k = 0; k < num_fiber_systems; ++k)
    {
       libMesh::System& system = equation_systems.get_system(fibers[k]);
       std::string name = system.name();
       std::cout << "Importing System: " << name << std::endl;
       int n_vars = system.n_vars();
       std::cout << "nvars = " << n_vars << std::endl;
       for (int l = 0; l < n_vars; ++l)
       {
           std::cout << "l=" << l << std::endl;
           std::string var_name = system.variable_name(l);
           auto elem_it = std::find(elem_first, elem_end, var_name);
           std::cout << "elem_it=" << * elem_it << std::endl;
           if (elem_it != elem_end)
           {
               std::cout << "\t elemental variable: " << *elem_it << std::endl;
               importer.copy_elemental_solution(system, *elem_it, *elem_it, step);
           }
       }
     }

    std::string output_path="fastl_fibers_lessvariables.e";
    //VTKIO(mesh).write_equation_systems("F72_fsi_fibers_named.pvtu", equation_systems);

        // Export fibers in nemesis format
    ExodusII_IO nemesis_exporter(mesh); // creates the exporter
    std::vector<std::string> output_variables(3); // creates a vector that stores the names of vectors
    output_variables[0] = "fibersx";
    output_variables[1] = "fibersy";
    output_variables[2] = "fibersz";
   /* output_variables[3] = "sheetsx";
    output_variables[4] = "sheetsy";
    output_variables[5] = "sheetsz";
    output_variables[6] = "xfibersx";
    output_variables[7] = "xfibersy";
    output_variables[8] = "xfibersz";
*/
    std::cout << "Export mesh" << std::endl;
    nemesis_exporter.write(output_path);
    nemesis_exporter.set_output_variables(output_variables);
    nemesis_exporter.write_equation_systems(output_path, equation_systems); 
    nemesis_exporter.write_element_data(equation_systems);

    return 0;

}
