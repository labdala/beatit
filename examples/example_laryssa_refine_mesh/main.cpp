/*
 * main_laryssa_refine_mesh
 *
 *  Created on: Oct, 2022
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
#include <cctype>
#include <cstring>
#include <string.h>

int main(int argc, char ** argv)
{

     // Import input file
    GetPot data = BeatIt::readInputFile(argc, argv);

    libMesh::LibMeshInit init(argc, argv, MPI_COMM_WORLD);


    // Create mesh:
    // one for the conforming mesh
    // one with the regions to be mapped
    libMesh::Mesh mesh(init.comm());
    libMesh::ExodusII_IO importer(mesh);

    
    std::string mesh_name = data("mesh","myocardium.e");
    std::cout << "Reading mesh "<< mesh_name << std::endl;

    // mesh.read(mesh_name);
    // mesh.print_info();
    importer.read(mesh_name);
    mesh.prepare_for_use();
    mesh.print_info();

    //Check if exodus has fiberZ or fiberz
    std::string lower_or_upper;
    std::cout << "Variables available in the exodus file:" << std::endl;
    auto& elem_var_names = importer.get_elem_var_names();
    for (auto && var : elem_var_names){
        std::cout << var << std::endl;
        lower_or_upper=var.back();
    }

    std::string var_vec;
    if (lower_or_upper=="z")
    {
        // strcpy (var_vec, "xyz");
        char xvar='x';
        char yvar='y';
        char zvar='z';
        var_vec = var_vec+xvar+yvar+zvar;
    } else if (lower_or_upper=="Z") {
        char xvar='X';
        char yvar='Y';
        char zvar='Z';
        var_vec = var_vec+xvar+xvar+yvar+zvar;
    }
    else
    {
        throw std::runtime_error("This mesh does not satisfy the assumption this mesh is three dimensional and has fiber fields");
    }

    libMesh::EquationSystems  es_fiber(mesh);

    std::cout << "Adding fibers to Eq Syst." << std::endl;
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& fiber_system = es_fiber.add_system<FiberSystem>("fibers");
    for( int i=0; i<3; i++)
    {
        std::string fib("fibers");
        fib.push_back(var_vec[i]);
        fiber_system.add_variable(fib, libMesh::CONSTANT, libMesh::MONOMIAL);
    }
    fiber_system.init();
    auto& ff = fiber_system.solution;

    std::cout << "Adding xfibers to Eq Syst." << std::endl;
    typedef libMesh::ExplicitSystem  xFiberSystem;
    xFiberSystem& xfiber_system = es_fiber.add_system<FiberSystem>("xfibers");
    for( int i=0; i<3; i++)
    {
        std::string xfib("xfibers");
        xfib.push_back(var_vec[i]); 
        xfiber_system.add_variable(xfib, libMesh::CONSTANT, libMesh::MONOMIAL);
    }
    xfiber_system.init();
    auto& ss = xfiber_system.solution;

    std::cout << "Adding sheets to Eq Syst." << std::endl;
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& sheets_system = es_fiber.add_system<FiberSystem>("sheets");
    for( int i=0; i<3; i++)
    {
        std::string sht("sheets");
        sht.push_back(var_vec[i]); 
        sheets_system.add_variable(sht, libMesh::CONSTANT, libMesh::MONOMIAL);
    }
    sheets_system.init();
    auto& nn = sheets_system.solution;

    std::cout << " +++++++++++++ EquationSystems setup in the mesh where we have fiber info \n";
    es_fiber.print_info();

    //Read fiber data into mesh_fiber
    std::cout << "Copying the exodus variables from into Eq Syst." << std::endl;
    unsigned int step=1;
    //Read fiber data into mesh_fiber
    for( int i=0; i<3; i++)
    {
        std::string fib("fibers");
        fib.push_back(var_vec[i]);
        importer.copy_elemental_solution(fiber_system, fib, fib, step);
    }
     for( int i=0; i<3; i++)
    {
        std::string xfib("xfibers");
        xfib.push_back(var_vec[i]);
        importer.copy_elemental_solution(xfiber_system, xfib, xfib, step);
    }
    for( int i=0; i<3; i++)
    {
        std::string fib("sheets");
        fib.push_back(var_vec[i]);
        importer.copy_elemental_solution(sheets_system, fib, fib, step);
    }
 

   // mesh_fiber.prepare_for_use();
    std::cout << "norm of fiber solution = " << fiber_system.solution->l2_norm() << std::endl;
    std::cout << "norm of fiber current local solution = " << fiber_system.current_local_solution->l2_norm() << std::endl;

    std::cout << "norm of xfiber solution = " << xfiber_system.solution->l2_norm() << std::endl;
    std::cout << "norm of xfiber current local solution = " << xfiber_system.current_local_solution->l2_norm() << std::endl;

    std::cout << "norm of fiber solution = " << sheets_system.solution->l2_norm() << std::endl;
    std::cout << "norm of fiber current local solution = " << sheets_system.current_local_solution->l2_norm() << std::endl;

    // We may want to refine this mesh n times
    // if the original mesh is too coarse
    // read the number of times we will refine the mesh, 0 by default
    int num_refs = data("refs", 0);;
    if(num_refs > 0) {
    std::cout << "Refining the mesh " << num_refs << " times. " << std::endl;
    //Refine the mesh
    libMesh::MeshRefinement(mesh).uniformly_refine(num_refs);
    // reinitialization is projecting the solution,
    // old_solution, etc... vectors from the old mesh to
    // the current one.
    es_fiber.reinit ();
    //Finalize the mesh to run simulations
    mesh.prepare_for_use();
    }
    
    mesh.print_info();

    std::string output_file = data("output_file","output.e");
    libMesh::ExodusII_IO exporter(mesh);
    exporter.write(output_file);
    exporter.write_equation_systems(output_file, es_fiber);
    exporter.write_element_data(es_fiber);


	std::cout << "Good luck!" << std::endl;
    return 0;

}


