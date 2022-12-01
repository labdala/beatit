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

int main(int argc, char ** argv)
{
    using namespace libMesh;
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    Mesh mesh_fiber(init.comm());
    ExodusII_IO importer(mesh_fiber);

    importer.read("RESILIENT_FIBERS_CLEAN.e");

    mesh_fiber.prepare_for_use();
    mesh_fiber.print_info();
    
    libMesh::EquationSystems  es_fiber(mesh_fiber);
    typedef libMesh::ExplicitSystem  FiberSystem;
    FiberSystem& fiber_system = es_fiber.add_system<FiberSystem>("fibers");
    fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
    fiber_system.init();
    auto& f_v = fiber_system.solution;

    unsigned int step=1;
    std::cout << "Available variables in exodus file" << std::endl;
    auto& elem_var_names = importer.get_elem_var_names();
    for (auto && var : elem_var_names) std::cout << var << std::endl;

    //Read fiber data into mesh_fiber
    importer.copy_elemental_solution(fiber_system, "fibersx", "fibersx", step);
    importer.copy_elemental_solution(fiber_system, "fibersy", "fibersy", step);
    importer.copy_elemental_solution(fiber_system, "fibersz", "fibersz", step);
    fiber_system.print_info();
    
    auto norm = fiber_system.solution->l2_norm();    
    std::cout << "norm = " << norm << std::endl; 

    ExodusII_IO exporter(mesh_fiber);
    exporter.write("mesh_fiber_1block.e");

    return 0;
}
