#include "libmesh/exodusII_io.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"

int main(int argc, char ** argv)
{
    using namespace libMesh;
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    Mesh mesh_fiber(init.comm());
    ExodusII_IO importer(mesh_fiber);

    importer.read("mesh_with_fibers.e");
    //mesh_fiber.read("mesh_with_fibers.e");
    mesh_fiber.prepare_for_use();
    mesh_fiber.print_info();

    // Looping over elements of the fiber mesh
    MeshBase::const_element_iterator el = mesh_fiber.elements_begin();
	const MeshBase::const_element_iterator end_el = mesh_fiber.elements_end();
    //MeshBase::const_element_iterator el = mesh_fiber.local_elements_begin();
	//const MeshBase::const_element_iterator end_el = mesh_fiber.local_elements_end();
    for(; el!=end_el; el++)  (*el)->subdomain_id()=1;
    
    ExodusII_IO exporter(mesh_fiber);
    exporter.write("mesh_fiber_1block.e");

    return 0;
}
