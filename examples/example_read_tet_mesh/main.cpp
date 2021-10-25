#include <iostream>

#include "Util/IO/io.hpp"


#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <string>

#include "libmesh/mesh.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/node.h"
#include "libmesh/point.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/unv_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/dof_map.h"

int main (int argc, char ** argv)
{
     //BeatIt::printBanner(std::cout);


    std::cout << " HI 0" <<std::endl;


      using namespace libMesh;
      // Initialize libraries, like in example 2.
      LibMeshInit init (argc, argv, MPI_COMM_WORLD);
      libMesh::Mesh mesh(init.comm());
      mesh.read("./bin/03patient/03meshLA.e", nullptr, true);
      std::ifstream fibers_list("./bin/03patient/03meshLA.lon");
      std::string line;
      double x,y,z;


      //create .e file using .elem and .pts
      if(false){
      mesh.set_mesh_dimension(3);
      mesh.set_spatial_dimension(3);
      mesh.reserve_elem(79067310);
      mesh.reserve_nodes(14352781);


      std::cout << " HI 0a" <<std::endl;


      std::ifstream nodes_list("./bin/02patient/02meshLA.pts");
      std::ifstream elements_list("./bin/02patient/02meshLA.elem");

      std::ofstream elements_file;
      elements_file.open("elements.txt");
      std::ofstream nodes_file;
      nodes_file.open("nodes.txt");


      // string that will store the line
      std::string line;

      long double x, y, z;


      std::cout << " HI 0b" <<std::endl;


      //int c = 0;

      // number of the node
      int n = 0;

      //get line returns true if not at the end of the file

      std::cout << " HI 1" <<std::endl;

      if(nodes_list.is_open())
      {
    	    getline (nodes_list,line); // just to skip the first line in the file
    	    std::cout << " HI 2" <<std::endl;
			while ( getline (nodes_list,line) )
			{
				// input string  > stream of strings (list that breaks down when there is a space)
			    std::istringstream ss(line);
				//line will have
			    ss >> x >> y >> z;
			    libMesh::Node node(x,y,z);
			    mesh.add_point(libMesh::Point(x,y,z),n);
                n++;
                elements_file << x << " " << y << " " << z << std::endl;
                 /*std::cout << line << std::endl;
                //test
                if(n > 10) 	break;*/
			}
      }


      std::cout << " HI 3, n=" << n <<std::endl;

      n=0;
      int n1, n2, n3, n0, nid;
      char T,t;
      if(elements_list.is_open())
      {
            getline (elements_list,line);
            while ( getline (elements_list,line) )
            {
                std::istringstream ss(line);


                libMesh::Elem * elem = mesh.add_elem(new libMesh::Tet4);
                // In libMesh:
                /*                *   TET4:
                *         3
                *         o                 zeta
                *        /|\                 ^
                *       / | \                |
                *      /  |  \               |
                *   0 o...|...o 2            o---> eta
                *      \  |  /                \
                *       \ | /                  \
                *        \|/                    xi (out of page)
                *         o
                *         1
                *
                */
                ss >> T >> t >> n0 >> n1 >> n2 >> n3 >> nid;
                elem->set_node(0) = mesh.node_ptr( n0 );
                elem->set_node(1) = mesh.node_ptr( n1 );
                elem->set_node(2) = mesh.node_ptr( n2 );
                elem->set_node(3) = mesh.node_ptr( n3 );
                elem->subdomain_id() = nid;

                n++;
                /*std::cout << line << std::endl;
                std::cout << " T="<< T <<" t="<<  t <<" n0="<< n0 <<" n1="<< n1 <<" n2="<< n2 <<" n3="<< n3 << std::endl;
                //test
                if(n > 10) 	break;*/

                nodes_file << n0<< " "<< n1<< " "<< n2<< " "<< n3<< std::endl;
                // TODO:: Check sides
                //mesh.boundary_info->add_side()
            }
      }

      std::cout << " HI 4, n="<< n <<std::endl;
      nodes_list.close();
      elements_list.close();

      mesh.prepare_for_use (/*skip_renumber =*/ true);
      }
      mesh.print_info();

//      libMesh::ExodusII_IO(mesh).write("./bin/02patient/02meshLA.e");

      std::cout << "CREATING EQUATION SYSTEMS"<< std::endl;
      libMesh::EquationSystems equation_systems(mesh);
      typedef libMesh::ExplicitSystem  FiberSystem;
         FiberSystem& fiber_system = equation_systems.add_system<FiberSystem>("fibers");
         fiber_system.add_variable("fibersx", libMesh::CONSTANT, libMesh::MONOMIAL);
         fiber_system.add_variable("fibersy", libMesh::CONSTANT, libMesh::MONOMIAL);
         fiber_system.add_variable("fibersz", libMesh::CONSTANT, libMesh::MONOMIAL);
         equation_systems.init();
         auto& f_v = fiber_system.solution;

         getline (fibers_list,line); // just to skip the first line in the file

         std::ofstream myfile;
         myfile.open("fibers.txt");

         int n = 1;
         for (const auto& elem: mesh.active_local_element_ptr_range()) {

        	 	 getline(fibers_list,line);

        	 	 std::vector<dof_id_type> fibers_dof_indices; //putting this back to the system so that we can export it
 			     std::istringstream ss(line);
 				 //line will have
 			     ss >> x >> y >> z;

                            /* x=1;
                             y=1;
                             z=1;*/
        	     // assign fiber coordinate to vector f0

 		         // prevent bug to happen - stop when we have gone through all lines in .lon
                 if(n==46154318){
                     //    z=0;
                	 libMesh::Point f0_new(x,y,z);
                     for (int idim = 0; idim < 2; ++idim)
                     {
                         fiber_system.get_dof_map().dof_indices(elem, fibers_dof_indices);
                         fiber_system.solution ->set(fibers_dof_indices[idim], f0_new(idim));
                         myfile << f0_new(idim)<< " ";
                     }
                     myfile << std::endl;
                     std::cout << "THIS IS THE VALUE OF N INSIDE THE LOOP" << n <<std::endl;
                     break;
                 }
                 else{
                	 libMesh::Point f0(x,y,z);
                     for (int idim = 0; idim < 3; ++idim)
                     {
                         fiber_system.get_dof_map().dof_indices(elem, fibers_dof_indices);
                         fiber_system.solution ->set(fibers_dof_indices[idim], f0(idim));
                         myfile << f0(idim)<< " ";
                     }
                     myfile << std::endl;
                 }

                 n++;
         }

         myfile.close();
         std::cout << "THIS IS THE VALUE OF N OUTSIDE THE LOOP" << n <<std::endl;

         std::cout << "READY TO WRITE THE FIBERS"<< std::endl;
         ExodusII_IO nemesis_exporter(mesh); // creates the exporter
         std::cout << "RIGHT BEFORE WRITE"<< std::endl;
         nemesis_exporter.write("03meshLA_fibers.e"); // saves the variables at the nodes
         std::cout << "RIGHT BEFORE WRITE ELEMENT DATA"<< std::endl;


        nemesis_exporter.write_element_data(equation_systems); // this saves the variables at the centroid of the elements

        std::cout << "fim" << std::endl;

        fibers_list.close();

      return 0;
}
