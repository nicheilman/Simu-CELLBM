
#include "header.hpp"
#include "LBM_node.hpp"
#include "IB_node.hpp"
#include "mesh_writer.hpp"
  
       void MeshWriter::writeVTK(const std::string &filename,
			 const int mesh[3],
                         const std::vector<std::shared_ptr<LBM_node>> &nodes)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            throw std::runtime_error("Unable to open file: " + filename);
        }

	ofs << "# vtk DataFile Version 4.2" << std::endl;
	ofs << "3D Mesh" << std::endl;
	ofs << "ASCII" << std::endl;
	ofs << "DATASET STRUCTURED_GRID" << std::endl;
	ofs << "DIMENSIONS " << mesh[0] << " " << mesh[1] << " " << mesh[2] << " " << std::endl;
       
        ofs << "POINTS " << nodes.size() << " float\n";
        ofs << std::scientific << std::setprecision(4);
/*        for (const auto &node : nodes) {
            ofs << node->get_pos()[0] << " " << node->get_pos()[1] << " " << node->get_pos()[2] << " " << std::endl;
        }
  */

	for(int k=0; k<mesh[2]; k++){
            for(int j=0; j<mesh[1]; j++){
        	for (int i=0; i<mesh[0]; i++){	
	    	    auto node = nodes[i*mesh[1]*mesh[2]+j*mesh[2]+k];
            	    ofs << node->get_pos()[0] << " " << node->get_pos()[1] << " " << node->get_pos()[2] << " " << std::endl;
        	}
	    }
	}
      
	ofs << "POINT_DATA " << nodes.size() << std::endl;
	ofs << "SCALARS " << "density " << "float" << std::endl;
	ofs << "LOOKUP_TABLE " << "default" <<std::endl;
        for(int k=0; k<mesh[2]; k++){
            for(int j=0; j<mesh[1]; j++){
                for (int i=0; i<mesh[0]; i++){
                    auto node = nodes[i*mesh[1]*mesh[2]+j*mesh[2]+k];
		    ofs << node->get_m()[0] << std::endl;
        	}
	    }
	}

        ofs << "SCALARS " << "internal " << "int" << std::endl;
        ofs << "LOOKUP_TABLE " << "default" <<std::endl;
        for(int k=0; k<mesh[2]; k++){
            for(int j=0; j<mesh[1]; j++){
                for (int i=0; i<mesh[0]; i++){
                    auto node = nodes[i*mesh[1]*mesh[2]+j*mesh[2]+k];
                    ofs << node->get_internal() << std::endl;
                }
            }
        }


        ofs << "VECTORS " << "Velocity " << "float" << std::endl;
        for(int k=0; k<mesh[2]; k++){
            for(int j=0; j<mesh[1]; j++){
                for (int i=0; i<mesh[0]; i++){
                    auto node = nodes[i*mesh[1]*mesh[2]+j*mesh[2]+k];
	    	    ofs << node->get_velocity()[0] << " " << node->get_velocity()[1] << " " << node->get_velocity()[2] << " " << std::endl;
        	}
	    }
	}


	ofs.close();
    };


	void MeshWriter::writeflagVTK(const std::string &filename,
                         const int mesh[3],
                         const std::vector<std::shared_ptr<IB_node>> &nodes)
    {
        std::ofstream ofs(filename);
        if (!ofs) {
            throw std::runtime_error("Unable to open file: " + filename);
        }

        ofs << "# vtk DataFile Version 4.2" << std::endl;
        ofs << "3D Mesh" << std::endl;
        ofs << "ASCII" << std::endl;
        ofs << "DATASET STRUCTURED_GRID" << std::endl;
        ofs << "DIMENSIONS " << mesh[1]+1 << " " << mesh[0]+1 << " " << mesh[2]+1 << " " << std::endl;

        ofs << "POINTS " << nodes.size() << " float\n";
        ofs << std::scientific << std::setprecision(4);
/*        for (const auto &node : nodes) {
            ofs << node->get_pos()[0] << " " << node->get_pos()[1] << " " << node->get_pos()[2] << " " << std::endl;
        }
  */


            for(int i=0; i<mesh[0]+1; i++){
                for (int j=0; j<mesh[1]+1; j++){
                    auto node = nodes[i*(mesh[1]+1)+j];
                    ofs << node->get_IB_pos()[0] << " " << node->get_IB_pos()[1] << " " << node->get_IB_pos()[2] << " " << std::endl;
                }
            }
        
}
