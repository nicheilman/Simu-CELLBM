//#include "header.hpp"
#include "IB.hpp"
#include "LBM_node.hpp"

IB::IB(

    int mesh_size[3],
    double com[3],
    double rad
){

double PI = 4*std::atan(1);
mesh_size_ = mesh_size;
com_ = com;
rad_ = rad;
for(int i=0; i<mesh_size_[0]+1; i++){
        for(int j=0; j<mesh_size_[1]+1; j++){
            for (int k=0; k<mesh_size_[2]+1; k++){
		IB_node_lst_.push_back(std::make_shared<IB_node>( i*(0.2/mesh_size[0])+com_[0]-0.0, rad_*cos(j*(2*PI/mesh_size[1]) )+com_[1], -1*rad_*sin(j*(2*PI/mesh_size[1]) )+com_[2], i*mesh_size_[1]+j));

            }
        }
    }
}

bool IB::in_bbox(std::shared_ptr<LBM_node> node, double mesh_space){
return (node->get_pos()[2] > com_[2] - rad_ - 2*mesh_space) && (node->get_pos()[2] < com_[2] + rad_ + 2*mesh_space)
        && (node->get_pos()[1] > com_[1] - rad_ - 2*mesh_space) && (node->get_pos()[1] < com_[1] + rad_ + 2*mesh_space);}

bool IB::is_internal(std::shared_ptr<LBM_node> node){
return ( ((std::pow( node->get_pos()[1] - com_[1],2) + std::pow(node->get_pos()[2] - com_[2],2) <= rad_*rad_) && !node->get_internal())
	|| ((std::pow( node->get_pos()[1] - com_[1],2) + std::pow(node->get_pos()[2] - com_[2],2) > rad_*rad_) && node->get_internal()) );}

