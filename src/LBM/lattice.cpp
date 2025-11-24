
#include "lattice.hpp"
//#include "node.hpp"

lattice::lattice(

	double domain_size[3],
        int mesh_size[3], 
	double dt

)

{
    dt_ = dt;
    c = (domain_size[0]/(dt_*mesh_size[0]));
    mesh_size_ = mesh_size;
    domain_size_ = domain_size;
    origin[0] = -3.0e-5; origin[1] = -3.0e-5; origin[2] = -3.0e-5;

    num_nodes = mesh_size[0]*mesh_size[1]*mesh_size[2];
    //node_lst_.reserve(domain_size[0]*domain_size[1]*domain_size[2]);

    double mesh_x = domain_size[0]/(mesh_size[0]-1);
    double mesh_y = domain_size[1]/(mesh_size[1]-1);
    double mesh_z = domain_size[2]/(mesh_size[2]-1);

    for(int i=0; i<mesh_size[0]; i++){
        for(int j=0; j<mesh_size[1]; j++){
            for (int k=0; k<mesh_size[2]; k++){
		if((i==0) || (j==0) || (k==0) || (i==(mesh_size[0]-1)) || (j==(mesh_size[1]-1)) || (k==(mesh_size[2]-1)) )
		    node_lst_.push_back(std::make_shared<LBM_node>(i*mesh_x+origin[0],j*mesh_y+origin[1],k*mesh_z+origin[2],1, false, i*mesh_size[1]*mesh_size[2]+j*mesh_size[2]+k));
		else
		    node_lst_.push_back(std::make_shared<LBM_node>(i*mesh_x+origin[0],j*mesh_y+origin[1],k*mesh_z+origin[2],0, false, i*mesh_size[1]*mesh_size[2]+j*mesh_size[2]+k));
            }
        }
    }

    //node_lst_ = node_lst;
}
//----------------------------------------------------------------------------------------//

void lattice::stream(){

int cx, cy, cz, nb_index;

#pragma omp parallel for private(cx, cy, cz, nb_index)
for(auto& node_ : node_lst_ ){

//#pragma omp parallel for private(cx, cy, cz, nb_index)
for(int i=1; i<num_dir; i++){
    
    if(node_->get_wallflag() == 0){
	    nb_index = (num_nodes + node_->get_idx() - c_i[i][0]*mesh_size_[1]*mesh_size_[2] - c_i[i][1]*mesh_size_[2] - c_i[i][2]) % num_nodes; //Stride length to neighboring nodes
		node_->set_f( node_lst_[nb_index], i ); 
    continue;
    }

if(node_->get_wallflag() == 1){ 

    if( (node_->get_pos()[0] == origin[0]) && (c_i[i][0] == 1) ){cx = (mesh_size_[0] - 1);}
    else if( (almost_equal(node_->get_pos()[0], domain_size_[0]+origin[0])) && (c_i[i][0] == -1) ){cx = -1*(mesh_size_[0] - 1);}
    else{cx = -1*c_i[i][0];}

    if( (node_->get_pos()[1] == origin[1]) && (c_i[i][1] == 1) ){cy = (mesh_size_[1] - 1);}
    else if( (almost_equal(node_->get_pos()[1], domain_size_[1]+origin[1])) && (c_i[i][1] == -1) ){cy = -1*(mesh_size_[1] - 1);}
    else{cy = -1*c_i[i][1];}

    if( (node_->get_pos()[2] == origin[2]) && (c_i[i][2] == 1) ){cz = (mesh_size_[2] - 1);}
    else if( (almost_equal(node_->get_pos()[2], domain_size_[2]+origin[2])) && (c_i[i][2] == -1) ){cz = -1*(mesh_size_[2] - 1);}
    else{cz = -1*c_i[i][2];}
 
    nb_index = (node_->get_idx() + cx*mesh_size_[1]*mesh_size_[2] + cy*mesh_size_[2] + cz);
    node_->set_f( node_lst_[nb_index], i ); 


    if(almost_equal(node_->get_pos()[0], origin[0])){
        if(c_i[i][0] == 1){node_->set_f_wall(i, bd_flip(i, 0), 4.0e-11*c_i[i][2] * c );} 
	 }
    if(almost_equal(node_->get_pos()[0], domain_size_[0]+origin[0])){
        if(c_i[i][0] == -1){node_->set_f_wall(i, bd_flip(i, 0), -4.0e-11*c_i[i][2] * c );} 
	 }

/*
//-------------
if(node_->get_pos()[1] == origin){
        if(c_i[i][1] == 1){node_->set_f_wall(i, bd_flip(i, 1), 0.0*c_i[i][2]/c );}
	if(c_i[i][1] == 0){node_->set_f_wall(i, bd_flip(i, 2), 0.0*c_i[i][2]/c );}
         }
    if(node_->get_pos()[1] == domain_size_[1]+origin){
        if(c_i[i][1] == -1){node_->set_f_wall(i, bd_flip(i, 1), -0.0*c_i[i][2]/c );}
	if(c_i[i][1] == 0){node_->set_f_wall(i, bd_flip(i, 2), -0.0*c_i[i][2]/c );}
         }
//-------------
*/
/*/-------------

if(node_->get_pos()[2] == origin[2]){
        if(c_i[i][2] == 1){node_->set_f_p(i, node_lst_[node_->get_idx()]);}
         }
if(node_->get_pos()[2] == domain_size_[2]+origin[2]){
        if(c_i[i][2] == -1){node_->set_f(node_lst_[node_->get_idx()-1], i);}
         }

//-------------*/

      }
   }
}
return;

}

