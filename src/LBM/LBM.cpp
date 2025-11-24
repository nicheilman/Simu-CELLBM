//#include "header.hpp"
//#include "lattice.hpp"
//#include "IB.hpp"
//#include "mesh_writer.hpp"
#include "LBM.hpp"
#include "../math_modules/vec3.hpp"
using namespace std;

LBM::LBM(double dt, std::vector<cell_ptr> cell_lst, double& min_edge_len) : dt_(dt), cell_lst_(cell_lst)//, mesh_space(min_edge_len)

{
    for(auto& node_ : L_->get_nodes()){
        node_->calc_eq();
	node_->set_meq();
        node_->ftom(L_, 1);
        for(int i=0; i<19; i++)node_->set_f_p(i, node_);
        node_->update_f();
    }
};

void LBM::run(double t, const std::vector<cell_ptr>& cell_lst){

cout << fixed << setprecision(10);
cout << "CellBM is running!" << endl;

array<double, 6> aabb;
//for(auto& cell_ : cell_lst) aabb = cell_->get_aabb();

vec3 F;

for(auto& node_ : L_->get_nodes() ){ //  

    node_->ftom(L_, 0);

    node_->calc_eq();

force_[0] = 0.; force_[1] = 0.; force_[2] = 0.; kernel = 0.;

//----------------Needs to be converted to own function(s)------------------------

//------------------------------
if(t == 0.){
for(auto& cell_ : cell_lst){
    aabb = cell_->get_aabb();


if(cell_->get_cell_type_id() != 2){
    if( (node_->get_pos()[0] >= aabb[0]-2*mesh_space)&&(node_->get_pos()[0] <= aabb[3]+2*mesh_space) &&
    (node_->get_pos()[1] >= aabb[1]-2*mesh_space)&&(node_->get_pos()[1] <= aabb[4]+2*mesh_space) &&
    (node_->get_pos()[2] >= aabb[2]-2*mesh_space)&&(node_->get_pos()[2] <= aabb[5]+2*mesh_space)){
         if(node_->is_internal(cell_,aabb)){
		node_->set_internal(1, velo); break;}
	 else continue;
    node_->set_internal(0, velo);};}
    else{
	if(node_->is_internal(cell_,aabb)) continue;   //node_->set_internal(0, velo); break;}
        else node_->set_internal(1, velo);
        };
    };

}
else if( int(t/dt_) % 10 == 0){
    for(auto& cell_ : cell_lst){
    aabb = cell_->get_aabb();

if(cell_->get_cell_type_id() != 2){
    if( (node_->get_pos()[0] >= aabb[0]-2*mesh_space)&&(node_->get_pos()[0] <= aabb[3]+2*mesh_space) &&
    (node_->get_pos()[1] >= aabb[1]-2*mesh_space)&&(node_->get_pos()[1] <= aabb[4]+2*mesh_space) &&
    (node_->get_pos()[2] >= aabb[2]-2*mesh_space)&&(node_->get_pos()[2] <= aabb[5]+2*mesh_space)){
         if(node_->is_internal(cell_,aabb)){
                node_->set_internal(1, velo); break;}
         else node_->set_internal(0, velo);
	    };
	};
    };

};

//------------------------------

	if(node_->get_internal()) continue;

// Only fluid nodes near the IB need to be looked at //
// Will be best to use Bounding Box //
for(auto& cell_ : cell_lst){
    aabb = cell_->get_aabb();

if(cell_->get_cell_type_id() == 2) continue;
// std::cout<< "Taking a look at Cell:"<< cell_->get_id() << std::endl;

    if( (node_->get_pos()[0] >= aabb[0]-2*mesh_space)&&(node_->get_pos()[0] <= aabb[3]+2*mesh_space) && 
    (node_->get_pos()[1] >= aabb[1]-2*mesh_space)&&(node_->get_pos()[1] <= aabb[4]+2*mesh_space) &&
    (node_->get_pos()[2] >= aabb[2]-2*mesh_space)&&(node_->get_pos()[2] <= aabb[5]+2*mesh_space)){
//  for(auto& cell_ : cell_lst){
//	if(t >= 0.){node_->set_internal(node_->is_internal(cell_,aabb), velo);}
    for(node& IB_node_ : cell_->node_lst_ ) {

        dist_x = (node_->get_pos()[0] - IB_node_.pos().dx()) ;
        dist_y = (node_->get_pos()[1] - IB_node_.pos().dy()) ;
        dist_z = (node_->get_pos()[2] - IB_node_.pos().dz()) ;

/* Internal nodes are treated as solid, Do not need IB consideration */
        //if(t == 0.){node_->set_internal(node_->is_internal(cell_,aabb), velo);}  //{if(node_->is_internal(cell_, aabb)) node_->set_internal(velo);}

        dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
        dist = sqrt(dist);
        if(dist <= 2*mesh_space ){
		//if(node_->get_internal()) continue;
            kernel = (1 - abs( dist_x/mesh_space/2)) * (1 - abs( dist_y/mesh_space/2)) * (1 - abs( dist_z/mesh_space/2));
            for(int i=0; i<3; i++) force_[i] = (IB_node_.momentum().to_array()[i]-1*node_->get_m()[i+1]/node_->get_m()[0]) * kernel * dt_ / mesh_space  ;
	    F = vec3(-1.*force_[0], -1.*force_[1], -1.*force_[2]);
	    IB_node_.add_force(F/dt_*mesh_space);
        }
    }
  }
}

//----------------------------------------------------------------------

    node_->collision(dt_, force_, node_->get_internal());

    node_->ftom(L_, 1);

}
    L_->stream();

for(auto& node_ : L_->get_nodes()){
    node_->update_f();
    }

};

void LBM::write(double t){


if(int(t/dt_) % 5 == -1){
MeshWriter::writeflagVTK("test2/test_flag_"+to_string((t))+".vtk", flag_mesh, flag_->get_IB_nodes() );
    }
if(int(t/dt_) % 10 == -1){
MeshWriter::writeVTK("test2/test_mesh_"+to_string(int(t/dt_/10))+".vtk", mesh, L_->get_nodes() );
    }


};
