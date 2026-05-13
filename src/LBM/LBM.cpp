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

//TIMER
double start, IB_start, NB_start, stream_start, coll_start; 
double end, IB_end, NB_end, stream_end, coll_end; 
start = omp_get_wtime(); 

cout << fixed << setprecision(10);
cout << "CellBM is running!" << endl;

array<double, 6> aabb;
double velo_interp[3];

vec3 F;
std::vector<double> zero = {0., 0., 0.};
NB_start = omp_get_wtime();
#pragma omp parallel for
for(auto& node_ : L_->get_nodes() ){ //  

    node_->ftom(L_, 0);
    
    node_->calc_eq();

//NB_end = omp_get_wtime();
//NB_accumulate = NB_end - NB_start;
//printf("Equil took %f seconds\n", NB_accumulate);

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
         if(node_->is_internal(cell_,aabb, t)){
		node_->set_internal(1);
		node_->set_velocity(zero); break;}
	 else continue;
    /*node_->set_internal(0);*/};}
    else{
	if(node_->is_internal(cell_,aabb, t)) continue;   //node_->set_internal(0, velo); break;}
        else node_->set_internal(1);
        };
    };

}
else if( int(t/dt_) % 10 == 0){
    for(auto& cell_ : cell_lst){
    aabb = cell_->get_aabb();

if(cell_->get_cell_type_id() != 2){
    if( (node_->get_pos()[0] <= aabb[0])&&(node_->get_pos()[0] >= aabb[3]) &&
    (node_->get_pos()[1] <= aabb[1])&&(node_->get_pos()[1] >= aabb[4]) &&
    (node_->get_pos()[2] <= aabb[2])&&(node_->get_pos()[2] >= aabb[5])){
         node_->set_internal(0);
	    };
	};
    };
};

//------------------------------

	if(node_->get_internal()){node_->set_velocity(zero); continue;}
}
NB_end = omp_get_wtime();
NB_accumulate += NB_end - NB_start;
printf("Equil took %f seconds\n", NB_accumulate);

//TIMER
IB_start = omp_get_wtime();

// Only fluid nodes near the IB need to be looked at //
// Will be best to use Bounding Box //
//#pragma omp parallel for
for(auto& cell_ : cell_lst){
    aabb = cell_->get_aabb();
    
if(cell_->get_cell_type_id() == 2) continue;

    for(node& IB_node_ : cell_->node_lst_ ) {
if(!IB_node_.is_used()) continue;

force_[0] = 0.; force_[1] = 0.; force_[2] = 0.;
for(int i =0; i<3; i++) velo_interp[i] = 0.0;

    for(auto& node_ : L_->IB_neighbors(IB_node_)){
	if(IB_node_.neighbor_base_check()){
//TIMER
//NB_start = omp_get_wtime();
	if(node_->is_internal(cell_,aabb, t)){
                node_->set_internal(1);}
         else node_->set_internal(0);
//NB_end = omp_get_wtime();
//NB_accumulate += NB_end - NB_start;
		}
if(node_->get_internal() == 1) continue;

        dist_x = (node_->get_pos()[0] - IB_node_.pos().dx()) ;
        dist_y = (node_->get_pos()[1] - IB_node_.pos().dy()) ;
        dist_z = (node_->get_pos()[2] - IB_node_.pos().dz()) ;

        dist = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;

        if(dist <= 4*mesh_space*mesh_space) {
            kernel = (1 - abs( dist_x/mesh_space/2)) * (1 - abs( dist_y/mesh_space/2)) * (1 - abs( dist_z/mesh_space/2));
        for(int i=0; i<3; i++) velo_interp[i] += node_->get_m()[i+1] / node_->get_m()[0] * kernel ;    
	}
    }

    for(int i=0; i<3; i++) force_[i] = ( (IB_node_.momentum()).to_array()[i] - velo_interp[i] );

	    F.reset(-1.*force_[0], -1.*force_[1], -1.*force_[2]);
	    IB_node_.add_force(F);
//(IB_node_.momentum()).print();
//std::cout << velo_interp[0] << ", " << velo_interp[1] << ", "<< velo_interp[2] << std::endl; 
//F.print();
F.reset();
//IB_node_.momentum().print();
    }
}
//TIMER
IB_end = omp_get_wtime();
IB_accumulate += IB_end - IB_start;
//printf("Internal took %f seconds\n", NB_accumulate);
printf("IB took %f seconds\n", IB_accumulate);

//----------------------------------------------------------------------
force_[0] = 0.; force_[1] = 0.; force_[2] = 0.;
coll_start = omp_get_wtime();
#pragma omp parallel for
for(auto& node_ : L_->get_nodes() ){

    node_->collision(dt_, force_, node_->get_internal());

    node_->ftom(L_, 1);

}
coll_end = omp_get_wtime();
coll_accumulate += coll_end - coll_start;
printf("Collision took %f seconds\n", coll_accumulate);

stream_start = omp_get_wtime();
    L_->stream();
stream_end = omp_get_wtime();
stream_accumulate += stream_end - stream_start;
printf("Stream took %f seconds\n", stream_accumulate);

for(auto& node_ : L_->get_nodes()){
    node_->update_f();
    node_->set_int_mem(0);
    }
end = omp_get_wtime();
//TIMER

accumulate += end - start;
printf("IB-LBM took %f seconds\n", accumulate);

};

void LBM::write(double t){


if(int(t/dt_) % 5 == -1){
MeshWriter::writeflagVTK("test2/test_flag_"+to_string((t))+".vtk", flag_mesh, flag_->get_IB_nodes() );
    }
if(int(t/dt_) % 10 == 0){
MeshWriter::writeVTK("/mnt/scratch/heilman7/simulation_results/cell_17_clustered/test_mesh_"+to_string(int(t/dt_/10))+".vtk", mesh, L_->get_nodes() );
    }


};
