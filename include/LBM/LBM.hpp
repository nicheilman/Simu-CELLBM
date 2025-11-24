#include "header.hpp"
#include "lattice.hpp"
#include "IB.hpp"
#include "mesh_writer.hpp"
#include "../mesh/cell.hpp"
using namespace std;

class LBM{

private:

    int mesh_size = 21;
    double dimensions[3] = {6.0e-5, 6.0e-5, 6.0e-5};
    double dt_;
 // ------------------------------------------- //

    double IB_CoM[3] = {1.0, 0.5*dimensions[1], 0.1*dimensions[2]};
    double IB_r = 0.15;
    double min_dim = min(min(dimensions[0], dimensions[1]), dimensions[2]);
    int mesh[3] = {mesh_size*int(dimensions[0]/min_dim), mesh_size*int(dimensions[1]/min_dim), mesh_size*int(dimensions[2]/min_dim)};
    int flag_mesh[3] = {0, 0, 0}; //2*mesh[0], 2*mesh[1], 0}; // Must be ~2x the fluid mesh spacing
    double force_[3]; // = {0., 0., 0.};
    vector<double> velo = {0., 0., 1.0};
    vector<double> shift = {0., 0., 0.};
    double mesh_space = dimensions[0]/mesh_size ;
    double dist, dist_x, dist_y, dist_z, kernel;

    shared_ptr<lattice> L_ = make_shared<lattice>(dimensions, mesh, dt_);
    shared_ptr<IB> flag_ = make_shared<IB>(flag_mesh, IB_CoM, IB_r);
    std::vector<cell_ptr> cell_lst_;

friend class cell;
friend class node;

public:

    LBM() = delete;                         //default constructor
    LBM(const LBM& c) = default;           //copy constructor
    LBM(LBM&& c) = delete;                 //move constructor
    LBM& operator=(const LBM& c) = delete; //copy assignment operator
    LBM& operator=(LBM&& c) = delete;      //move assignment operator 

    LBM(double dt, std::vector<cell_ptr> cell_lst, double& min_edge_len);

    void run(double t, const std::vector<cell_ptr>& cell_lst);
    void write(double t);

};

