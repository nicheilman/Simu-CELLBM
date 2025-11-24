
#include "header.hpp"


class IB_node: public std::enable_shared_from_this<IB_node>{

    protected:
    
    double dx_, dy_, dz_;
    double p_, u_[3]; 
    int idx_;


    public:

    IB_node() = default;                         //default constructor
    IB_node(const IB_node& c) = default;           //copy constructor
    IB_node(IB_node&& c) = default;                 //move constructor
    IB_node& operator=(const IB_node& c) = default; //copy assignment operator
    IB_node& operator=(IB_node&& c) = default;      //move assignment operator 


    IB_node( double dx, double dy, double dz, int idx);

    std::vector<double> get_IB_pos(){return {dx_, dy_, dz_};};
    int get_idx(){return idx_;};

    void set_pos(std::vector<double> pos){dx_ = pos[0]; dy_ = pos[1]; dz_ = pos[2]; return;}
    void update_node_pos(std::vector<double> shift){dx_ = dx_+shift[0]; dy_ = dy_+shift[1]; dz_ = dz_+shift[2]; return;}

};



