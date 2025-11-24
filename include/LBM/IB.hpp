
#include "header.hpp"
#include "IB_node.hpp"

class LBM_node;

class IB{

    protected:
    
    std::vector<std::shared_ptr<IB_node>> IB_node_lst_ = {};
    int* mesh_size_;
    double* com_;
    double rad_;

    friend class IB_node;

    public:

    IB() = delete;                         //default constructor
    IB(const IB& c) = default;           //copy constructor
    IB(IB&& c) = delete;                 //move constructor
    IB& operator=(const IB& c) = delete; //copy assignment operator
    IB& operator=(IB&& c) = delete;      //move assignment operator 

    IB(int mesh_size[3], double com[3], double rad);

    std::vector<std::shared_ptr<IB_node>> get_IB_nodes(){return IB_node_lst_;};
    bool in_bbox(std::shared_ptr<LBM_node> node, double mesh_space);
    bool is_internal(std::shared_ptr<LBM_node> node);

    void update_pos(std::vector<double> pos){
	for(int i=0; i<3; i++) com_[i] = com_[i] + pos[i];
	for(auto& IB_node_ : IB_node_lst_) IB_node_->update_node_pos(pos);
    }

};



