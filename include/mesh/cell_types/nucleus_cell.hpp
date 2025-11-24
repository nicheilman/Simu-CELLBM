#ifndef DEF_NUCLEUS_CELL
#define DEF_NUCLEUS_CELL

#include <cassert>

#include "../cell.hpp"


/*
    Use polymorphism to create the different cell types. Even though using IF statements is more efficient,
    polymorphism creates code that is more readable and maintainable. The IF statements are also more error prone.
*/



class nucleus_cell: public cell {

    private: 
        friend class cell_divider;

    protected:

	cell_ptr paired_cell_ptr_;

    public:

        //Make sure the nucleus_cell cannot be default instantiated
        nucleus_cell() = delete;                        //default constructor
        nucleus_cell(const nucleus_cell& c) = default;           //copy constructor
        nucleus_cell(nucleus_cell&& c) = delete;                //move constructor
        nucleus_cell& operator=(const nucleus_cell& c) = delete;//copy assignment operator
        nucleus_cell& operator=(nucleus_cell&& c) = delete;     //move assignment operator 


        //Construct the nucleus_cell with a node_lst and a face_lst
        nucleus_cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_lst, face_lst, cell_id, cell_type_parameters){}


        //Construct the nucleus_cell based on the node position and the node  ids of the faces
        nucleus_cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_position, face_node_ids, cell_id, cell_type_parameters){}


        //Convert a mesh object into a nucleus_cell object
        nucleus_cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(m, cell_id, cell_type_parameters){}

        cell_ptr get_cell_same_type(const mesh& m) noexcept(false) override{
            //Return a copy of the cell
            return std::make_shared<nucleus_cell>(m, cell_id_, cell_type_);
        }

        //Divide the cell when its volume has reached the maximum volume
        bool is_ready_to_divide() const noexcept override {
            assert(cell_type_ != nullptr);
            return volume_ >= division_volume_; 
        }
	
	cell_ptr get_paired_cell_ptr() const noexcept override {return paired_cell_ptr_;}
	void set_paired_cell_ptr(cell_ptr partner_cell_ptr) noexcept override {paired_cell_ptr_ = partner_cell_ptr; return;}
	void update_paired_cell_ptr(std::vector<cell_ptr> cell_lst) noexcept override {
		for(auto& cell_ : cell_lst){
			for(auto face_ : cell_->get_face_lst()){
			if(!face_.is_used()) continue;
			if(!(cell_->point_is_on_positive_side_of_face(face_.get_local_id(), centroid_))){ 
					set_paired_cell_ptr(cell_);
					break;
					};
				};
			};
		};
};

#endif
