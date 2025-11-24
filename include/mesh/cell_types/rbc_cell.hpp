#ifndef DEF_RBC_CELL
#define DEF_RBC_CELL

#include <cassert>

#include "../cell.hpp"

#include "utils.hpp"


/*
    Use polymorphism to create the different cell types. Even though using IF statements is more efficient,
    polymorphism creates code that is more readable and maintainable. The IF statements are also more error prone.
*/


class rbc_cell: public cell {

    private: 
        friend class cell_divider;
        friend class cell_tester;

        
    protected:


    public:

        //------------------------------------------------------------------------------------------------------
	
	

	//------------------------------------------------------------------------------------------------------

        //------------------------------------------------------------------------------------------------------
        //Make sure the rbc_cell cannot be default instantiated
        rbc_cell() = delete;                        //default constructor
        rbc_cell(const rbc_cell& c) = default;           //copy constructor
        rbc_cell(rbc_cell&& c) = delete;                //move constructor
        rbc_cell& operator=(const rbc_cell& c) = delete;//copy assignment operator
        rbc_cell& operator=(rbc_cell&& c) = delete;     //move assignment operator 
        //------------------------------------------------------------------------------------------------------


        //------------------------------------------------------------------------------------------------------ 
        //Construct the rbc_cell with a node_lst and a face_lst
        rbc_cell(
            const std::vector<node>& node_lst,
            const std::vector<face>& face_lst,
            unsigned int cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_lst, face_lst, cell_id, cell_type_parameters){}


        //Construct the rbc_cell based on the node position and the node  ids of the faces
        rbc_cell(
            const std::vector<double>&  node_position,
            const std::vector<unsigned>& face_node_ids,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(node_position, face_node_ids, cell_id, cell_type_parameters){}


        //Convert a mesh object into a rbc_cell object
        rbc_cell(
            const mesh& m,
            unsigned cell_id,
            cell_type_param_ptr cell_type_parameters
        ) noexcept : cell(m, cell_id, cell_type_parameters){}
        //------------------------------------------------------------------------------------------------------


        //------------------------------------------------------------------------------------------------------
        //Divide the cell when its volume has reached the maximum volume
        bool is_ready_to_divide() const noexcept override {
            assert(cell_type_ != nullptr);
            return volume_ >= division_volume_; 
        }
        //------------------------------------------------------------------------------------------------------

        //------------------------------------------------------------------------------------------------------
        cell_ptr get_cell_same_type(const mesh& m) noexcept(false) override{
            //Return a copy of the cell
            return std::make_shared<rbc_cell>(m, cell_id_, cell_type_);
        }
        //------------------------------------------------------------------------------------------------------

};

#endif
