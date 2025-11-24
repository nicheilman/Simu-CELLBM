#include "global_configuration.hpp"
#include "cell.hpp"

#include <utility>

/* class OurContactModel {

public:

OurContactModel() noexcept(false);

std::tuple<double, double> get_contact_strength(cell_type_param_ptr c1_cell_type, cell_type_param_ptr c2_cell_type) noexcept;

}; */

struct OurContactModel {
 public:

  static std::pair<double, double> get_contact_strength(const short int cell_type_id1, const short int cell_type_id2) {
     return contact_strengths_.at(std::make_pair(cell_type_id1, cell_type_id2));
  }

private:
  static inline std::map<std::pair<const short int, const short int>, std::pair<double, double>> contact_strengths_ = 
  {
    {std::make_pair(0, 0), std::make_pair(1.0e11, 1.0e9)},
    {std::make_pair(0, 2), std::make_pair(0.0e8, 0.0e8)},
    {std::make_pair(2, 2), std::make_pair(0.0, 0.0)},
    {std::make_pair(2, 0), std::make_pair(0.0e8, -5.0e8)},
    {std::make_pair(1, 1), std::make_pair(0.0e7, 1.0e8)},
    {std::make_pair(0, 1), std::make_pair(0.0e8, 1.0e8)},
    {std::make_pair(1, 0), std::make_pair(0.0e8, 1.0e8)},
    {std::make_pair(1, 2), std::make_pair(0.0e8, 0.0e8)},
    {std::make_pair(2, 1), std::make_pair(0.0e8, -5.0e8)},
    
    {std::make_pair(3, 0), std::make_pair(0.0e8, -1.0e9)},
    {std::make_pair(0, 3), std::make_pair(0.0e9, 0.0e9)},
    {std::make_pair(3, 1), std::make_pair(0.0e8, 0.0e8)},
    {std::make_pair(1, 3), std::make_pair(0.0e8, 0.0e8)},
    {std::make_pair(3, 2), std::make_pair(0.0e8, 0.0e8)},
    {std::make_pair(2, 3), std::make_pair(0.0e8, 0.0e8)},
    {std::make_pair(3, 3), std::make_pair(0.0e7, 1.0e10)}
  };  // (celltype id , cell type id) -> (adherence_strength, repulsion_strength)
};  // OurContactModel

