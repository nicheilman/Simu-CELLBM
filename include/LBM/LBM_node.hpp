
#include "header.hpp"
#include "../mesh/cell.hpp"
#include "../math_modules/vec3.hpp"

class lattice;
class IB;

class LBM_node: public std::enable_shared_from_this<LBM_node>{

    protected:
    
    const int velo_dim = 19;
    std::array<double, 19> f_, f_tmp, m_, meq, lambda, f_eq; 
    double dx_, dy_, dz_;
    double p_, u_[3]; 
    double m2;
    bool internal_;
    int wallflag_, idx_;

    static constexpr double wk[19] = {3., 18., 18., 18., 18., 18., 18., 36., 36., 36., 36., 36., 36., 36., 36., 36., 36., 36., 36.};

    public:

    LBM_node() = default;                         //default constructor
    LBM_node(const LBM_node& c) = default;           //copy constructor
    LBM_node(LBM_node&& c) = default;                 //move constructor
    LBM_node& operator=(const LBM_node& c) = default; //copy assignment operator
    LBM_node& operator=(LBM_node&& c) = default;      //move assignment operator 


    LBM_node( double dx, double dy, double dz, int wallflag, bool internal, int idx);

    auto get_f(){return f_;};
    auto get_m(){return m_;};
    auto get_f_eq(){return f_eq;};
    std::vector<double> get_pos(){return {dx_, dy_, dz_};};
    std::vector<double> get_velocity(){return {m_[1]/m_[0], m_[2]/m_[0], m_[3]/m_[0]};};
    auto get_wallflag(){return wallflag_;};
    auto get_internal(){return internal_;};
    int get_idx(){return idx_;};

    void set_f(std::shared_ptr<LBM_node> node_, int i){f_tmp[i] = node_->get_f()[i]; return;};
    void set_f_wall(int i, int j, double c_i){f_tmp[i] = f_[j] + (c_i / wk[i]); return;}; 
    void set_f_p(int i, std::shared_ptr<LBM_node> node_){f_tmp[i] = f_eq[i] /*node_->get_f_eq()[i]*/; return;};
    void update_f(){for(int i=1; i<velo_dim; i++)f_[i] = f_tmp[i]; return;};
    void set_velocity(std::vector<double> velo){m_[1]=velo[0]*m_[0]; m_[2]=velo[1]*m_[0]; m_[3]=velo[2]*m_[0]; return;}
    void set_internal(bool onoff, std::vector<double> velo){ internal_ =  onoff;  return;}

    void ftom(shared_ptr<lattice> L, bool flag);
    void calc_eq();
    void collision(double dt, double fext[3], bool internal);

bool is_internal(std::shared_ptr<cell> cell_ptr, std::array<double, 6> aabb);
double SignedVolume(vec3 a, vec3 b, vec3 c, const vec3 d);

void set_meq(){ m_ = meq; return;}

};



