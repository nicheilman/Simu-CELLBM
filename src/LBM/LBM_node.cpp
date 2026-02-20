
#include "lattice.hpp"
#include "IB.hpp"
#include "../math_modules/vec3.hpp"

LBM_node::LBM_node(const double dx, const double dy, const double dz, int wall_flag, bool internal, int idx){

f_ = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
m_ = {1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
f_tmp = f_;
lambda = {0, 0, 0, 0, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
dx_ = dx;
dy_ = dy;
dz_ = dz;
wallflag_ = wall_flag;
internal_ = internal;
idx_ = idx;

}

//--------------------------------------------------------------------------------//

void LBM_node::ftom(shared_ptr<lattice> L, bool flag){
if(!flag){
std::fill(std::begin(m_), std::end(m_), 0);

//#pragma omp parallel for schedule(static)
for(int i=0; i<19; i++){
	for(int j=0; j<19; j++){
		m_[i] += L->evector[i][j] * f_[j];
	}
    }
}

else{
std::fill(std::begin(f_), std::end(f_), 0); 
//#pragma omp parallel for
for(int i=0; i<19; i++){
        for(int j=0; j<19; j++){
                f_[i] +=   (L->inv_evector[i][j] * m_[j]) / 144;
if(f_[i] < 0.0) f_[i] = 0.0;
        }
    }
std::fill(std::begin(f_eq), std::end(f_eq), 0);
//#pragma omp parallel for
for(int i=0; i<19; i++){
        for(int j=0; j<19; j++){
                f_eq[i] +=   (L->inv_evector[i][j] * meq[j]) / 144;
        }
    }
}

return;
}
//--------------------------------------------------------------------//

void LBM_node::calc_eq(){

            m2 = m_[1] * m_[1] + m_[2] * m_[2] + m_[3] * m_[3];
            meq[0] = m_[0];
            meq[1] = m_[1];
            meq[2] = m_[2];
            meq[3] = m_[3];
            meq[4] = m2 / meq[0];
            meq[5] = (3 * m_[1] * m_[1] - m2) / meq[0];
            meq[6] = (m_[2] * m_[2] - m_[3] * m_[3]) / meq[0];
            meq[7] = m_[1] * m_[2] / meq[0];
            meq[8] = m_[3] * m_[2] / meq[0];
            meq[9] = m_[3] * m_[1] / meq[0];
            for (int q = 10; q < velo_dim; ++q)
                meq[q] = 0.0;

return;
}
//--------------------------------------------------------------------//

void LBM_node::collision(double dt, double fext[3], bool internal){
//#pragma omp parallel for
for(int i=4; i<velo_dim; i++)
	m_[i] += -1.0 * (m_[i] - meq[i]) * lambda[i];

//if(internal) for(int i=0; i<3; i++) m_[i+1] = 0.0; 

double f_extx = fext[0];
double f_exty = fext[1];
double f_extz = fext[2];

//if(dz_ == -4.0e-5) meq[3] = 1.0e-15;
//if(dz_ == -3.0e-5) meq[2] =  1.0e-11;

if(!internal){
      m_[1] += f_extx;
      m_[2] += f_exty;
      m_[3] += f_extz;
      m_[4] += 2.0*(m_[1]*f_extx+m_[2]*f_exty+m_[3]*f_extz)/meq[0];
      m_[5] += (4.0*m_[1]*f_extx-2.0*(m_[2]*f_exty+m_[3]*f_extz))/meq[0];
      m_[6] += 2.0*(m_[2]*f_exty-m_[3]*f_extz)/meq[0];
      m_[7] += (m_[1]*f_exty+m_[2]*f_extx)/meq[0];
      m_[8] += (m_[3]*f_exty+m_[2]*f_extz)/meq[0];
      m_[9] += (m_[3]*f_extx+m_[1]*f_extz)/meq[0];
}

if(internal) for(int i=0; i<3; i++) m_[i+1] = 0.0;

return;

}


bool LBM_node::is_internal(std::shared_ptr<cell> cell_ptr, array<double, 6> aabb){

vec3 test_point(dx_, dy_, dz_);

aabb = cell_ptr->get_aabb();
vec3 arrow = ((test_point - cell_ptr->compute_centroid() ).normalize()); //*(4.0*max(max(aabb[3]-aabb[0], aabb[4]-aabb[1]), aabb[5]-aabb[2]));
vec3 final_point = test_point - arrow*(4.0*max(max(aabb[3]-aabb[0], aabb[4]-aabb[1]), aabb[5]-aabb[2]));

auto node_lst_ = cell_ptr->get_node_lst();
int count = 0;
int face_counter = 0;
vec3 a, b, c;

for(auto& face : cell_ptr->get_face_lst() ){
	if(!face.is_used()) continue;
	a = node_lst_[face.n1_id()].pos();
	b = node_lst_[face.n2_id()].pos();
	c = node_lst_[face.n3_id()].pos();

//TEST
if(!octant_check(arrow, a, test_point) && !octant_check(arrow, b, test_point) && !octant_check(arrow, c, test_point)) continue;

//face_counter += 1;

	     if( ((SignedVolume(a, b, c, test_point)>0) != (SignedVolume(a, b, c, final_point)>0)) ) continue;
        else if( ((SignedVolume(test_point, final_point, a, b)>0) == (SignedVolume(test_point, final_point,b, c)>0) ) &&
                 ((SignedVolume(test_point, final_point, a, b)>0) == (SignedVolume(test_point, final_point,c, a)>0) ) ) {count += 1;}
    
	}

//printf("Internal used %d faces\n", face_counter);
if(count % 2 == 1) return 1;
else return 0;

//--------------------------------------------

}

double LBM_node::SignedVolume(vec3 a, vec3 b, vec3 c, const vec3 d){
	return /*1.0/6.0 */ (a-d).dot((b-d).cross(c-d));}





