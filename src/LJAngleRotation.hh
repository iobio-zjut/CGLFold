
#ifndef INCLUDED_protocols_abinitio_LJAngleRotation_hh
#define INCLUDED_protocols_abinitio_LJAngleRotation_hh

#include<vector>
#include<map>
#include<string>
#include <core/types.hh>
#include <math.h>

using namespace std;
using namespace core;


namespace protocols {
namespace abinitio {
  
class LJAngleRotation {
public:
	LJAngleRotation() {};
	
	void init(vector<vector<Real> > &init_points, vector<vector<Real> > &rotation_axis);
	
	void get_parameters(map<string,Real> &parametersMap);
	
	void rotation(vector<Real> &rotationAngles);
	
	Real calculate_bais();
	
	Real score(vector<Real> &rotationAngles);
	
	Real difference(vector<Real> &Angles_1, vector<Real> &Angles_2);
	
	vector<vector<Real> > angle_learning();
	
	Size dimension(){
		return rotationAxis.size();
	};
	
	void rotation_parameters(vector<vector<Real> > &init_points, vector<vector<Real> > &rotation_axis);
  
private:
	///@note 旋转点的坐标集
	vector<vector<Real> > targetPoints;
	///@note 旋转点的初始坐标集
	vector<vector<Real> > initPoints;
	///@note 旋转轴集
	vector<vector<Real> > rotationAxis;
	
private:
	Size NP2;
	Size G2;
	Real F;
	Real CR;
	Real KTl;
	Real KT_reciprocal;
//	Real Driver_Angle;
	Real Max_Disturbance;
	Real Use_best;
	Size N_Candidate;
	Real Greedy_strategy;
//	Size Convergence_G;
	
public:
	///@brief 相似性判断函数，作为泛型算法 find_if 的谓词
	static Real TrialScore;
	static bool isSimilar(Real &score){
		return ( fabs(score - TrialScore) < 0.00001 );
	}
};

} //abinitio
} //protocols

#endif
