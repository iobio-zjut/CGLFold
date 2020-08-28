
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
	
	void get_parameters(map<string,Real> &parametersMap);
	
	void rotation(vector<Real> &rotationAngles);
	
	Real calculate_bais();
	
	Real score(vector<Real> &rotationAngles);
	
	Real difference(vector<Real> &Angles_1, vector<Real> &Angles_2);
	
	vector<vector<Real> > angle_learning();
	
	Size dimension(){
		return rotationAxis.size();
	};
	
	void rotation_parameters(vector<vector<Real> > &init_points, vector<pair<vector<Real>, vector<Real> > > &rotation_axis);
  
private:
	
	vector<vector<Real> > targetPoints;
	
	vector<vector<Real> > initPoints;
	
	vector<pair<vector<Real>, vector<Real> > > rotationAxis;
	
private:
	Size NP2;
	Size G2;
	Real F;
	Real CR;
	Real KTl;
	Real KT_reciprocal;
	Real Max_Disturbance;
	Real Use_best;
	Size N_Candidate;
	Real Greedy_strategy;
	
public:
	
	static Real TrialScore;
	static bool isSimilar(Real &score){
		return ( fabs(score - TrialScore) < 0.00001 );
	}
};

} //abinitio
} //protocols

#endif
