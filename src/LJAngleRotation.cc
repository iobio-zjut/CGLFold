
#include <protocols/abinitio/LJAngleRotation.hh>
#include <numeric/random/random.hh>

#include <map>

#include<iostream>
#include<fstream>
#include<sstream>
#include<utility>
#include<iomanip>
#include<math.h>
#include<algorithm>
#include<stdlib.h>

using namespace std;


namespace protocols {
namespace abinitio {
  
Real LJAngleRotation::TrialScore = 0.0;

void 
LJAngleRotation::init(vector< vector< Real > >& init_points, vector< vector< Real > >& rotation_axis){
	initPoints = init_points;
	rotationAxis = rotation_axis;
//	read_parameters();
}

void 
LJAngleRotation::rotation_parameters(vector< vector< Real > >& init_points, vector< vector< Real > >& rotation_axis){
	initPoints = init_points;
	rotationAxis = rotation_axis;
}


void 
LJAngleRotation::get_parameters(map<string,Real>& parametersMap){
	
	if (parametersMap.find("NP2") != parametersMap.end()){
	    NP2 = parametersMap["NP2"];
	}
	else{
	    cout << "====================================parameter 'NP2' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("G2") != parametersMap.end())
	    G2 = parametersMap["G2"];
	else{
	    cout << "====================================parameter 'G2' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("CR") != parametersMap.end())
	    CR = parametersMap["CR"];
	else{
	    cout << "====================================parameter 'CR' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("F") != parametersMap.end())
	    F = parametersMap["F"];
	else{
	    cout << "====================================parameter 'F' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("KTl") != parametersMap.end()){
	    KTl = parametersMap["KTl"];
	    KT_reciprocal = (Real)1/KTl;
	    cout << "Loop Search 1/KT = " << KT_reciprocal << endl;
	}
	else{
	    cout << "====================================parameter 'KT' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	Max_Disturbance = 5;
	Use_best = 1;
	N_Candidate = 10;
	Greedy_strategy = 1;
/*	
	if (parametersMap.find("Max_Disturbance") != parametersMap.end())
	    Max_Disturbance = parametersMap["Max_Disturbance"];
	else{
	    cout << "====================================parameter 'Max_Disturbance' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("Use_best") != parametersMap.end())
	    Use_best = parametersMap["Use_best"];
	else{
	    cout << "====================================parameter 'Use_best' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("N_Candidate") != parametersMap.end())
	    N_Candidate = parametersMap["N_Candidate"];
	else{
	    cout << "====================================parameter 'N_Candidate' is not found in parameters!!!" << endl;
	    exit(0);
	}
	
	if (parametersMap.find("Greedy_strategy") != parametersMap.end())
	    Greedy_strategy = parametersMap["Greedy_strategy"];
	else{
	    cout << "====================================parameter 'Greedy_strategy' is not found in parameters!!!" << endl;
	    exit(0);
	}
*/
}

vector<vector<Real> > 
LJAngleRotation::angle_learning(){
///===========================================================================================================================================================================	

	Size Angle_Dimensions( dimension() );
	vector<vector<Real> > Population_Angles;
	vector<Real> Population_Energys;
	
	for (Size p = 1; p <= NP2; ++p){
		vector<Real> Angles;
		for (Size d = 0; d < Angle_Dimensions; ++d){
			Angles.push_back( numeric::random::rg().uniform() * 2 *Max_Disturbance - Max_Disturbance );
		}
		Population_Angles.push_back( Angles );
		Population_Energys.push_back( score( Angles ) );
	}
	
	Size accept_num( 0 );
	for (Size g = 1; g <= G2; ++g){
	//	cout << "=============================================================== g = " << g << endl;
		for (Size p = 0; p < Population_Angles.size(); ++p){
			Size base(0);
			if ( Use_best && numeric::random::rg().uniform() <= 0.5 )
				base = distance(Population_Energys.begin(), min_element( Population_Energys.begin(), Population_Energys.end() ) );
			else
				base = numeric::random::rg().random_range(0, Population_Angles.size() - 1);
			
			Size rand1( numeric::random::rg().random_range(0, Population_Angles.size() - 1) );
			Size rand2( numeric::random::rg().random_range(0, Population_Angles.size() - 1) );
			while (rand1 == base)
				rand1 = numeric::random::rg().random_range(0, Population_Angles.size() - 1);
			while (rand2 == base || rand2 == rand1)
				rand2 = numeric::random::rg().random_range(0, Population_Angles.size() - 1);
			///@note Mutation
			vector<Real> Mutant_Angles;
			for (Size d = 0; d < Angle_Dimensions; ++d){
				Mutant_Angles.push_back( Population_Angles[base][d] + F * ( Population_Angles[rand1][d] - Population_Angles[rand2][d] ) );
			}
			///@note Crossover
			vector<Real> Cross_Angles;
			Size rand_d( numeric::random::rg().random_range(0, Angle_Dimensions - 1) );
			for (Size d = 0; d < Angle_Dimensions; ++d){
				if ( numeric::random::rg().uniform() <= CR || d == rand_d )
					Cross_Angles.push_back( Mutant_Angles[d] );
				else
					Cross_Angles.push_back( Population_Angles[p][d] );
			}
			///@note Selection
			Real targetScore( Population_Energys[p] );
			Real trialScore( score( Cross_Angles ) );
			
			bool success( false );
			if ( trialScore <= targetScore )
				success = true;
			else{
				if ( exp( -(trialScore - targetScore) * KT_reciprocal ) < numeric::random::rg().uniform() )
					success = false;
				else
					success = true;
			}
			if ( success ){
			//	++accept_num;
				TrialScore = trialScore;
				if ( find_if(Population_Energys.begin(), Population_Energys.end(), isSimilar) == Population_Energys.end() ){
					Population_Angles[p] = Cross_Angles;
					Population_Energys[p] = trialScore;
				}
				else{
					while ( 1 ){
						vector<Real> Angles;
						for (Size d = 0; d < Angle_Dimensions; ++d){
							Angles.push_back( numeric::random::rg().uniform() * 2 *Max_Disturbance - Max_Disturbance );
						}
						TrialScore = score( Angles );
						if ( find_if(Population_Energys.begin(), Population_Energys.end(), isSimilar) == Population_Energys.end() ){
							Population_Angles[p] = Angles;
							Population_Energys[p] = TrialScore;
							
							break;
						}
					}
					  
				}
			}
			else{
				if ( Greedy_strategy ){
					
					Size biggest_index( distance(Population_Energys.begin(), max_element( Population_Energys.begin(), Population_Energys.end() ) ) );
					if ( trialScore < Population_Energys[biggest_index] ){
					//	++accept_num;
						TrialScore = trialScore;
						if ( find_if(Population_Energys.begin(), Population_Energys.end(), isSimilar) == Population_Energys.end() ){
							Population_Angles[biggest_index] = Cross_Angles;
							Population_Energys[biggest_index] = trialScore;
						}
						else{
							while ( 1 ){
								vector<Real> Angles;
								for (Size d = 0; d < Angle_Dimensions; ++d){
									Angles.push_back( numeric::random::rg().uniform() * 2 *Max_Disturbance - Max_Disturbance );
								}
								TrialScore = score( Angles );
								if ( find_if(Population_Energys.begin(), Population_Energys.end(), isSimilar) == Population_Energys.end() ){
									Population_Angles[biggest_index] = Angles;
									Population_Energys[biggest_index] = TrialScore;
									
									break;
								}
							} 
						}
					}
				}
			}
			
		}
	}
//	cout << accept_num << endl;
	map<Real, vector<Real> > Bais_Angle_Population;

	for (Size p = 0; p < Population_Angles.size(); ++p){
		Bais_Angle_Population.insert( make_pair( Population_Energys[p], Population_Angles[p] ) );
	}
	
	vector<vector<Real> > Candidate_Disturbance_Angle;
	Size number = 1;
	for ( map<Real, vector<Real> >::iterator iter = Bais_Angle_Population.begin(); number <= N_Candidate; ++iter, ++number ){
		Candidate_Disturbance_Angle.push_back( iter->second );
		
	//	cout << right << setw(3) << number
	//	    << right << setw(10) << setprecision(3) << fixed << iter->first;
	//	for (Size d = 0; d < Angle_Dimensions; ++d)
	//		cout << right << setw(8) << setprecision(3) << fixed << iter->second[d];
	//	cout << endl;
	}
	
	return Candidate_Disturbance_Angle;
}

Real 
LJAngleRotation::difference(vector< Real >& Angles_1, vector< Real >& Angles_2){
	Real dist( 0 );
	for (Size d = 0; d < dimension(); ++d)
		dist += pow( Angles_1[d] - Angles_2[d], 2 );
	dist = sqrt( dist );
	
	return dist;
}

Real 
LJAngleRotation::score(std::vector<Real>& rotationAngles){
	rotation(rotationAngles);
	return calculate_bais();
}

void 
LJAngleRotation::rotation(std::vector<Real>& rotationAngles){
	targetPoints = initPoints;
	
	///@brief Matrix of rotation
	vector<vector<Real> > Matrix_rotation(3, vector<Real>(3, 0));
	for (Size i = 0; i < rotationAxis.size(); ++i){
		///@note axis
		Real x( rotationAxis[i][0] );
		Real y( rotationAxis[i][1] );
		Real z( rotationAxis[i][2] );
		///@note angle
		Real angle( rotationAngles[i] * 3.14159 / 180 );
		///@note matrix
		vector<vector<Real> > matrix(3, vector<Real>(3, 0));
		matrix[0][0] = x*x*(1 - cos(angle)) + cos(angle);
		matrix[0][1] = x*y*(1 - cos(angle)) - z*sin(angle);
		matrix[0][2] = x*z*(1 - cos(angle)) + y*sin(angle);
		matrix[1][0] = x*y*(1 - cos(angle)) + z*sin(angle);
		matrix[1][1] = y*y*(1 - cos(angle)) + cos(angle);
		matrix[1][2] = y*z*(1 - cos(angle)) - x*sin(angle);
		matrix[2][0] = x*z*(1 - cos(angle)) - y*sin(angle);
		matrix[2][1] = y*z*(1 - cos(angle)) + x*sin(angle);
		matrix[2][2] = z*z*(1 - cos(angle)) + cos(angle);
		
		if ( i == 0 )
			Matrix_rotation = matrix;
		else{
			vector<vector<Real> > new_matrix(3, vector<Real>(3, 0));
			new_matrix[0][0] = Matrix_rotation[0][0] * matrix[0][0] + Matrix_rotation[0][1] * matrix[1][0] + Matrix_rotation[0][2] * matrix[2][0];
			new_matrix[0][1] = Matrix_rotation[0][0] * matrix[0][1] + Matrix_rotation[0][1] * matrix[1][1] + Matrix_rotation[0][2] * matrix[2][1];
			new_matrix[0][2] = Matrix_rotation[0][0] * matrix[0][2] + Matrix_rotation[0][1] * matrix[1][2] + Matrix_rotation[0][2] * matrix[2][2];
			
			new_matrix[1][0] = Matrix_rotation[1][0] * matrix[0][0] + Matrix_rotation[1][1] * matrix[1][0] + Matrix_rotation[1][2] * matrix[2][0];
			new_matrix[1][1] = Matrix_rotation[1][0] * matrix[0][1] + Matrix_rotation[1][1] * matrix[1][1] + Matrix_rotation[1][2] * matrix[2][1];
			new_matrix[1][2] = Matrix_rotation[1][0] * matrix[0][2] + Matrix_rotation[1][1] * matrix[1][2] + Matrix_rotation[1][2] * matrix[2][2];
			
			new_matrix[2][0] = Matrix_rotation[2][0] * matrix[0][0] + Matrix_rotation[2][1] * matrix[1][0] + Matrix_rotation[2][2] * matrix[2][0];
			new_matrix[2][1] = Matrix_rotation[2][0] * matrix[0][1] + Matrix_rotation[2][1] * matrix[1][1] + Matrix_rotation[2][2] * matrix[2][1];
			new_matrix[2][2] = Matrix_rotation[2][0] * matrix[0][2] + Matrix_rotation[2][1] * matrix[1][2] + Matrix_rotation[2][2] * matrix[2][2];
			
			Matrix_rotation = new_matrix;
		}
	}
	///@brief Point of rotation
	for (Size j = 0; j < targetPoints.size(); ++j){
		Real x_old( targetPoints[j][0] );
		Real y_old( targetPoints[j][1] );
		Real z_old( targetPoints[j][2] );
		
		///@note rotation
		Real x_new( Matrix_rotation[0][0] * x_old + Matrix_rotation[0][1] * y_old + Matrix_rotation[0][2] * z_old );
		Real y_new( Matrix_rotation[1][0] * x_old + Matrix_rotation[1][1] * y_old + Matrix_rotation[1][2] * z_old );
		Real z_new( Matrix_rotation[2][0] * x_old + Matrix_rotation[2][1] * y_old + Matrix_rotation[2][2] * z_old );
		
		targetPoints[j][0] = x_new;
		targetPoints[j][1] = y_new;
		targetPoints[j][2] = z_new;
	}
}

Real 
LJAngleRotation::calculate_bais(){
	Real total_bias( 0 );
	for (Size k = 0; k < initPoints.size(); ++k){
		Real bias( sqrt( pow( (initPoints[k][0] - targetPoints[k][0]), 2 ) + pow( (initPoints[k][1] - targetPoints[k][1]), 2 ) + pow( (initPoints[k][2] - targetPoints[k][2]), 2 ) ) );
		total_bias += sqrt(k + 1) * bias;
	}
	return total_bias;
}


} //abinitio
} //protocols
