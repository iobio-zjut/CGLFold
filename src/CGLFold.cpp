#include<cstdlib>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include <algorithm>
#include <math.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<iomanip>
#include<unistd.h>

using namespace std;


struct ContactPair
{
	int r1;
	int r2;
	double confidence;
	
	ContactPair(){r1=0; r2=0; confidence=0;}
};
//sequence
string fasta;
//sequence route
string fasta_route;
//seqence length
int fasta_length;
//contact map route
string cmap_route;
//contact map used to build the Cscore model
vector<ContactPair> cmap;

//native_pdb route
string native_route;
//fragment library route
string frag_lib_3_route;
string frag_lib_9_route; 

//check whether fastaFile is input
bool fFile;
//check whether cmapFile is input
bool cFile;
//check whether frag_library are input
bool f3File;
bool f9File;
//check whether ntive structure is input
bool nFile;


void getParameter(int argc, char ** argv, int & i);
void extractParameters(int argc, char ** argv);
bool is_amino(char &i);
// read the fasta file
void getFasta();
//filter and extract original contact
void getCmap();

void apply();

void cluster();

int main(int argc, char ** argv){
	cout << endl;
	cout << "###########################################################################" << endl;
	cout << "#                                CGLFold                                  #" << endl;
	cout << "#                              Version: 1.0                               #" << endl;
	cout << "#          A contact-assisted de novo protein structure prediction        #" << endl;
	cout << "#    using global exploration and loop perturbation sampling algorithm    #" << endl;
	cout << "#                         Copyright (C) Guijun Zhang                      #" << endl;
	cout << "#                     College of Information engineering                  #" << endl;
	cout << "#            Zhejiang university of technology, Hangzhou, China           #" << endl;
	cout << "###########################################################################" << endl;
	
	
	cout << access("output_files", 00);
	
	if ( !access("output_files", 00) )
		system("rm -r output_files");
	
//	system("mkdir -p output_files/PDB");
	system("mkdir -p output_files/SPICKER_data");
	
	extractParameters(argc, argv);
	getFasta();
	getCmap();
	
	apply();
	
	cluster();
	
	cout << endl;
	cout << "###########################################################################" << endl;
	cout << "#                                  END                                    #" << endl;
	cout << "###########################################################################" << endl;
  
	return 0;
}

void getFasta(){
    ifstream fasta_file( fasta_route.c_str() );
    
    if ( fasta_file.is_open() ){
	string line;
	while( getline(fasta_file, line) ){
	    if ( line[0] != '>' )
		fasta += line;
	}
	fasta_file.close();
	
	for (int i = 0; i < fasta.size(); ++i){
	    if ( fasta[i] == ' ' ){
		    fasta.erase(i,1);
		    --i;
	    }
	    else{
		if ( !is_amino( fasta[i] ) ){
		    cout << "Error! Invalid amino acid " << fasta[i] << endl;
		    exit(0);
		}
	    }
	}
	fasta_length = fasta.size();
	cout << "#  " << fasta << endl;
	cout << "#  fasta length: " << fasta_length << endl;
    }
    else{
	cout << "Error! fasta file can not open " << fasta_route << endl;
	exit(0);
    }
}

bool is_amino(char &r){
	if (r == 'A')
            return true;
        else if (r == 'C')
            return true;
        else if (r == 'D')
            return true;
        else if (r == 'E')
            return true;
        else if (r == 'F')
            return true;
        else if (r == 'G')
            return true;
        else if (r == 'H')
            return true;
        else if (r== 'I')
            return true;
        else if (r == 'K')
            return true;
        else if (r == 'L')
            return true;
        else if (r == 'M')
            return true;
        else if (r == 'N')
            return true;
        else if (r == 'P')
            return true;
        else if (r == 'Q')
            return true;
        else if (r == 'R')
            return true;
        else if (r == 'S')
            return true;
        else if (r == 'T')
            return true;
        else if (r == 'V')
            return true;
        else if (r == 'W')
            return true;
        else if (r == 'Y')
            return true;
        else 
	    return false;
}

void getCmap(){
    ifstream contact_file( cmap_route.c_str() );
    
    if( contact_file.is_open() ){
	string line;
	//extract top L/2 filtered contact
	int num_contact( fasta_length / 2 + fasta_length % 2 );
//	cout << "num_contact: " << num_contact << endl;
//	int count_contact( 0 );
	while( getline(contact_file, line) ){
	    istringstream line_data( line );
	    ContactPair contact;
	    
//	    int zero, eight;
//	    line_data >> contact.r1 >> contact.r2 >> zero >> eight >> contact.confidence;
	    
	    line_data >> contact.r1 >> contact.r2 >> contact.confidence;
	    
//	    if ( zero == 0 && eight == 8 && 
	    if (contact.r1 > 0 && contact.r1 <= fasta_length &&
		contact.r2 > 0 && contact.r1 <= fasta_length && 
		contact.confidence > 0 && contact.confidence <= 1 )
	    {
		bool save = true;
		for ( int j = 0; j < cmap.size(); ++j ){
		    if ( sqrt( pow((contact.r1 - cmap[j].r1), 2) + pow((contact.r2 - cmap[j].r2), 2) ) <= 2 ){
			save = false;
			break;
		    }
		}
		if ( save )
			cmap.push_back(contact);
	    }
	    if ( cmap.size() >= num_contact )
		break;
	}
	cout << "#  extracted contacts number: " << cmap.size() << endl;
    }
    else{
	cout << "Error! contact map file can not open " << cmap_route << endl;
	exit(0);
    }
    
    //wright to file
    ofstream filtered_cmap("output_files/filtered_cmap.txt");
    for (int i = 0; i < cmap.size(); ++i)
	filtered_cmap << right << setw(5) << cmap[i].r1
	    << right << setw(5) << cmap[i].r2
	    << right << setw(12) << cmap[i].confidence
	    << endl;
    filtered_cmap.close();
}

void extractParameters(int argc, char ** argv) {
    int i = 1;
//    cout << argc << endl;
    while (i < argc)
        getParameter(argc, argv, i);
    if (!fFile) {
        cout << endl;
        cout << "Error! fasta File must be provided" << endl << endl;
        exit(0);
    }
    if (!cFile) {
        cout << endl;
        cout << "Error! contact map File must be provided" << endl << endl;
        exit(0);
    }
    if (!f3File) {
        cout << endl;
        cout << "Error! 3-mer fragment library File must be provided" << endl << endl;
        exit(0);
    }
    if (!f9File) {
        cout << endl;
        cout << "Error! 9-mer fragment library File must be provided" << endl << endl;
        exit(0);
    }
//     if (!nFile) {
//         cout << endl;
//         cout << "Error! native structure File must be provided" << endl << endl;
//         exit(0);
//     }
}

void getParameter(int argc, char ** argv, int & i){
	string flag( argv[i] );
	if ( flag == "-f" ){
		if (argc < i + 2){
			cout << endl;
			cout << "Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -c 		cmap   			: contact map file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -n 		native 			: native structure for comparison (optional)" << endl;
			exit(0);
		}
		fasta_route = argv[++i];
	//	cout << fasta_route << endl;
		fFile = true;
	}
	else if( flag == "-c" ){
		if (argc < i + 2){
			cout << endl;
			cout << "Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -c 		cmap   			: contact map file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -n 		native 			: native structure for comparison (optional)" << endl;
			exit(0);
		}
		cmap_route = argv[++i];
	//	cout << cmap_route << endl;
		cFile = true;
	}
	else if( flag == "-frag3" ){
		if (argc < i + 2){
			cout << endl;
			cout << "Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -c 		cmap   			: contact map file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -n 		native 			: native structure for comparison (optional)" << endl;
			exit(0);
		}
		frag_lib_3_route = argv[++i];
	//	cout << frag_lib_3_route << endl;
		f3File = true;
	}
	else if( flag == "-frag9" ){
		if (argc < i + 2){
			cout << endl;
			cout << "Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -c 		cmap   			: contact map file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -n 		native 			: native structure for comparison (optional)" << endl;
			exit(0);
		}
		frag_lib_9_route = argv[++i];
	//	cout << frag_lib_9_route << endl;
		f9File = true;
	}
	else if( flag == "-n" ){
		if (argc < i + 2){
			cout << endl;
			cout << "Error! No fasta file provided" << endl << endl;
			cout << "   -f 		fasta  			: fasta file  " << endl;
			cout << "   -c 		cmap   			: contact map file " << endl;
			cout << "   -frag3	3mer_fragment_library  	: fragment library with fragment length 3" << endl;
			cout << "   -frag9	9mer_fragment_library  	: fragment library with fragment length 9" << endl;
			cout << "   -n 		native 			: native structure for comparison (optional)" << endl;
			exit(0);
		}
		native_route = argv[++i];
	//	cout << native_route << endl;
		nFile = true;
	}/*
	else{
		cout << "error!!!!!!!!!!" << endl;
		exit(0);
	}*/
	++i;
}

void apply(){
	stringstream command;
	command << "~/rosetta_src_2018.33.60351_bundle/main/source/bin/AbinitioRelax.linuxgccrelease "
		<< "-in:file:fasta " << fasta_route << " "
	//	<< "-in:file:native " << native_route << " "
		<< "-in:file:frag3 " << frag_lib_3_route << " "
		<< "-in:file:frag9 " << frag_lib_9_route << " "
		<< "-abinitio:increase_cycles " << 1 << " "
		<< "-nstruct " << 1 << " "
		<< "-out:pdb";
	string run_rosetta;
	getline(command, run_rosetta);
	cout << run_rosetta << endl;
		
	system(run_rosetta.c_str());
}

void cluster(){
	ofstream output("output_files/SPICKER_data/cluster.sh");
	output << "cd output_files/SPICKER_data/" << endl;
	output << "./../../../bin/spicker" << endl;
	output << "./../../../bin/pulchra -p combo1.pdb" << endl;
	output << "mv pul_combo1.pdb ../model_1.pdb" << endl;
	output << "./../../../bin/pulchra -p closc1.pdb" << endl;
	output << "mv pul_closc1.pdb ../model_2.pdb" << endl;
	output << "./../../../bin/pulchra -p combo2.pdb" << endl;
	output << "mv pul_combo2.pdb ../model_5.pdb" << endl;
	output << "cd ../../" << endl;
	output << "rm S_00000001.pdb score.fsc default.out" << endl;
	output << "rm -r output_files/SPICKER_data" << endl;
//	output << "rm output_files/SPICKER_data/cluster.sh" << endl;
	output.close();
	system("chmod 777 output_files/SPICKER_data/cluster.sh");
	system("./output_files/SPICKER_data/cluster.sh");
}

