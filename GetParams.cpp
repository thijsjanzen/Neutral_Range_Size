#include "GetParams.h"
#include <fstream>
#include <iostream>

//This cpp file was adapted from
//H.Hildenbrandt 2007




GetParams::GetParams() {
	seed = 1;
	std::string emp_file_name = "noPelagicEgg_HighMobility.csv";
	n_lineages = 5000;
	softness = 0.0;
	sampling = 0.8;
	speciationRate = 0.0005;
	SoftBorder = 0;
	
	re = 0.5;
	ce = 7.5;

	Xmean = 0.3;
	alpha = 2;
	
	numParticles = 1000;
	fitting = 1;
	replicates = 100;
	
}

void GetParams::readFromIni( const char * filename ) {
	
	//locate file and tranfer text to stringstream
	std::ifstream ifs( filename ); 
	std::stringstream ss;

	if( ifs ) { //only for succesfully created ifstream: otherwise null-pointer?
		ss << ifs.rdbuf(); //config.ini content is transferred to stringstream: easier to search?
	}
	else {
		throw "Can't locate file";
		//std::cout << "Can't locate file" << std::endl;
	}
	while( ss.good() ) {	
				readNameValuePair( ss,  "seed", seed);
				readNameValuePair( ss,  "emp_file_name", emp_file_name);
				readNameValuePair( ss,  "n_lineages", n_lineages);
                readNameValuePair( ss,  "numParticles",numParticles);
                readNameValuePair( ss, "replicates",replicates);
                readNameValuePair( ss, "fitting", fitting);
				readNameValuePair( ss,  "speciationRate", speciationRate);
                readNameValuePair( ss,  "protractedNess", protractedNess);
				readNameValuePair( ss,  "alpha",alpha);
				readNameValuePair( ss,  "Xmean",Xmean);
                readNameValuePair( ss,  "sampling", sampling);
	}
}

template <typename T>
void GetParams::readNameValuePair( std::stringstream& ss, std::string iniName, T& value ) {
	std::string name;
	char sign;
	ss >> name; //>> copies ss content to string until white space is encountered
	if( name != iniName )
		throw "expect parameter";
	ss >> sign; //copies ss content to character 
	if( sign != '=' )
		throw "text format of ini file is not compatible";
    ss >> value;
	std::cout << iniName << ": " << value << std::endl;
//	std::cerr << iniName << ": " << value << std::endl;
}


