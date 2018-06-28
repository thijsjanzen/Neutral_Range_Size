#ifndef GETPARAMS_H_INCLUDED
#define GETPARAMS_H_INCLUDED

#include <sstream>
#include <string>
#include <vector>


class GetParams {
public:
	GetParams();
	void readFromIni( const char * filename );
	
	template <typename T>
	void readNameValuePair( std::stringstream& ss, std::string iniName, T& value );
	
	int seed;
	std::string emp_file_name;
    int emp_numSpecies;

	int n_lineages;
	double softness;
	double sampling;
	double speciationRate;
    double protractedNess;
	int SoftBorder;

	int numParticles;

	double ce;
	double re;
	
	int singleRun;
	double Xmean;
	double alpha;

	int fitting;
	int replicates;
	std::string sDataFileName;
	int focalReplicate;

    bool custom_mask;
    std::string mask_file_name;
    std::vector<bool> custom_mask_vector;
};

#endif //GETPARAMS_H_INCLUDED
