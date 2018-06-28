#include <iostream>
#include "randomc.h"
#include "Particle.h"
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#include <random>

#include <unistd.h> //this is only for on macs!

#include "Simulation.h"
#include "GetParams.h"

void progressBar(double percent);

void only_infer_number_species(const GetParams& P);
void simulate_without_fitting(GetParams P);
void fit_to_data(const GetParams& P);

void macstart(const char * argv[]);


int main(int argc, const char * argv[]) {

	//this is some MAC code that changes the working directory to the directory in which the executable is placed, can be commented away on windows/linux
    macstart(argv);


	//read parameters from parameter file
	//the parameter file is named "config.ini" and contains the following parameters:
	//	seed = 2330
	//	emp_file_name = noPelagicEgg_LowMobility.csv
	//	n_lineages = 50000
	//	softness = 0.0
	//	sampling = 1.0
	//	speciationRate = 0.0005
	//	suboptimal = 0.0
	//	numParticles = 1000
	//	alpha = 0.5
	//	Xmean = 0.5
	//  abunSubop = 0.5
	//  fitting = 0
	//  replicates = 100

	GetParams P;
	P.readFromIni("config.ini");

    if(P.custom_mask) {
        read_mask(P.custom_mask_vector, P.mask_file_name, P.n_lineages);
    }
	
	//hardcoded parameter values:
	/*P.seed = 1;
	P.emp_file_name = "noPelagicEgg_LowMobility.csv";
	P.n_lineages = 50000;
	P.softness = 0.0;
	P.sampling = 1.0;
	P.speciationRate = 0.0005;
	P.alpha = 3.0;
	P.Xmean = 0.01;
	P.replicates = 100;
	P.fitting = 0;
	P.re = 0.5;
    P.ce = 7.5;
	P.SoftBorder = 0;*/

    set_seed(P.seed);

    if(P.fitting == 2) {
        only_infer_number_species(P);
    }

	if(P.fitting == 0)  {
        simulate_without_fitting(P);
    }

    if(P.fitting == 1) { //perform ABC fitting
        fit_to_data(P);
	}

	return 0;
}

void only_infer_number_species(const GetParams& P) {
    std::ofstream outFile("only_num_species.txt", std::ios::app);
    for(int r = 0; r < P.numParticles; ++r) {
        if(r % 1000 == 0) {
            std::cout << r << "\n";
        }
        particle candidate; //candidate parameter combination
        candidate.getFromPrior(); //get combination from the prior
        doSimulation(candidate,P, candidate.numSpecies); //use the parameter combination to run the simulation
        outFile << candidate << "\n"; outFile.flush();  //store the particle in a file
    }
    outFile.close();
    return;
}

void simulate_without_fitting(GetParams P) {
    std::cout<<"please enter name for your outputs\n";
    std::cout<<"Recomended order: no individuals-speciation-date DDMMYY\n";
    std::cout<<"file name : ";
    char sDataFileName[80] = "test";
    //std::cin >> sDataFileName;
    std::stringstream ss;  ss << sDataFileName;   //some tricks to store the output file name into a string I HATE STRINGS (in C++)
    ss  >> P.sDataFileName;

    //P.sDataFileName = std::string(sDataFileName);

    char sDataFileName1[80] = "";
    char sDataFileName2[80] = "";
    char sDataFileName3[80] = "";
    char sDataFileName4[80] = "";



    strncpy(sDataFileName1, sDataFileName,80); // Raw data
    strncpy(sDataFileName2, sDataFileName,80); // Raw Metrics data
    strncpy(sDataFileName3, sDataFileName,80); // ICD of theoretical data
    strncpy(sDataFileName4, sDataFileName,80); // outcome results
    const std::string DataName2 = std::string(sDataFileName2) + "Metric_Data.txt";
    const std::string DataName3 = std::string(sDataFileName3) + "ICD_Data.txt";
    const std::string DataName4 = std::string(sDataFileName4) + "Outcome_Data.txt";

    //const char * name1 = DataName1.c_str();
    const char * name2 = DataName2.c_str();
    const char * name3 = DataName3.c_str();
    const char * name4 = DataName4.c_str();

    //write the headers to a file!

    std::ofstream ofDataFile2(name2,std::ios_base::app);
    ofDataFile2<<"alpha,xmean,speciation,sampling,SoftBorder,replicate,speciesID,Abundance,Rel_Range"<<std::endl;
    ofDataFile2.close();
    std::ofstream ofDataFile3(name3,std::ios_base::app);
    ofDataFile3<<"alpha,xmean,speciation,sampling,SoftBorder,replicate,x,y"<<std::endl;
    ofDataFile3.close();
    std::ofstream ofDataFile4(name4,std::ios_base::app);
    ofDataFile4<<"alpha,xmean,speciation,sampling,SoftBorder,Richness_emp,Richness_theo,LS"<<std::endl;
    ofDataFile4.close();

    std::cout<<"STARTING...\n";



    theta init(P.alpha, P.Xmean, P.speciationRate, P.protractedNess, P.sampling); //you can write a loop around this one to loop through different parameter values
    particle candidate(init);

    for(int r = 0; r < P.replicates; ++r) {
        std::cout << "replicate: \t" << r << "\n";
        P.focalReplicate = r;
        int numSpecies;
        std::vector<double> result = doSimulation(candidate, P, numSpecies); //use the parameter combination to run the simulation
    }
    std::cout << "DONE!\n";
    return;
}

void fit_to_data(const GetParams& P) {
    // for stopping and starting simulations
    // it is important to always have a new seed:

    std::random_device rdev;
    unsigned int chosen_seed = rdev();

    // std::random_device uses the internal randomizer
    // to randomly generate a seed.

    set_seed(chosen_seed);
    std::ofstream seedFile("seed.txt",std::ios::app);
    seedFile << chosen_seed << "\n";
    seedFile.close();

    std::ofstream errorFile("log.txt");
    std::cout << "reading empirical data" << "\n";
    errorFile << "reading empirical data\n"; errorFile.flush();
    // read the empirical curve from the file specified in the config.ini file
    std::vector<double> EMP = readEmpiricalData(P.emp_file_name);
    //    setTargetNumSpecies(P.emp_file_name, P.emp_numSpecies);


    // maximum number of iterations of the SMC algorithm
    int maxIter = 50;
  //  int maxIter = 1; // for Travis

    double sigma = 0.05;

    std::vector<double> Thresholds; //thresholds used for the ABC algorithm
    for(int i = 0; i < (maxIter-1); ++i)    {
        Thresholds.push_back(20 * exp(-0.25*i));
    }

    int step = P.numParticles / 100; //for updating the output
    if(step < 1) step = 1;

    int iter = 0;
    int numberAccepted = 0;
    int numberProposed = 0;

    errorFile << "reading old acceptance file\n"; errorFile.flush();
    std::cout << "reading old acceptance file\n";
    std::ifstream acceptFile("acceptFile.txt");
    if(file_exists("acceptFile.txt")) {
        while(!acceptFile.eof()) {
            std::string s_iter;
            acceptFile >> s_iter;
            iter = atoi(s_iter.c_str());

            std::string s_numberAccepted;
            acceptFile >> s_numberAccepted;
            numberAccepted = atoi(s_numberAccepted.c_str());

            std::string s_numberProposed;
            acceptFile >> s_numberProposed;
            numberProposed = atoi(s_numberProposed.c_str());

            std::string s_temp;
            acceptFile >> s_temp;
            std::cout << iter << "\t" << numberAccepted << "\t" << numberProposed << "\t" << s_temp << "\n";
        }
    }
    acceptFile.close();

    // checking previous ABC output
    std::cout << "checking previous ABC output" << "\n";
    errorFile << "checking previous ABC output" << "\n"; errorFile.flush();
    for(iter = 50; iter >= 0; iter--)   {
        std::vector<particle> temp;
        readParticles(iter,temp,P.numParticles);
        if(temp.size() > 0) {
            std::cout << "Found previous results, continuing at t = " << iter << "\n";
            break;
        }
    }
    errorFile << "finished checking previous ABC output" << "\n";
    errorFile << "starting ABC\n";
    if(iter < 0) iter = 0;

    // iterating through the ABC thing
    for(; iter < maxIter; ++iter) {
        std::vector < particle > previousParticles;
        if(iter > 0) readParticles(iter-1,previousParticles, P.numParticles); //read previous output
        double maxW = 1;
        if(iter > 0) normalize(previousParticles,maxW); //normalize the weights of the accepted particles

        // set up file name
        std::stringstream ss;
        ss << iter;
        std::string time = ss.str();
        std::string f_name = "particles_t=" + time + ".txt";
        std::ofstream out_particle(f_name.c_str(), std::ios::app);

        std::string f_name2 = "curve_t=" + time + ".txt";
        std::ofstream out_curve(f_name2.c_str(), std::ios::app);

        // check if a previously unfinished run has left some already computed particles:
        std::vector< particle > temp;
        readParticles(iter,temp,P.numParticles);
        if(temp.size() > 0) {
            numberAccepted = (int)temp.size();
            std::cout << "found previous particles, starting now at iter: " << iter << " particle: " << numberAccepted << "\n";
            errorFile << "found previous particles, starting now at iter: " << iter << " particle: " << numberAccepted << "\n";

        }

        for(; numberAccepted <  P.numParticles  ; ++numberProposed) { //go collect more particles
            particle candidate; //candidate parameter combination
            if(iter == 0) candidate.getFromPrior(); //get combination from the prior
            if(iter > 0) {
                candidate.getFromPrevious(previousParticles,maxW); //draw a combination proportional to the weight from the previous iteration
                candidate.perturb(sigma); //change the values a little bit
            }

            //candidate.Params.sampling = P.sampling;

            if(candidate.withinPrior()) {
                // std::cout << candidate.Params << "\n";
                std::vector<double> result;
                if(iter > 0) {
                    result = doSimulation(candidate,P, candidate.numSpecies); //use the parameter combination to run the simulation
                    //    candidate.calcFitness(EMP,result, P.emp_numSpecies); //calculate the fitness
                    candidate.calcFitness(EMP, result);
                }

                if(candidate.accept(Thresholds[iter],iter)) { //check whether the simulation outcome is sufficiently similar
                    if(iter > 0) {
                        candidate.calculateWeight(previousParticles, sigma); //calculate the weight of the particle

                        if(isnan(candidate.weight)) {
                            std::cout << "ERROR\t" << candidate << "\n";
                        } else {
                            for(auto it = result.begin(); it != result.end(); ++it) {
                                out_curve << (*it) << "\t";
                            }
                            out_curve << "\n";
                        }
                    }

                    out_particle << candidate << "\n"; out_particle.flush();  //store the particle in a file
                    numberAccepted++;
                    errorFile << "Accepted " << numberAccepted << " particles of " << P.numParticles << "\n";
                }

                if(numberProposed % 1000 == 0) {
                    // update acceptance rate and write to file in case another simulation will have to restart after this one
                    double acceptRate = 1.0 * numberAccepted/ numberProposed;
                    std::ofstream acceptFile("acceptFile.txt",std::ios::app);
                    acceptFile << iter << "\t" << numberAccepted << "\t" << numberProposed << "\t" << acceptRate << "\n";
                    acceptFile.close();
                    if(1.0 * numberProposed / numberAccepted > 1e6 && numberProposed > 1000) {
                        // accept rate drops below 1 in a million
                        std::ofstream stopFile("stopFile.txt");
                        stopFile << "Acceptance rate below: " << numberAccepted  << " in " << numberProposed << "particles\n";
                        stopFile.close();
                        return;
                    }
                }
            }
        }

        std::cout << "\niteration " << iter << "done, acceptrate: " << "\t" << 1.0 * numberAccepted / numberProposed << "\n";
        out_particle.close();
        numberProposed = 0;
        numberAccepted = 0;
    }
    errorFile.close();
    std::cout << "Done!\n";
}

void progressBar(double percent)
{
    if(percent == 0) std::cout << "\n";

    int number = (int)(1.0 * percent / 100 * 20);

    std::cout << "\r"; //clear line
    for(int i = 0; i <= number; ++i)
    {
        std::cout << "=";
    }
    std::cout << " " << percent << "%";
    std::cout.flush();
    return;
}

void macstart(const char * argv[])  {
#ifdef __APPLE__
    {
        //   char *dirsep = strrchr(argv[0], '/');
        //   if (dirsep != NULL) *dirsep = 0;
        int changeDir = chdir(argv[0]);
        std::cout << "Changing Dir: " << changeDir << "\n";
        std::string cwd = getcwd(NULL, 0);
        std::cout << cwd << "\n";
        std::cout << "Starting simulation\n";
    }
#endif
}


