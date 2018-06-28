//
//  Adriana_simulation.cpp
//  SMC_Adriana
//
//  Created by Thijs Janzen on 03/12/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#include "Simulation.h"
#include <vector>
#include <string>
#include <string.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cmath>
#include <vector>
#include <time.h>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <boost/tokenizer.hpp>
#include <numeric>

//**********************************************************************************************************************************
//														Simulation function
//***************************************************************************************************************************


std::vector<double> doSimulation(particle candidate,
                                 const GetParams& P,
                                 int& numSpecies) {
	int species = 1;
	int n_lineages = P.n_lineages;
    double sampling = candidate.Params.sampling; // percentage of individuals sampled

	std::vector<int> position;
	std::vector<int> descendant;
	std::vector<int> index(n_lineages);
	std::vector<int> result(n_lineages);
	std::vector<bool> speciation(n_lineages);
	std::vector<bool> mask; // mask where individuals can live
	int nExtraZeros = 0;

	//initialize vectors
	//create the mask, using exp function for softborders
	
	//--------------------20-04-2016---------------------
	//      add sampling
	//------------------------------------------

    initialize_vectors(P.mask_file_name, index, result, mask,   position,
                       descendant, sampling, nExtraZeros,
                       n_lineages, P.custom_mask);

	//--------------------------------------------------------------------
	//   end editing 20-04-2016
	//-------------------------------------------------------------------

    simulate_model(candidate,
                   position, descendant, result, index, speciation,
                   species, n_lineages, P.custom_mask, mask);

	/****************************************************************************************
	 PROCESSING THE DATA
	 *****************************************************************************************/
	
	//speciate the last lineage!
    if(descendant.size() > 0) {
        if(descendant[0] > (int)result.size() || descendant[0] > (int)speciation.size()) {
            std::cout << "\n\nERROR ERROR ERROR\n";
            std::cout << "descendant[0] outside result.size or speciation.size\n";
            std::cout << "ERROR ERROR ERROR\n\n";
            std::vector<double> curve = calcCurve(result, n_lineages, numSpecies);
            return curve;
        } else {
            result[descendant[0]] = species;
            speciation[descendant[0]]=true;
        }
    } else {
        std::cout << "\n\nERROR ERROR ERROR\n";
        std::cout << "descendant[0] outside descendant.size\n";
        std::cout << "ERROR ERROR ERROR\n\n";
        std::vector<double> curve = calcCurve(result,n_lineages, numSpecies);
        return curve;
    }

	// addition 20-04-2016
	int nZeros = ((int)speciation.size() - (species + nExtraZeros) ); //CHECK 040714
	while(nZeros > 0){
		for(int i = 0; i != (int)result.size(); ++i){
			if(result[i] != -2) {
				if(speciation[i] == 0) {
					if((speciation[i] == 0) &&
                       (speciation[result[i]] == 1)) {
						nZeros --;
					}
					speciation[i]=speciation[result[i]];
					result[i]=result[result[i]];
				}// end if statement
			}
		}//end for loop
	} //end while loop
 
	std::vector<double> curve = calcCurve(result,n_lineages, numSpecies);
	
	if(P.fitting == 0) {
        write_to_file(P,  result, n_lineages,  curve, species);
    }

	return curve;
}

void simulate_model(particle candidate,
                    std::vector<int>& position,
                    std::vector<int>& descendant,
                    std::vector<int>& result,
                    std::vector<int>& index,
                    std::vector<bool>& speciation,
                    int& species,
                    const int& n_lineages,
                    bool custom_mask,
                    const std::vector<bool>& mask) {

    int remain_lineages = (int)position.size(); // will follow only sampled lineagues
    std::vector<double> recordYpos;

    double GenerationTime = 0;
    double gtminus;
    int steps = 0;
    double D; // distance of the random walk
    double Xmin; //Xmin of the pareto distribution, it is defined below in terms of Xmean

     double speciationRate = candidate.Params.speciationRate;
    int duration_speciation = candidate.Params.protractedNess; //  generations for speciation to happen
    double Xmean = candidate.Params.Xmean * n_lineages;

    while (remain_lineages > 1) { // do the coalescence until there are at least 2 ancestors
        steps ++;
        int Y0 = random_number(remain_lineages); //Random death, the function Rand now makes sure that Y0 is never equal to remain_lineages.
        int Y1 = 0;
        int lastPosition = (int)position.size()-1;
        int lastDescendant = (int)descendant.size()-1;

        if ((GenerationTime >= duration_speciation) &&
            (uniform() < speciationRate)) {
            speciation [descendant[Y0]] = true;
            result[descendant[Y0]] = species;
            index[position[lastPosition]] = Y0;
            index[position[Y0]] = -1;
            position[Y0] = position[lastPosition];
            descendant[Y0] = descendant[lastDescendant];
            position.pop_back();
            descendant.pop_back();
            ++ species;
        } else {
            int birth = 1;
            while (birth == 1) {
                double U = uniform(); // random variate U drawn from the uniform distribution on the unit interval (0, 1]
                Xmin = ceil(((candidate.Params.alpha - 1) / candidate.Params.alpha)* Xmean); // Xmin in terms of Xmean
                D = Xmin / powf((float)U, (float)(1.0f / (float)candidate.Params.alpha));
                double direction = uniform(); // here we only need to know whether the new individual comes from north or south //LOOK AT THIS!!

                if (direction < 0.5) {
                    Y1 = position[Y0] - D;
                } else {
                    Y1 = position[Y0] + D;
                }

                //********
                // Birth different scenarios
                // If there are soft borders affecting probability of reproduction softness must be set > 0.0 (the proportion of the lattice where there is a decrease of reprod.)
                if((Y1 > 0) && (Y1 <= (n_lineages-1))) { // is it still in the area, the sample must be from the whole result vector, does it need to be >=0??
                    birth = 0;
                    if(custom_mask == true)  {
                        if(mask[(int)Y1] == false) {
                            birth = 1; // if you can't live there, birth doesn't happpen
                        }
                    }
                }
            }

            /*************************************************************************************
             UPDATE THE VECTORS
             **************************************************************************************/
            // here I have to find out whether there is the number Y1 in my position vector OR noT
            int i = index[Y1];
            if (i > -1) {
                result[descendant[Y0]]=descendant[i];
                index[position[lastPosition]]=Y0; // changed order with next line
                index[position[Y0]] = -1;
                position[Y0]=position[lastPosition];
                descendant[Y0]=descendant[lastDescendant];
                position.pop_back();
                descendant.pop_back();
            } else  {
                index[position[Y0]] = -1;
                index[Y1] = Y0;
                position[Y0] = Y1;
            }//end else statement
        } // end else statement

        remain_lineages = (int)position.size(); // count the number of remaining lineages
        gtminus = 2 / double (remain_lineages);
        GenerationTime = GenerationTime + gtminus;
    }
    return;
}

void initialize_vectors(std::string mask_file_name,
                        std::vector<int>& index,
                        std::vector<int>& result,
                        std::vector<bool>& mask,
                        std::vector<int>& position,
                        std::vector<int>& descendant,
                        double sampling,
                        int& nExtraZeros,
                        int n_lineages,
                        bool custom_mask) {

    for(int i = 0; i < (int)index.size(); i ++) {
        index[i] = -2;
        result[i] = -2;
    }

    if(custom_mask == true) {
        read_mask(mask, mask_file_name, n_lineages);
        for (int i = 0; i < (int)mask.size(); ++i) {// synchronize the other linked vectors
            if (mask[i] == true) {
                if(uniform() <= sampling) {
                    position.push_back(i);
                    descendant.push_back(i);
                    index[i] = (int)position.size()-1;
                } else {
                    nExtraZeros++;
                }
            } else {
                nExtraZeros++;
            }
        }
    }

    if(custom_mask == false) {
        for (int i = 0; i < n_lineages; ++ i) {
            if(uniform() <= sampling) {
                mask.push_back(true);
            } else {
                mask.push_back(false);
                nExtraZeros ++; // check later
            }
        }

        for (int i = 0; i < (int)mask.size(); ++i) {// synchronize the other linked vectors
            if (mask[i] == true) {
                position.push_back(i);
                descendant.push_back(i);
                index[i] = (int)position.size()-1;
            }
        }
    }
    return;
}



void write_to_file(const GetParams& P,
                   const std::vector<int>& result,
                   int n_lineages,
                   const std::vector<double>& curve,
                   int species) {

    char sDataFileName1[80] = "";
    char sDataFileName2[80] = "";
    char sDataFileName3[80] = "";
    char sDataFileName4[80] = "";

    const char * sDataFileName = P.sDataFileName.c_str();

    strncpy(sDataFileName1, sDataFileName,80); // Raw data
    strncpy(sDataFileName2, sDataFileName,80); // Raw Metrics data
    strncpy(sDataFileName3, sDataFileName,80); // ICD of theoretical data
    strncpy(sDataFileName4, sDataFileName,80); // outcome results
    const std::string DataName1 = std::string(sDataFileName1) + "Raw_results_vector.txt";
    const std::string DataName2 = std::string(sDataFileName2) + "Metric_Data.txt";
    const std::string DataName3 = std::string(sDataFileName3) + "ICD_Data.txt";
    const std::string DataName4 = std::string(sDataFileName4) + "Outcome_Data.txt";

    const char * name1 = DataName1.c_str();
    const char * name2 = DataName2.c_str();
    const char * name3 = DataName3.c_str();
    const char * name4 = DataName4.c_str();

    output_all(result,
               n_lineages,
               name1,
               name2,
               name3,
               name4,
               P.Xmean, P.alpha, P.speciationRate, P.sampling, P.SoftBorder,
               P.focalReplicate);

}

//************************************************************************************************************************************************
//							Min function
// **********************************************************************************************************************************************
int mint (int value1, int value2) {
	if(value1 > value2) {
		return value2;
	} else {
		return value1;
	}
}
//************************************************************************************************************************************************
//							Mean function
// **********************************************************************************************************************************************

template <typename T>
double calculateMean(const std::vector<T>& v) {
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = 1.0 * sum / v.size();
	return(mean);
}

template <typename T>
std::vector <T> getUnique(const std::vector<T>& V) {
	std::vector<T> output;
	std::unique_copy(V.begin(),V.end(),std::back_inserter(output));
	return output;
}

template <typename T>
std::vector<int> findMatches(const std::vector<T>& V, T target) {
	std::vector<int> indices;
	for(int i = 0; i < (int)V.size(); ++i) {
		if(V[i] == target) indices.push_back(i);
	}
	return indices;
}


//****************************************************** Interpolate function thijs********************************************************

double interpolate(double x,
                   const std::vector< std::pair<double, double> >& table){
	static const double INF = 1.e100;
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	
	if (x > table.back().first) return 0;
	if (x < table[0].first) return 1;
	std::vector<std::pair<double, double> >::const_iterator it, it2;
	
	it = std::lower_bound(table.begin(), table.end(), std::make_pair(x, -INF));
	// Corner case
	if (it == table.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

//**********************************************************************************************************************************
//														Read Empirical Data Function
//***************************************************************************************************************************

//**********************************************************************************************************************************
//														Read Line Function that can read text files originating from multiple operating systems
//***************************************************************************************************************************

std::istream& safeGetline(std::istream& is,
                          std::string& t) {
    t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
            case '\n':
                return is;
            case '\r':
                if(sb->sgetc() == '\n')
                    sb->sbumpc();
                return is;
            case EOF:
                // Also handle the case when the last line has no line ending
                if(t.empty())
                    is.setstate(std::ios::eofbit);
                return is;
            default:
                t += (char)c;
        }
    }
}

std::vector<double> readEmpiricalData(std::string file_name) {
	//std::cout << "Reading data from file: " << file_name.c_str() << "\n";
	std::vector<double> output;
	double vectorX;
	double vectorY;
	
	std::vector<double> xvals;
	std::vector<double> yvals;
	
	std::ifstream myfile(file_name.c_str());
	
	if(!myfile.is_open()) {
		std::cout << "ERROR! Can't find empirical data!\n";
		return output;
	}
	
	if (myfile.is_open()) {
		std::string line;
        safeGetline(myfile, line); //we throw away the first line which is "X,Y"

		while (safeGetline(myfile, line)) {
			typedef boost::escaped_list_separator<char> Separator;
			typedef boost::tokenizer<Separator> Tokenizer;
			
			Tokenizer tokenizer(line);
			int count = 0;
			for (Tokenizer::iterator iter = tokenizer.begin();
                 (iter != tokenizer.end()) && (count < 3);
                 ++iter) {
				std::string::size_type sz;
				
				if (count == 0) {
					vectorX = std::stod (*iter,&sz); // well... at the end we dont need this
					xvals.push_back(vectorX);
				} else {
					vectorY = std::stod (*iter,&sz);
					yvals.push_back(vectorY);
				}
				++count;
			}
		}
	}
	
	if(xvals.size() > 0) {
	//	std::cout << "Empirical file Reading was a success!\n";
	}
	std::vector<std::pair<double, double> > emp;
	
	emp.clear(); //clear the vector, just in case there was still some data in it (shouldn't happen)
	for(int i = 0; i < (int)xvals.size(); ++i) {
		emp.push_back(std::pair<double, double>(xvals[i], yvals[i]));
	}
	
	std::vector<double> xlInterp;
	double aa = 0.5;
	for(int i = 1; i <= 200; ++ i) {
		xlInterp.push_back(aa);
		aa = aa + 0.5;
	}

    for(int i = 0; i < (int)xlInterp.size(); ++i) {// the LS is base on the empirical data range. usually Y values reach 0 before X is 100.
		double xval = xlInterp[i];
		double y_emp = interpolate(xval,emp); // emp[i];
		output.push_back(y_emp);
	}
	
	return output;
}

void read_mask(std::vector<bool>& mask,
               std::string mask_file_name,
               int n_lineages) {

    std::ifstream input_file(mask_file_name.c_str());
    while(!input_file.eof()) {
        int temp;
        input_file >> temp;
        mask.push_back(temp);
    }
    while(mask.size() > n_lineages) {
        mask.pop_back(); // extra entries can occur due to weird eof behaviour.
    }
    return;
}


//**********************************************************************************************************************************
//														Function to make curve output (for graphs later)
//***************************************************************************************************************************

std::vector<double> calcCurve(const std::vector<int>& FinalResults,
                              int gridSize,
                              int& numSpecies){
	if(FinalResults[0] == -500) {
		//pre-emptive exit
		std::vector<double> output(200,0.0);
		return output;
	}

	std::vector<int> temp = FinalResults;            // make a copy of spID vector
	std::vector<int> spID_uniques;
	std::vector<int> abundances;			// vector to store abundances as an output
	std::vector<int> Abs_Range;			// vector to store Latitudinal extent as an output
	std::vector<double> Rel_Range;       // vector to store ranges in percentage, useful for comparing with the empirical data.
	
	std::sort(temp.begin(),temp.end()); //sort spID vector... it is useful for deleting equal values to get the unique sp below.
	// it is better to sort the copy, because we need the actual order to link with XY coords
	std::unique_copy(temp.begin(),temp.end(), std::back_inserter(spID_uniques));

    numSpecies = (int)spID_uniques.size();

	/*************************************************************************************************************************/
	// Now we need to compare the uniques vector with the Full vector in order to get usefull information: abundances, range size metrics
	
	//spID_uniques.erase(spID_uniques.begin()+0); // remove the no sampled cells (-2) BETTER REMOVE IT MANUALLY
	//for (int i = 0; i < (int)spID_uniques.size(); ++ i)
	for(std::vector<int>::iterator i = spID_uniques.begin(); i != spID_uniques.end(); ++i) {
		int counter = 0;				// to sum abundances for each sp
		std::vector<int> Ycoord_focalsp;		// coordinate of each individual for each focal sp
		// extract XY coordinate for the focal species and stored in 2 new vectors
		
		// extract Y coordinate for the focal species and stored in 2 new vectors
		for(int j = 0; j < (int)FinalResults.size(); ++ j) {
			//if(spID_uniques[i] == FinalResults[j])
			if((*i) == FinalResults[j]) {
				++ counter;                    // sum the abundances
				Ycoord_focalsp.push_back(j);  // fill coordinates for a focal species i
			}
		}
		
		abundances.push_back(counter);        // fill abundances vector
		
		// estimate latitudinal and longitudinal extent using iterators
		std::vector<int>::iterator maxY, minY;
		maxY = max_element(Ycoord_focalsp.begin(),Ycoord_focalsp.end());
		minY = min_element(Ycoord_focalsp.begin(),Ycoord_focalsp.end());
		int range = *maxY - *minY;
		//double relRange = (double(range) / double(gridSize)) * 100;
		Abs_Range.push_back(range); // fill range size vectors
		//Rel_Range.push_back(relRange);
		
	} // end for loop for spID_uniques
	
	/***********************************************************************************************************************************
	 convert output in data for Inverse Cummulative Distribution
	 ***********************************************************************************************************************************/
	
	std::vector<int> l2;
    for(auto it = Abs_Range.begin(); it != Abs_Range.end(); ++it) if((*it) >= 0) {
        l2.push_back((*it)); // *it>=0 to include singletons
    }
	
	std::sort(l2.begin(),l2.end());
	
	std::vector<double> cuml; // dataICD; Y
	int N = (int)l2.size();
	for(int i = N; i >= 1; -- i) {
        double a = double(i); // / double(N);
		cuml.push_back(a);
	}

	std::vector<double> ALL_CDF_lat;
	std::vector<int>    new_xl;
	//first we need the unique range levels
	std::vector<int> unique_xl = getUnique(l2);
	
	//then we iterate
    for(std::size_t i = 0; i < unique_xl.size(); ++i) {
		std::vector<int> matches = findMatches(l2,unique_xl[i]);
		if(matches.size() == 1) {
			ALL_CDF_lat.push_back(cuml[matches[0]]);
		}
		
        if(matches.size() > 1){ //only store the lowest value (why?)

			double cdf = 1e20;
			
            for(std::size_t j = 0; j < matches.size(); ++j){
				if(cuml[matches[j]] < cdf) cdf = cuml[matches[j]];
			}
			ALL_CDF_lat.push_back(cdf);
		}
	}
	
	cuml = ALL_CDF_lat;
	
	//now we have two vectors, one with all the ranges, and one with the inverse cumulative distribution
	std::vector<double> xl; //lets make the ranges relative to the maximum range found
	int maxXL = 0;
	for(auto it = unique_xl.begin(); it != unique_xl.end(); ++it){
		if((*it) > maxXL) maxXL = (*it);
	}
	
	for(auto it = unique_xl.begin(); it != unique_xl.end(); ++it) {
		xl.push_back(100.0 * (*it) / maxXL);
	}

	// The xl values should be re-organize to have pre-defined intervals to be able to compare with empirical data using the LS
	std::vector<double> xlInterp;
	double aa = 0.5;
	for(int i = 1; i <= 200; ++ i) {
		xlInterp.push_back(aa);
		aa = aa + 0.5;
	}
	
	std::vector< std::pair<double, double> > sim;
	for(int i = 0; i < (int)xl.size(); ++i) {
		sim.push_back(std::pair<double,double>(xl[i],cuml[i]));
	}
	
	std::vector<double> output;
	
    for(int i = 0; i < (int)xlInterp.size(); ++i) {// the LS is base on the empirical data range. usually Y values reach 0 before X is 100.
		double xval = xlInterp[i];
		double y_sim = interpolate(xval,sim); // emp[i];
		output.push_back(y_sim);
	}
	return output;
}

//**********************************************************************************************************************************
//														Output function
//***************************************************************************************************************************

void output_all(std::vector<int> FinalResults,
                int gridSize,
                const char * name1,
                const char * name2,
                const char * name3,
				const char * name4,
                double Xmean,
                double alpha,
                double speciation,
                double samp,
                int SoftBorder,
                int replicate)

{
    std::ofstream ofDataFile1(name1, std::ios_base::app);
	std::ofstream ofDataFile2(name2, std::ios_base::app);
	std::ofstream ofDataFile3(name3, std::ios_base::app);
	std::ofstream ofDataFile4(name4, std::ios_base::app);
	
	/**********************************************************************************************************
	 OUTPUT & PRINTING
	 ***********************************************************************************************************/

    /**********************************************************************************************************
     OUTPUT RAW Results vector
     ***********************************************************************************************************/

    for(auto it = FinalResults.begin(); it != FinalResults.end(); ++it) {
        ofDataFile1 << (*it) << "\t";
    }
    ofDataFile1 << "\n";
    ofDataFile1.close();

	
	/***************************************************************************************************************
	 METRIC OUTPUTS
	 Abundance, latitudinal extent
	 ****************************************************************************************************************/
	std::cout<<"outputs \n";
	std::vector<int> spID_uniques;		// vector to store species ID
	std::vector<int> abundances;			// vector to store abundances as an output
	std::vector<int> Abs_Range;			// vector to store Latitudinal extent as an output
	std::vector<double> Rel_Range;       // vector to store ranges in percentage, useful for comparing with the empirical data.
	
	
	spID_uniques = FinalResults;            // make a copy of spID vector
	
	std::sort(spID_uniques.begin(),spID_uniques.end()); //sort spID vector... it is useful for deleting equal values to get the unique sp below.
	// it is better to sort the copy, because we need the actual order to link with XY coords
	
	// get unique values or species identities
    for (std::size_t spu = 0; spu < (spID_uniques.size() - 1); ++ spu) {
		if(spID_uniques[spu] == spID_uniques[spu + 1]) {
			spID_uniques.erase(spID_uniques.begin()+ spu); //erase the element i+1
			spu -= 1;	// if if is evaluated true, it has to be keep the same for continuing comparisons
		}
	}
	
	/*************************************************************************************************************************/
	// Now we need to compare the uniques vector with the Full vector in order to get usefull information: abundances, range size metrics
	
	//spID_uniques.erase(spID_uniques.begin()+0); // remove the no sampled cells (-2) BETTER REMOVE IT MANUALLY
	
    for (std::size_t i = 0; i < spID_uniques.size(); ++ i) {
		int counter = 0;				// to sum abundances for each sp
		std::vector<int> Ycoord_focalsp;		// coordinate of each individual for each focal sp
		// extract XY coordinate for the focal species and stored in 2 new vectors

		// extract Y coordinate for the focal species and stored in 2 new vectors
        for(std::size_t j = 0; j < FinalResults.size(); ++ j) {
			if(spID_uniques[i] == FinalResults[j]) {
				++ counter;                    // sum the abundances
				Ycoord_focalsp.push_back((int)j);  // fill coordinates for a focal species i
			}
		}
		
		abundances.push_back(counter);        // fill abundances vector
		
		// estimate latitudinal and longitudinal extent using iterators
		std::vector<int>::iterator maxY, minY;
		maxY = max_element(Ycoord_focalsp.begin(),Ycoord_focalsp.end());
		minY = min_element(Ycoord_focalsp.begin(),Ycoord_focalsp.end());
		int range = *maxY - *minY;
		range = range; //to avoid ranges of size 0
		double relRange = (double(range) / double(gridSize)) * 100;
		Abs_Range.push_back(range); // fill range size vectors
		Rel_Range.push_back(relRange);
		
	} // end for loop for spID_uniques
	
	/*****************************************************************************************************************************
	 save second output as a txt file:
	 *****************************************************************************************************************************/
	
    for (std::size_t i = 0; i < spID_uniques.size(); ++ i) {
		//ofDataFile2<<alpha<<","<<Xmean<<","<<spID_uniques[i]<<","<<abundances[i]<<","<<Lat_ext[i]<<","<<Lon_ext[i]<<","<<hullArea[i]<<endl;
		//ofDataFile2<<alpha<<","<<Xmean<<","<<repl<<","<<spID_uniques[i]<<","<<abundances[i]<<","<<Abs_Range[i]<<endl;
		//ofDataFile2<<alpha<<","<<Xmean<<","<<speciation<<","<<samp<<","<<sub_optimal<<","<<abundSubOptimal<<","<<spID_uniques[i]<<","<<abundances[i]<<","<<Rel_Range[i]<<std::endl;
		ofDataFile2<<alpha<<","<<Xmean<<","<<speciation<<","<<samp<<","<<SoftBorder<<","<<replicate<<","<<spID_uniques[i]<<","<<abundances[i]<<","<<Abs_Range[i]<<std::endl;

		//cout<<"species ID = "<<spID_uniques[i];
	}
	//ofDataFile1.close();
	ofDataFile2.close();

    int numSpecies;
    std::vector<double> output = calcCurve(FinalResults, gridSize, numSpecies);

    for(std::size_t sp4 = 0; sp4<output.size()-1;++sp4) {
		ofDataFile3<<alpha<<","<<Xmean<<","<<speciation<<","<<samp<<","<<SoftBorder<<","<<replicate<<","<<output[sp4]<<std::endl;
	}

	ofDataFile3.close();
    ofDataFile4.close();

	return;
}
