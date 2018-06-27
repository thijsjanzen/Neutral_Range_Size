//
//  Particle.h
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#ifndef __SecondaryContact__Particle__
#define __SecondaryContact__Particle__

#include <vector>
#include "randomc.h"
#include <string>
#include <cmath>

//parameter structure
struct theta
{
	double alpha;
	double Xmean;
    double speciationRate;
    double protractedNess;
    double sampling;

    theta() : alpha(1.5),
              Xmean(0.5),
              speciationRate(1.0),
              protractedNess(0.0),
              sampling(1.0) {

	}
	
	theta(double Alpha, double xMean,
          double specR, double tau,
          double sample) : alpha(Alpha),
                           Xmean(xMean),
                           speciationRate(specR),
                           protractedNess(tau),
                           sampling(sample) {

	}
	
	//prior functions
	void getFromPrior() {
		alpha = 1.01 + uniform() * (6 - 1.01);
		Xmean = uniform() * 0.56;
        speciationRate = -5 + uniform() * 5; // 0.0001 to 1
        speciationRate = powf(10,speciationRate);
        protractedNess = random_number(1e5);
        sampling = 1e-4 + uniform() * (1-1e-4); //to ensure that it doesn't get too small

		return;
	}

    bool withinPrior()  {
        if(alpha < 1.01) return false;
        if(alpha > 6) return false;
        if(Xmean < 0) return false;
        if(Xmean > 0.56) return false;

        if(speciationRate < 1e-5) return false;
        if(speciationRate > 1) return false;

        if(protractedNess < 0) return false;
        if(protractedNess > 1e5) return false;

        if(sampling < 1e-4) return false;
        if(sampling > 1-1e-4) return false;

        return true;
    }
	
	std::vector<double> operator-(const theta& other);

	void Perturb(double sigma);
};


//particle structure, stores parameter values, weight, and fit
struct particle
{
	theta Params;
	double weight;
	double fitness;
    int numSpecies;
	
	
	
    particle() : fitness(-1),
                 weight(1),
                 numSpecies(0) {
		Params = theta();
	}
	
    particle(const theta& init) : Params(init),
                                  fitness(-1),
                                  weight(1),
                                  numSpecies(0)
	{
	}
	
	particle(const particle& other);
	particle & operator=(const particle& other);

	void getFromPrevious(const std::vector< particle >& P, double maxW);
	void perturb(double sigma);
	void calculateWeight(const std::vector<particle>& P, double sigma);
	bool accept(double threshold, int iter);
	
    void calcFitness(const std::vector<double>& emp, const std::vector<double>& Sim);
	
	void getFromPrior() {
		Params.getFromPrior();
	}

    bool withinPrior() {
        return(Params.withinPrior());
    }

};

bool file_exists(const std::string& file_name);
void normalize(std::vector< particle >& P, double& maxW);

std::istream& operator >> (std::istream& is, particle& p);
std::ostream& operator << (std::ostream& os, const particle& p);
std::istream& operator >> (std::istream& is, theta& t);
std::ostream& operator << (std::ostream& os, const theta& t);
void readParticles(int t, std::vector<particle>& parts, int numberParticles);

#endif /* defined(__SecondaryContact__Particle__) */
