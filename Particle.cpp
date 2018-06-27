//
//  Particle.cpp
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//
#include "randomc.h"
#include "Particle.h"
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////  Specific parameters for the model being fitted ////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


void theta::Perturb(double sigma)
{
    alpha = log10(alpha);
    alpha += normal(0,sigma);
    alpha = pow(10,alpha);
//	if(alpha < 1.01) alpha = 1.01;
//	if(alpha > 6) alpha = 6;


    Xmean = log10(Xmean);
	Xmean += normal(0,sigma);
    Xmean = pow(10,Xmean);
//	if(Xmean < 1e-6) Xmean = 1e-6;
//	if(Xmean > 0.56) Xmean = 0.56;


    speciationRate = log10(speciationRate);
    speciationRate += normal(0,sigma);
    speciationRate = pow(10,speciationRate);


    double temp = log10(protractedNess);
    temp += normal(0,sigma);
    protractedNess = (int)pow(10,temp);

   // if(protractedNess < 0) protractedNess = 0;


    sampling = log10(sampling);
    sampling += normal(0,sigma);
    sampling = pow(10,sampling);

  //  if(sampling > 1) sampling = 1;
  //  if(sampling < 1e-4) sampling = 1e-4;


	return;
}

std::vector<double> theta::operator-(const theta& other)
{
	std::vector<double> output;
	output.push_back(log10(other.alpha) - log10(alpha));
	output.push_back(log10(other.Xmean) - log10(Xmean));
    output.push_back(log10(other.speciationRate) - log10(speciationRate));
    output.push_back(log10(other.protractedNess) - log10(protractedNess));
    output.push_back(log10(other.sampling) - log10(sampling));
	return output;
}


std::istream& operator >> (std::istream& is, theta& t)
{
	is >> t.alpha;
	is >> t.Xmean;
    is >> t.speciationRate;
    is >> t.protractedNess;
    is >> t.sampling;
	return is;
}

std::ostream& operator << (std::ostream& os, const theta& t)
{
	os << t.alpha << "\t";
	os << t.Xmean << "\t";
    os << t.speciationRate << "\t";
    os << t.protractedNess << "\t";
    os << t.sampling << "\t";
	return os;
}


bool particle::accept(double threshold, int iter)
{
	if(iter == 0) return true;
	if(fitness < threshold) return true;
	return false;
}


//void particle::calcFitness(std::vector<double> emp, std::vector<double> Sim, int emp_num_species)
void particle::calcFitness(const std::vector<double>& emp, const std::vector<double>& Sim)
{
	fitness = 0.0;
    for(std::size_t i = 0; i < emp.size(); ++i)
	{
		double diff = emp[i] - Sim[i];
		fitness += diff * diff;
	}

    // calculate difference of number of species
    //double diff = (emp_num_species - numSpecies) * (emp_num_species - numSpecies) ;

  //  diff = 1.0 * diff / (57.78812 * 57.78812); // 57.78812 is the standard deviation of number of species in the empirical data. This is to "normalize" the difference
   // diff = 1.0 * diff / (30 * 30);  //when generating from the prior, the standard deviation in number of species (excluding extremely high numbers) is ~30
  //  fitness += diff;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//// General functions for SMC particle behaviour //////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void particle::getFromPrevious(const std::vector<particle>& P, double maxW)
{
    if(maxW < 0) {
        std::cout << "ERROR MAXW < 0\n";
        exit(1);
    }


	int i = 0;
	int maxTries = 1e9;
	while(i < maxTries) {
		int index= random_number((int)P.size());
		double acceptRate = 1.0 * P[index].weight / maxW;
		if(uniform() < acceptRate)
		{
			(*this) = P[index];
			return;
		}
		i++;
	}
	
	
	(*this) = P[0]; //failsafe, should NEVER be reached!
	return;
}


void particle::perturb(double sigma)
{
	Params.Perturb(sigma);
	return;
}


void particle::calculateWeight(const std::vector<particle>& P, double sigma)
{
	double sum = 0.0;
	static double SQRT2PI = sqrt(2*3.14159265359);
	static double mult = 1.0 / (sigma * SQRT2PI);
	static double divisor = 2 * sigma * sigma;
	for(std::vector<particle>::const_iterator p = P.begin(); p != P.end(); ++p)
	{
		
		std::vector<double> diff = Params - (*p).Params;
				
		double prod = (*p).weight;
        for(std::size_t i = 0; i < diff.size(); ++i)
		{
			prod *= 1.0 * mult * exp(-1.0 * (diff[i]*diff[i]) / divisor);
		}
		
		sum += prod;
	}
	
	weight = 1.0 / sum;
}

void normalize(std::vector< particle >& P, double& maxW)
{
	maxW = -5;
	double sum = 0.0;
    int counter = 0;
	for(std::vector< particle>::iterator it = P.begin(); it != P.end(); ++it)   {
		sum += (*it).weight;
        counter++;
	}
	
	double mult = 1.0 / sum;
	
	for(std::vector< particle>::iterator it = P.begin(); it != P.end(); ++it)
	{
		(*it).weight= (*it).weight * mult;
		if((*it).weight > maxW) maxW = (*it).weight;
	}
	return;
}


bool file_exists(const std::string& file_name)
{
	std::ifstream f(file_name.c_str());
	if (f.good()) {
		f.close();
		return true;
	} else {
		f.close();
		return false;
	}
}

void readParticles(int t, std::vector<particle>& particles, int numberParticles)
{
	//std::string tt = boost::lexical_cast<std::string>(t);
	std::stringstream ss;
	ss << t;
	std::string tt = ss.str();
	std::string file_name = "particles_t=" + tt + ".txt";

    if(file_exists(file_name)) {

        std::ifstream read_part(file_name.c_str());

        std::cout << "attempting to read \t" << file_name << "\n";

        if(!read_part.is_open()) {
            std::cout << "No input file of particles found!!!!\n";
            return;
        }
        int counter = 0;
        while(!read_part.eof()) {
            particle temp;
            read_part >> temp;

            particles.push_back((temp));
            counter++;

            if(counter > numberParticles + 1) {
                break;
            }
        }
    }
    std::cout << "Done reading particles\n";
    while((int)particles.size() > numberParticles) {
        particles.pop_back();
    }

	return;
}



particle::particle(const particle& other)
{
	weight = other.weight;
	Params = other.Params;
	fitness = other.fitness;
    numSpecies = other.numSpecies;
}


particle& particle::operator=(const particle& other)
{
	if(this == &other) return *this;
	
	weight = other.weight;
	Params = other.Params;
	fitness = other.fitness;
    numSpecies = other.numSpecies;
	
	return *this;
}



std::istream& operator >> (std::istream& is, particle& p)
{
	is >> p.Params;
	is >> p.weight;
	is >> p.fitness;
    is >> p.numSpecies;
	return is;
}

std::ostream& operator << (std::ostream& os, const particle& p)
{
	os << p.Params << "\t";
	os << p.weight   << "\t";
	os << p.fitness << "\t";
    os << p.numSpecies << "\t";
	return os;
}


