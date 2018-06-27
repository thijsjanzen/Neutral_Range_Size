//
//  Adriana_simulation.h
//  SMC_Adriana
//
//  Created by Thijs Janzen on 03/12/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#ifndef __SMC_Adriana__Adriana_simulation__
#define __SMC_Adriana__Adriana_simulation__

#include <stdio.h>
#include <vector>
#include "Particle.h"
#include "GetParams.h"

int mint (int value1, int value2);
double interpolate(double x, const std::vector< std::pair<double, double> >& table);

std::vector<double> readEmpiricalData(std::string file_name);

double calcFitness(std::vector<double> emp, std::vector<double> Sim);

std::vector<double> calcCurve(const std::vector<int>& FinalResults, int gridSize, int& numSpecies);
std::vector<double> doSimulation(particle candidate, const GetParams& P, int& numSpecies);

void initialize_vectors(std::vector<int>& index,
                        std::vector<int>& result,
                        std::vector<bool>& mask,
                        std::vector<int>& position,
                        std::vector<int>& descendant,
                        double sampling,
                        int& nExtraZeros,
                        int n_lineages);

void simulate_model(particle candidate,
                    const std::vector<int>& position,
                    std::vector<int>& descendant,
                    std::vector<int>& result,
                    const std::vector<int>& index,
                    const std::vector<bool>& speciation,
                    int& species,
                    const int& n_lineages);

void write_to_file(const GetParams& P,
                   const std::vector<int>& result,
                   int n_lineages,
                   const std::vector<double>& curve,
                   int species);

void output_all(std::vector<int> FinalResults,int gridSize, const char * name2,const char * name3,
                const char * name4, double Xmean, double alpha, double speciation,double samp, int SoftBorder, int richness, std::string emp_file_name, std::vector<double> curve, int replicate);


template <typename T>
double calculateMean(const std::vector<T>& v);
template <typename T>
std::vector<int> findMatches(const std::vector<T>& V, T target);
template <typename T>
std::vector <T> getUnique(const std::vector<T>& V);










#endif /* defined(__SMC_Adriana__Adriana_simulation__) */
