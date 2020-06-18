
//               hintzelab.msu.edu
//       for MABE documentation:
//               github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//       to view the full license, visit:
//               github.com/Hintzelab/MABE/wiki/License

#pragma once

#include "../AbstractWorld.h"

#include <cstdlib>
#include <thread>
#include <vector>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>

struct RankEpistasisData{
    size_t offset;
    size_t rank;
    double score_original;
    double score_mutant;
};


class NKWorld : public AbstractWorld {

public:

    // Parameters
    static std::shared_ptr<ParameterLink<int>> nPL;
    static std::shared_ptr<ParameterLink<int>> kPL;
    static std::shared_ptr<ParameterLink<int>> evaluationsPerGenerationPL;

    static std::shared_ptr<ParameterLink<bool>> readNKTablePL;
    static std::shared_ptr<ParameterLink<std::string>> inputNKTableFilenamePL;
    
    static std::shared_ptr<ParameterLink<bool>> writeNKTablePL;

    static std::shared_ptr<ParameterLink<bool>> treadmillPL;
    static std::shared_ptr<ParameterLink<double>> velocityPL; 

    static std::shared_ptr<ParameterLink<bool>> outputRankEpistasisPL; 
    static std::shared_ptr<ParameterLink<std::string>> outputRankEpistasisFilenamePL; 
    static std::shared_ptr<ParameterLink<int>> outputRankEpistasisIntervalPL; 
    static std::shared_ptr<ParameterLink<int>> outputEditDistanceMetricPL; 
    
    static std::shared_ptr<ParameterLink<bool>> outputMutantFitnessPL; 
    static std::shared_ptr<ParameterLink<std::string>> outputMutantFitnessFilenamePL; 
    static std::shared_ptr<ParameterLink<int>> outputMutantFitnessIntervalPL; 
    
    static std::shared_ptr<ParameterLink<std::string>> groupNamePL;
    static std::shared_ptr<ParameterLink<std::string>> brainNamePL;

    
    int N;
    int K;

    // Rank Epistasis output variables
    bool output_rank_epistasis;
    std::string output_rank_epistasis_filename;    
    int output_rank_epistasis_interval;
    int edit_distance_metric;
    std::stringstream output_string_stream;
    std::vector<size_t> rank_vec_original;  
    std::vector<size_t> rank_vec_mutated;  
    std::vector<RankEpistasisData> mutant_data_vec;  
    // Mutant fitness variables
    bool output_mutant_fitness;
    std::string output_mutant_fitness_filename;    
    int output_mutant_fitness_interval;

    std::vector<std::vector<std::pair<double,double>>> NKTable;

    NKWorld(std::shared_ptr<ParametersTable> PT_ = nullptr);
    virtual ~NKWorld() = default;

    // NK-specific functions
    double triangleSin(double x);

    void recordRankEpistasis(std::map<std::string, std::shared_ptr<Group>> &groups);
    void recordMutantFitness(std::map<std::string, std::shared_ptr<Group>> &groups);

    // evaluate functions
    double evaluateData(const std::vector<uint8_t>& data);
    double evaluateBrain(std::shared_ptr<AbstractBrain>& brain);
    void evaluateSolo(std::shared_ptr<Organism> org, int analyze, int visualize, int debug);
   
    virtual void evaluate(std::map<std::string, std::shared_ptr<Group>> &groups,
                                                int analyze, int visualize, int debug) {
        int popSize = groups[groupNamePL->get(PT)]->population.size();
        for (int i = 0; i < popSize; i++) {
            evaluateSolo(groups[groupNamePL->get(PT)]->population[i], analyze,
                                     visualize, debug);
        }
        if(output_rank_epistasis && Global::update % output_rank_epistasis_interval == 0)
            recordRankEpistasis(groups);
        if(output_mutant_fitness && Global::update % output_mutant_fitness_interval == 0)
            recordMutantFitness(groups);
    }

    virtual std::unordered_map<std::string, std::unordered_set<std::string>>
    requiredGroups() override {
        return {{groupNamePL->get(PT),
                         {"B:" + brainNamePL->get(PT) + ",0," +
                            std::to_string(nPL->get(PT))}}};
        // requires a root group and a brain (in root namespace) and no addtional
        // genome,
        // the brain must have 1 input, and the variable numberOfOutputs outputs
    }
};

