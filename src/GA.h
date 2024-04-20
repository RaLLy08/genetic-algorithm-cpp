#pragma once

#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <functional>
#include <sstream>
#include <random>

// #include "random.h"

namespace Random {
    int randInt(int min, int max);
    double randDouble(double min, double max);
    double norm_dist(double mean, double std_dev);
    void seed(int seed);
};

using namespace std;


class GA {
    private:
        const float mutationChance;

        const size_t chromosomeLength;
        const size_t maxGenerations;

        const size_t parentSize;
        const size_t offspringSize;
        const size_t populationSize;

        const size_t childCount;

        const size_t keepWorstSize;
        const size_t eliteParentSize;
        size_t generation;

    public:
        vector<vector<double>> population;

        function<double(double)> mutateGeneFunc;
        function<double(vector<double>&)> fitnessFunc;

        function<double(void)> initRandFunc;
        function<void(vector<vector<double>>&, size_t)> onGenerationEndFunc;
        function<void()> onFinish;


    GA(
        const size_t maxGenerations,
        const size_t chromosomeLength, 
        const size_t parentSize, 
        const float elitePercent = 0.1,
        const float mutationChance = 0.05,
        const float keepWorstPercent = 0.2,
        const size_t childCount = 1
    );

    void validateParameters();

    void breed();

    void replaceWeakToPopulationEnd();

    void run();

    private:
        vector<double> crossover(
            vector<double> *chromosome1, 
            vector<double> *chromosome2
        );

        void mutateChromosome(vector<double> *chromosome);

        void mutateChildren();
        
        vector<double> createChromosome();

        void initializeParents();

        void printParameters();
};
