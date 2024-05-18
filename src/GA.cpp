#include "GA.h"

namespace Random {
    mt19937 gen;

    int randInt(int min, int max) {
        uniform_int_distribution<int>distribution(min, max);
        return distribution(gen);
    }

    double randDouble(double min, double max) {
        uniform_real_distribution<double> distribution(min, max);
        return distribution(gen);
    }

    double norm_dist(double mean, double std_dev) {
        normal_distribution<double> dist(mean, std_dev);

        return dist(gen);
    };

    void seed(int seed) {
        gen.seed(seed);
    }
}

GA::GA(
    const size_t maxGenerations,
    const size_t chromosomeLength, 
    const size_t parentSize, 
    const float elitePercent,
    const float mutationChance,
    const float keepWorstPercent,
    const size_t childCount
): 
    chromosomeLength(chromosomeLength), 
    maxGenerations(maxGenerations),
    parentSize(parentSize),
    childCount(childCount),
    offspringSize(parentSize * childCount),
    eliteParentSize(parentSize * elitePercent),
    keepWorstSize(offspringSize * keepWorstPercent),
    populationSize(parentSize + offspringSize),
    mutationChance(mutationChance),
    population(populationSize)
{}

void GA::validateParameters() {
    stringstream ss;

    if (initRandFunc == nullptr) {
        throw invalid_argument("initRandFunc is not set");
    }

    if (this->fitnessFunc == nullptr) {
        throw invalid_argument("fitnessFunc is not set");
    }

    if (keepWorstSize > parentSize - eliteParentSize) {
        ss << "keepWorstSize(" << keepWorstSize << ") is less than parentSize(" << parentSize << ") - eliteParentSize(" << eliteParentSize << ")";

        throw invalid_argument(ss.str());
    }
}

void GA::breed() {
    // variable child count
    
    for (size_t j = 0; j < parentSize; j++) {
        size_t parent2Index = Random::randInt(0, parentSize - 1);

        vector<double> *parent1 = &population[j];
        vector<double> *parent2 = &population[parent2Index];

        vector<double> child = crossover(parent1, parent2);

        population[parentSize + j] = child;
    }
}

void GA::replaceWeakToPopulationEnd() {
    shuffle(population.begin() + parentSize, population.end(), Random::gen);

    for (size_t wI = 0; wI < keepWorstSize; wI++) {
        // starts from end of population -> new children size
        size_t parentIndex = (parentSize - 1) - wI;
        size_t offspringRandIndex = parentSize + wI;

        swap(
            population[parentIndex], 
            population[offspringRandIndex]
        );
    }
}

void GA::run() {
    validateParameters();

    initializeParents();

    for (generation = 0; generation < maxGenerations; generation++) {
        breed();
        mutateChildren();

        sort(population.begin(), population.end(), [this](vector<double> a, vector<double> b) {
            return fitnessFunc(a) < fitnessFunc(b);
        });

        if (onGenerationEndFunc != nullptr) {
            onGenerationEndFunc(population, generation);
        }

        replaceWeakToPopulationEnd();
    }

    onFinish();
    printParameters();
}

vector<double> GA::crossover(
    vector<double> *chromosome1, 
    vector<double> *chromosome2
) {
    int type = Random::randInt(0, 1);

    vector<double> child = *chromosome1;

    if (type == 0) {
        size_t crossoverPoint = Random::randInt(0, chromosomeLength - 1);

        for (size_t i = crossoverPoint; i < chromosomeLength; i++) {
            child[i] = (*chromosome2)[i];
        }
    }

    if (type == 1) {
        for (size_t i = 0; i < chromosomeLength; i++) {
            if (i % 2 == 0) {
                child[i] = (*chromosome2)[i];
            }
        }    
    }

    return child;
}

void GA::mutateChromosome(vector<double> *chromosome) {
    size_t mutationPoint = Random::randInt(0, chromosomeLength - 1);
    auto* mutValueRef = &(*chromosome)[mutationPoint];

    *mutValueRef = mutateGeneFunc(*mutValueRef);
}

void GA::mutateChildren() {
    for (size_t i = parentSize; i < populationSize; i++) {
        if (Random::randDouble(0.0, 1.0) < mutationChance) {
            mutateChromosome(&population[i]);
        }
    }
}

vector<double> GA::createChromosome() {
    vector<double> chromosome(chromosomeLength);

    for (double& gene : chromosome) {
        gene = initRandFunc();
    }

    return chromosome;
}

void GA::initializeParents() {
    for (size_t i = 0; i < parentSize; i++) {
        vector<double> chromosome = createChromosome();

        population[i] = chromosome;
    }
}

void GA::printParameters() {
    cout << endl << "--------------------------------" << endl;
    cout << "Max generations: " << maxGenerations << endl;
    cout << "Chromosome length: " << chromosomeLength << endl;
    cout << "Population size: " << populationSize << endl;
    cout << endl;
    cout << "Offspring size: " << offspringSize << endl;
    cout << "Parent size: " << parentSize << endl;
    cout << endl;
    cout << "Elite Parent Size: " << eliteParentSize << endl;
    cout << "Keep worst percent: " << keepWorstSize << endl;
    cout << "Mutation chance: " << mutationChance << endl;

    cout << "Size of population in bytes: " << (population.size() * chromosomeLength * sizeof(double))  << endl;

    cout << "--------------------------------" << endl;
}
