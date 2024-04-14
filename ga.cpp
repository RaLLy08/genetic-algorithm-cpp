#include <stdlib.h>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <sys/resource.h>
#include <cassert>

using namespace std;

double norm_dist(double mean, double std_dev);
template <typename T>
T find_median(T arr[], size_t arr_size);


mt19937 gen;


static int randInt(int min, int max) {
    uniform_int_distribution<int>distribution(min, max);
    return distribution(gen);
}

static double randDouble(double min, double max) {
    uniform_real_distribution<double> distribution(min, max);
    return distribution(gen);
}

template <typename E>
class GA {
    const float mutationChance;

    const size_t chromosomeLength;
    const size_t maxGenerations;

    const size_t parentSize;
    const size_t offspringSize;
    const size_t populationSize;

    const size_t childCount;

    const size_t keepWorstSize;
    const size_t eliteParentSize;

    vector<double> fitness;
    vector<vector<E>> population;

    function<E(void)> initRandFunc;
    function<E(E)> mutateRandFunc;

    public:
        void setInitRandFunc(function<E(void)> initRandFunc) {
            this->initRandFunc = initRandFunc;
        }

        void setMutateRandFunc(function<E(E)> mutateRandFunc) {
            this->mutateRandFunc = mutateRandFunc;
        }


        void validateParameters() {
            if (this->initRandFunc == nullptr) {
                throw invalid_argument("initRandFunc is not set");
            }
        }

        /*
            @param parentSize: number of parents
            @param chromosomeLength: length of chromosome
            @param maxGenerations: number of generations
            @param childCount: number of children for each parent
        */
        GA(
            const size_t parentSize, 
            const size_t chromosomeLength, 
            const size_t maxGenerations,
            const float elitePercent = 0.1,
            const float mutationChance = 0.05,
            const float keepWorstPercent = 0.2,
            const size_t childCount = 1
        ): 
            chromosomeLength(chromosomeLength), 
            maxGenerations(maxGenerations),
            parentSize(parentSize),
            childCount(childCount),
            offspringSize(parentSize * childCount),
            // fitness(parentSize),
            eliteParentSize(parentSize * elitePercent),
            keepWorstSize(offspringSize * keepWorstPercent),
            populationSize(parentSize + offspringSize),
            mutationChance(mutationChance),
            population(populationSize)
        {}

        double fitnessFunction(vector<E> chromosome) {
            double sum = 0;

            for (E gene : chromosome) {
                sum += gene;
            }

            return sum;
        }

        void breed() {
            // variable child count
            for (size_t j = 0; j < parentSize; j++) {
                size_t parent2Index = randInt(0, parentSize - 1);

                vector<E> parent1 = population[j];
                vector<E> parent2 = population[parent2Index];

                vector<E> child = crossover(parent1, parent2);

                population[parentSize + j] = child;
            }
        }

        void replaceWeakToPopulation() {
            assert((
                "keepWorstSize must be less than parentSize - eliteParentSize", 
                keepWorstSize < parentSize - eliteParentSize
            ));

            shuffle(population.begin() + parentSize, population.end(), gen);

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

        // void replaceWeakToPopulationRandom 

        void run() {
            validateParameters();

            initializePopulation();
            printPopulation();

            for (size_t i = 0; i < maxGenerations; i++) {
                breed();
                mutateChildren();

                // calculate fitness separately
                sort(population.begin(), population.end(), [this](vector<E> a, vector<E> b) {
                    return fitnessFunction(a) > fitnessFunction(b);
                });

                // float randomWorstPercent = 0.2;
                replaceWeakToPopulation();
            }

            cout << "Final population" << endl;
            printPopulation();

            cout << "Best chromosome: " << endl;
            printChromosome(population[0]);

            printParameters();
        }

    private:


        vector<E> crossover(vector<E> chromosome1, vector<E> chromosome2) {
            size_t crossoverPoint = randInt(0, chromosomeLength - 1);
            vector<E> child = chromosome1;

            for (size_t i = crossoverPoint; i < chromosomeLength; i++) {
                child[i] = chromosome2[i];
            }

            return child;
        }

        void mutateChromosome(vector<E> chromosome) {
            size_t mutationPoint = randInt(0, chromosomeLength - 1);

            chromosome[mutationPoint] = getMutationValue(chromosome[mutationPoint]);
        }

        void mutateChildren() {
            for (size_t i = parentSize; i < populationSize; i++) {
                if (randDouble(0.0, 1.0) < mutationChance) {
                    mutateChromosome(population[parentSize]);
                }
            }
        }

        E getMutationValue(E value) {
            // change depenging on fitness etc.
            return mutateRandFunc(value);
        }

        vector<E> createChromosome() {
            vector<E> chromosome(chromosomeLength);

            for (E& gene : chromosome) {
                gene = initRandFunc();
            }

            return chromosome;
        }

        void initializePopulation() {
            cout << "Initializing population" << endl;

            for (size_t i = 0; i < parentSize; i++) {
                vector<E> chromosome = createChromosome();

                population[i] = chromosome;
            }
        }

        void printChromosome(vector<E> chromosome) {
            for (auto it = chromosome.begin(); it != chromosome.end(); ++it) {
                cout << *it << " ";
            }
        }

        void printPopulation() {
            for (auto c = population.begin(); c != population.end(); ++c) {
                vector<E> chromosome = *c;

                printf("Chromosome: %ld:\t", c - population.begin());

                printChromosome(chromosome);

                cout << endl;
            }
        }
        
        void printParameters() {
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

            cout << "Size of population in bytes: " << (population.size() * chromosomeLength * sizeof(E))  << endl;

            cout << "--------------------------------" << endl;

            
        }
};

void printMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);

    cout << "Memory usage: " << usage.ru_maxrss << " kb" << endl;
}


int main() {
    // gen.seed(random_device()());
    // srand(time(NULL));

    GA<int> ga(10, 10, 20);

    ga.setInitRandFunc([]() -> int {
        return randInt(0, 1);
    });

    ga.setMutateRandFunc([](int value) -> int {
        //  value + norm_dist(0, 0.1);
        return value == 0 ? 1 : 0;
    });

    ga.run();

    printMemoryUsage();

    return 0;
}

double norm_dist(double mean, double std_dev) {
    normal_distribution<double> dist(mean, std_dev);

    return dist(gen);
};

template <typename T>
T find_median(T arr[], size_t arr_size) {
    T* arr_copy = new T[arr_size];
    for (int i = 0; i < arr_size; i++) {
        arr_copy[i] = arr[i];
    }

    sort(arr_copy, arr_copy + arr_size);

    T median;
    
    if (arr_size % 2 == 0) {
        median = (arr_copy[arr_size / 2 - 1] + arr_copy[arr_size / 2]) / 2;
    } else {
        median = arr_copy[arr_size / 2];
    }

    delete[] arr_copy;

    return median;
}

