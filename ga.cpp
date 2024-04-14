#include <stdlib.h>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <sys/resource.h>
#include <cassert>

#define MUTATION_CHANCE_PERCENT 5

double norm_dist(double mean, double std_dev);
template <typename T>
T find_median(T arr[], size_t arr_size);

// class Rand

using namespace std;
mt19937 gen;

template <typename E>
class GA {
    size_t populationSize;
    size_t chromosomeLength;
    size_t maxGenerations;
    vector<double> fitness;
    vector<vector<E>> population;

    function<E(void)> chromoInitRandLambda;

    public:
        static int randInt(int min, int max) {
            uniform_int_distribution<int> distribution(min, max);
            return distribution(gen);
        }

        static double randDouble(double min, double max) {
            uniform_real_distribution<double> distribution(min, max);
            return distribution(gen);
        }

        void setChromoInitRandLambda(function<E(void)> chromoInitRandLambda) {
            this->chromoInitRandLambda = chromoInitRandLambda;
        }

        GA(
            size_t populationSize, 
            size_t chromosomeLength, 
            size_t maxGenerations
        ): 
            populationSize(populationSize), 
            chromosomeLength(chromosomeLength), 
            maxGenerations(maxGenerations),
            fitness(populationSize, 0),
            population(populationSize * 2)
        {}

        double fitnessFunction(vector<E> chromosome) {
            double sum = 0;

            for (E gene : chromosome) {
                sum += gene;
            }

            return sum;
        }

        void run() {
            validateParameters();
            printParameters();

            initializePopulation();
            printPopulation();

            size_t childNumber = 1;
            size_t eliteSize = 2;

            for (size_t i = 0; i < maxGenerations; i++) {
                
                for (size_t j = 0; j < populationSize; j++) {
                    size_t parent2Index = rand() % populationSize;

                    vector<E> parent1 = population[j];
                    vector<E> parent2 = population[parent2Index];

                    vector<E> child = crossover(parent1, parent2);

                    if (rand() % 100 < MUTATION_CHANCE_PERCENT) {
                        mutate(child);
                    }

                    population[populationSize + j] = child;
                }

                // calculate fitness separately
                sort(population.begin(), population.end(), [this](vector<E> a, vector<E> b) {
                    return fitnessFunction(a) > fitnessFunction(b);
                });

                // float randomWorstPercent = 0.2;

                size_t pickMaxWorstSize = 3;
    
                assert(("keepWorstMax must be less than populationSize - eliteSize", pickMaxWorstSize < populationSize - eliteSize));
                // unique random
                shuffle(population.begin() + populationSize, population.end(), gen);

                for (size_t wI = 0; wI < pickMaxWorstSize; wI++) {
                    // starts from end of population -> new children size
                    size_t parentIndex = (populationSize - 1) - wI;
                    size_t offspringRandIndex = populationSize + wI;

                    swap(
                        population[parentIndex], 
                        population[offspringRandIndex]
                    );
                }
            }

            cout << "Final population" << endl;
            printPopulation();

            cout << "Best chromosome: " << endl;
            printChromosome(population[0]);
        }

    private:
        void validateParameters() {
            if (this->chromoInitRandLambda == nullptr) {
                throw invalid_argument("chromoInitRandLambda is not set");
            }
        }

        vector<E> crossover(vector<E> chromosome1, vector<E> chromosome2) {
            size_t crossoverPoint = rand() % chromosomeLength;
            vector<E> child = chromosome1;

            for (size_t i = crossoverPoint; i < chromosomeLength; i++) {
                child[i] = chromosome2[i];
            }

            return child;
        }

        void mutate(vector<E> chromosome) {
            size_t mutationPoint = rand() % chromosomeLength;

            chromosome[mutationPoint] = getMutationValue(chromosome[mutationPoint]);
        }

        E getMutationValue(E value) {
            // change depenging on fitness etc.
            return value + norm_dist(0, 0.1);
        }

        vector<E> createChromosome() {
            vector<E> chromosome(chromosomeLength);

            for (E& gene : chromosome) {
                gene = chromoInitRandLambda();
            }

            return chromosome;
        }

        void initializePopulation() {
            cout << "Initializing population" << endl;

            for (size_t i = 0; i < populationSize; i++) {
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
            cout << "Population size: " << populationSize << endl;
            cout << "Chromosome length: " << chromosomeLength << endl;
            cout << "Max generations: " << maxGenerations << endl;
            cout << "Size of population in bytes: " << (population.size() * chromosomeLength * sizeof(E))  << endl;
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

    ga.setChromoInitRandLambda([]() -> int {
        return rand() % 2;
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

