#include <stdlib.h>
#include <iostream>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <sys/resource.h>
#include <cassert>
#include <sstream>

#include <thread>
#include <chrono>

using namespace std;

mt19937 gen;

static int randInt(int min, int max) {
    uniform_int_distribution<int>distribution(min, max);
    return distribution(gen);
}

static double randDouble(double min, double max) {
    uniform_real_distribution<double> distribution(min, max);
    return distribution(gen);
}


static double norm_dist(double mean, double std_dev) {
    normal_distribution<double> dist(mean, std_dev);

    return dist(gen);
};

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
    size_t generation = 0;

    public:
        vector<vector<E>> population;

        function<E(E)> mutateGeneFunc;
        function<E(vector<E>*)> fitnessFunc;

        function<E(void)> initRandFunc;
        function<void(vector<vector<E>>*, size_t)> onGenerationEndFunc;
        function<void()> onFinish;


        void validateParameters() {
            stringstream ss;


            if (this->initRandFunc == nullptr) {
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

        /*
            @param parentSize: number of parents
            @param chromosomeLength: length of chromosome
            @param maxGenerations: number of generations
            @param childCount: number of children for each parent
            @param elitePercent: percentage of best parents that will not be mutated and replaced
        */
        GA(
            const size_t maxGenerations,
            const size_t chromosomeLength, 
            const size_t parentSize, 
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

        void breed() {
            // variable child count
            for (size_t j = 0; j < parentSize; j++) {
                size_t parent2Index = randInt(0, parentSize - 1);

                vector<E> *parent1 = &population[j];
                vector<E> *parent2 = &population[parent2Index];

                vector<E> child = crossover(parent1, parent2);

                population[parentSize + j] = child;
            }
        }

        void replaceWeakToPopulationEnd() {
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

            initializeParents();

            for (generation = 0; generation < maxGenerations; generation++) {
                breed();
                mutateChildren();

                sort(population.begin(), population.end(), [this](vector<E> a, vector<E> b) {
                    return fitnessFunc(&a) < fitnessFunc(&b);
                });


                if (onGenerationEndFunc != nullptr) {
                    onGenerationEndFunc(&population, generation);
                }

                replaceWeakToPopulationEnd();
            }

            onFinish();
            printParameters();
        }

    private:
        vector<E> crossover(
            vector<E> *chromosome1, 
            vector<E> *chromosome2
        ) {
            int type = randInt(0, 1);

            vector<E> child = *chromosome1;

            if (type == 0) {
                size_t crossoverPoint = randInt(0, chromosomeLength - 1);

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

        void mutateChromosome(vector<E> *chromosome) {
            size_t mutationPoint = randInt(0, chromosomeLength - 1);
            auto* mutValueRef = &(*chromosome)[mutationPoint];
    
            *mutValueRef = mutateGeneFunc(*mutValueRef);
        }

        void mutateChildren() {
            for (size_t i = parentSize; i < populationSize; i++) {
                if (randDouble(0.0, 1.0) < mutationChance) {
                    mutateChromosome(&population[i]);
                }
            }
        }

        vector<E> createChromosome() {
            vector<E> chromosome(chromosomeLength);

            for (E& gene : chromosome) {
                gene = initRandFunc();
            }

            return chromosome;
        }

        void initializeParents() {
            for (size_t i = 0; i < parentSize; i++) {
                vector<E> chromosome = createChromosome();

                population[i] = chromosome;
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

double rosensbrock(double a, double b) {
    return 100 * pow(b - a * a, 2) + pow(1 - a, 2);
}

double sphere(double a, double b) {
    return a * a + b * b;
}

double deJong(double a, double b) {
    return 3905.93 - 100 * (a * a - b) * (a * a - b) - (1 - a) * (1 - a);
}

double venter(double a, double b) {
    return 100 * pow(b - a * a, 2) + pow(1 - a, 2);
}

double ackley(double a, double b) {
    return -20 * exp(-0.2 * sqrt(0.5 * (a * a + b * b))) - exp(0.5 * (cos(2 * M_PI * a) + cos(2 * M_PI * b))) + 20 + M_E;
}

double schwefel(double a, double b) {
    return 418.9829 * 2 - (a * sin(sqrt(abs(a))) + b * sin(sqrt(abs(b))));
}

double griewank(double a, double b) {
    return 1 + (a * a + b * b) / 4000 - cos(a) * cos(b / sqrt(2));
}

double styblinski(double a, double b) {
    return 0.5 * (a * a * a * a - 16 * a * a + 5 * a + b * b * b * b - 16 * b * b + 5 * b);
}

// test equations
double equation(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j) {
    // return a * x * x * x + b * x * x + c * x;

    // return a * a + b * b + c * c + d * d;

    return rosensbrock(a, b) + sphere(c, d) + deJong(e, f) + venter(g, h) + ackley(i, j) + schwefel(a, b) + griewank(c, d) + styblinski(e, f);
}


double desiredResult = equation(4, 10, 6, 7, 8, 9, 10, 11, 12, 13);

void printChromosome(vector<double> *chromosome) {
    for (auto it = chromosome->begin(); it != chromosome->end(); ++it) {
        cout << *it << " ";
    }

    cout << endl;
}

double calcResult(vector<double> *chromosome) {
    return equation(
        chromosome->at(0),
        chromosome->at(1),
        chromosome->at(2),
        chromosome->at(3),
        chromosome->at(4),
        chromosome->at(5),
        chromosome->at(6),
        chromosome->at(7),
        chromosome->at(8),
        chromosome->at(9)
    );

}

double fitnessFunc(vector<double> *chromosome) {
    double result = calcResult(chromosome);

    return fabs(result - desiredResult);
}


int main() {
    int seed = random_device()();
    
    // gen.seed(seed);

    GA<double> ga(
        300, // max generations
        10, 
        2000, // parent size
        0.1, // elite percent
        0.2, // mutation chance
        0.3 // keep worst percent
    );

    ga.initRandFunc = []() -> double {
        return norm_dist(10, 10);
    };

    ga.mutateGeneFunc = [](double value) -> double {
        //  value + norm_dist(0, 0.1);
        return value + norm_dist(10, 10);
    };

    ga.fitnessFunc = fitnessFunc;

    auto start = chrono::high_resolution_clock::now();

    ga.onGenerationEndFunc = [&](auto* population, size_t generation) {
        auto duration = chrono::duration_cast<chrono::milliseconds>(
            chrono::high_resolution_clock::now() - start
        );

        cout << "*-------------------------------------------*" << endl;

        cout << "Generation: " << generation << " Between Calls: " << duration.count() << " ms" << endl;

        vector<double> bestChromosome = population->at(0);

        float bestFitness = fitnessFunc(&bestChromosome);
        float result = calcResult(&bestChromosome);
        
        cout << "*-------------------------------------------*" << endl;
        cout << "Best fitness: " << bestFitness << endl;
        cout << "Desired result: " << "Result" << endl;
        cout << desiredResult << " | " << result << endl;

        cout << "Best chromosome: ";
        printChromosome(&bestChromosome);

        // this_thread::sleep_for(chrono::milliseconds(200));
        start = chrono::high_resolution_clock::now();
    };

    ga.onFinish = [&]() {
        cout << "Generation end" << endl;
        cout << "Seed: " << seed << endl;
    };

    ga.run();


    cout << "Desired result: " << desiredResult << endl;

    printMemoryUsage();

    return 0;
}



