#include "examples.h"

void static printChromosome(vector<double>& chromosome) {
    for (auto it = chromosome.begin(); it != chromosome.end(); ++it) {
        cout << *it << " ";
    }

    cout << endl;
}


namespace EquationsGA {
    double equation(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j) {
        return rosensbrock(a, b) + 
            sphere(c, d) + 
            deJong(e, f) + 
            venter(g, h) + 
            ackley(i, j) + 
            schwefel(a, b) + 
            griewank(c, d) +
            styblinski(e, f);
    }

    double desiredResult = equation(4, 10, 6, 7, 8, 9, 10, 11, 12, 13);

    double calcResult(vector<double>& chromosome) {
        return equation(
            chromosome.at(0),
            chromosome.at(1),
            chromosome.at(2),
            chromosome.at(3),
            chromosome.at(4),
            chromosome.at(5),
            chromosome.at(6),
            chromosome.at(7),
            chromosome.at(8),
            chromosome.at(9)
        );
    }

    double fitnessFunc(vector<double>& chromosome) {
        double result = calcResult(chromosome);

        return fabs(result - desiredResult);
    }

    void runGA() {
        int seed = random_device()();

        Random::seed(seed);
        
        GA ga(
            300, // max generations
            10, 
            2000, // parent size
            0.1, // elite percent
            0.2, // mutation chance
            0.3 // keep worst percent
        );

        ga.initRandFunc = []() -> double {
            return Random::norm_dist(10, 10);
        };

        ga.mutateGeneFunc = [](double value) -> double {
            //  value + norm_dist(0, 0.1);
            return value + Random::norm_dist(10, 10);
        };

        ga.fitnessFunc = fitnessFunc;

        auto start = chrono::high_resolution_clock::now();

        ga.onGenerationEndFunc = [&](auto population, size_t generation) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - start
            );

            cout << "*-------------------------------------------*" << endl;

            cout << "Generation: " << generation << " Between Calls: " << duration.count() << " ms" << endl;

            vector<double> bestChromosome = population.at(0);

            float bestFitness = fitnessFunc(bestChromosome);
            float result = calcResult(bestChromosome);
            
            cout << "*-------------------------------------------*" << endl;
            cout << "Best fitness: " << bestFitness << endl;
            cout << "Desired result: " << "Result" << endl;
            cout << desiredResult << " | " << result << endl;

            cout << "Best chromosome: ";
            printChromosome(bestChromosome);

            // this_thread::sleep_for(chrono::milliseconds(200));
            start = chrono::high_resolution_clock::now();
        };

        ga.onFinish = [&]() {
            cout << "Generation end" << endl;
            cout << "Seed: " << seed << endl;
        };

        ga.run();
    }
}


namespace DampedWavesGA {

    class SineWave {
        public:
            double x = 0;

            void addWave(double amplitude, double frequency, double phase, double decay) {
                waves.push_back(
                    make_tuple(amplitude, frequency, phase, decay)
                );
            }

            vector<
                tuple<double, double, double, double>
            > waves;


            SineWave(double amplitude, double frequency, double phase, double decay = 0) {
                addWave(
                    amplitude, frequency, phase, decay
                );
            }

            // SineWave& operator+(const SineWave& other) {
            //     for (auto wave : other.waves) {
            //         waves.push_back(wave);
            //     }
            //     return *this;
            // }

            double sumAmplitude() {
                double sum = 0;

                for (auto wave : waves) {
                    double amplitude = get<0>(wave);

                    sum += fabs(amplitude);
                }

                return sum;
            }

            double operator()(double x) {
                double y = 0;

                for (auto wave : waves) {
                    double amplitude = get<0>(wave);
                    double frequency = get<1>(wave);
                    double phase = get<2>(wave);
                    double decay = get<3>(wave);

                    y += amplitude * sin(2 * M_PI * frequency * x + phase) * exp(-decay * x);
                }

                return y;
            }
    };
    // double sineWave(double &x, double amplitude, double frequency, double phase, double decay) {
    //     return amplitude * sin(2 * M_PI * frequency * x + phase) * exp(-decay * x);
    // }

    void printWavePoint(int y, int width) {
        string waveStr = "";

        if (y > 0) {
            for (int i = 0; i < width; i++) {
                waveStr += ' ';
            }

            for (int i = 0; i < y; i++) {
                waveStr += 'o';
            }
        } else {
            for (int i = 0; i < width + y; i++) {
                waveStr += ' ';
            }

            for (int i = 0; i < -y; i++) {
                waveStr += 'o';
            }
        }

        cout << waveStr << endl;
    }

    auto printSineWave() {
        SineWave wave = SineWave(10, 2, 20, 0.4);

        wave.addWave(-40, 2, 0, 0.4);
        float maxAmplitude = wave.sumAmplitude();

        // wave.addWave(10, 2, 0, 0.4);
        // wave.addWave(10, 3, 0, 0.4);
        // wave.addWave(20, 4, 2, 0.8);

        double y = 0;
        double x = 0;
            
        while (true) {
            y = wave(x);
            x += 0.03;

            printWavePoint(y, maxAmplitude);

            this_thread::sleep_for(chrono::milliseconds(100));
        }
    }




    void runGA() {
        printSineWave();

        // sineWave1();
    }
}