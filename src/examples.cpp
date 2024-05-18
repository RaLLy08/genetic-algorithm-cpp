#include "examples.h"

string static chromosomeToString(vector<double>& chromosome) {
    string str = "";

    for (auto it = chromosome.begin(); it != chromosome.end(); ++it) {
        str += to_string(*it) + " ";
    }

    return str;
}


namespace EquationsGA {
    double equation(double x, double y) {
        return ackley(x, y);
    }

    double desiredResult = 0;

    double calcResult(vector<double>& chromosome) {
        return equation(
            chromosome.at(0),
            chromosome.at(1)
        );
    }

    double fitnessFunc(vector<double>& chromosome) {
        double result = calcResult(chromosome);

        return fabs(result - desiredResult);
    }

    void writeLog(ofstream& logFile, vector<vector<double>>& population, int generation) {
        string separator = ";";

        if (generation == 0) {
            string header = "Generation" + separator + "x" + separator + "y" + separator + "z";
            logFile << header << endl;
        }

        // vector<double> bestChromosome = population.at(0);
        // float x = bestChromosome.at(0);
        // float y = bestChromosome.at(1);
        // float z = equation(x, y);
        // logFile << generation << separator << x << separator << y << separator << z << endl;
        
        for (auto chromosome : population) {
            float x = chromosome.at(0);
            float y = chromosome.at(1);
            float z = equation(x, y);

            logFile << generation << separator << x << separator << y << separator << z << endl;
        }
        
    }

    void runGA() {
        int seed = random_device()();

        Random::seed(-1497870209);
        
        ofstream logFile;

        logFile.open("log.csv", std::ios::out | std::ios::trunc);

        GA ga(
            150, // max generations
            2, // chromosome size
            10, // parent size
            0.1, // elite percent
            0.4, // mutation chance
            0.3// keep worst percent
        );

        ga.initRandFunc = []() -> double {
            return Random::norm_dist(0, 10);
        };

        ga.mutateGeneFunc = [](double value) -> double {
            //  value + norm_dist(0, 0.1);
            return value + Random::norm_dist(0.5, 2);
        };

        ga.fitnessFunc = fitnessFunc;


        auto start = chrono::high_resolution_clock::now();

        ga.onGenerationEndFunc = [&](auto population, size_t generation) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(
                chrono::high_resolution_clock::now() - start
            );

            cout << endl << "*-------------------------------------------*" << endl;

            cout << "Generation: " << generation << " Between Calls: " << duration.count() << " ms" << endl;

            vector<double> bestChromosome = population.at(0);

            float bestFitness = fitnessFunc(bestChromosome);
            float result = calcResult(bestChromosome);
            
            cout << "*-------------------------------------------*" << endl;
            cout << "Best fitness: " << bestFitness << endl;
            cout << "Desired result: " << "Result" << endl;
            cout << desiredResult << " | " << result << endl;

            cout << "Best chromosome: " << endl;
            cout << chromosomeToString(bestChromosome);

            writeLog(logFile, population, generation);

            // this_thread::sleep_for(chrono::milliseconds(20));
            start = chrono::high_resolution_clock::now();
        };

        ga.onFinish = [&]() {
            cout << "Generation end" << endl;
            cout << "Seed: " << seed << endl;
            logFile.close();
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