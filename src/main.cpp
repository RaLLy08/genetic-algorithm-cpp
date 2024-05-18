#include <iostream>
// #include <sys/resource.h>

#include "examples.h"

using namespace std;

// void printMemoryUsage() {
//     struct rusage usage;
//     getrusage(RUSAGE_SELF, &usage);

//     cout << "Memory usage: " << usage.ru_maxrss << " kb" << endl;
// }


int main() {
    EquationsGA::runGA();
    // auto sineWave1 = printSineWave(30, 1, 0, 0.4);

    // sineWave1();

    // printMemoryUsage();

    return 0;
}



