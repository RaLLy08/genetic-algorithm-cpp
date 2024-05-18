#pragma once
#include <stdlib.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <functional>
#include <sstream>

#include <thread>
#include <chrono>
#include "GA.h"
#include "equations.h"

#include <fstream>

namespace EquationsGA {


    using namespace std;

    void runGA();
}

namespace DampedWavesGA {
    using namespace std;

    void runGA();
}