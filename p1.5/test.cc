#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;

int main() {
    int N = 10;
    std::vector<ITensor> funcs(N);
    for (int i = 0; i < N; i++) {
        funcs[i] = ITensor(Index("name", 2, Site), Index("name", 2, Site));
        randomize(funcs[i]);
        funcs[i] /= norm(funcs[i]);
    }
    for (int j = 0; j < N; j++) {
        printfln("Norm, %.10f, %.10f", j, norm(funcs[j]));
    }
}
