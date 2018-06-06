#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;

int main() {
    int N = 10;

    //Create ITensor
    std::vector<ITensor> mps(N);
    Index prevRight = Index("virtual1", Link);
    mps[0] = ITensor(Index("physical1", 2, Site), prevRight);
    randomize(mps[0]);
    for (int i = 2; i < N - 1; i++) {
        Index left = Index("virtual2", 2, Link);
        Index physical = Index("physical", 2, Site);
        mps[i] = ITensor(prevRight, physical, left);
        prevRight = left;
        randomize(mps[i]);
    }
    mps[N - 1] = ITensor(prevRight, Index("physical", 2, Site));
    randomize(mps[N-1]);

    //Canonical Form - initialize so that orthogonality center is site 1
    for (int i = N - 1; i > 1; i--) {
        ITensor V = mps[i]; //So V gets same indices
        ITensor D, U;
        svd(mps[i], U, D, V, {"Cutoff", 0});
        mps[i] = V;
        mps[i - 1] = mps[N-2]*U*D;
}



}
