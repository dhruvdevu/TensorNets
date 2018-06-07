#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;

int main() {
    int N = 10;
    int steps = 10;
    maxm = 15;
    cutoff = 0;
    //Coupling constants
    float J = 1.0;
    float h = 1.0;
    Real T = 0.001;
    SiteSet sites = SpinHalf(N);
    OptSet opts;
    opts.add("Cutoff",cutoff);
    opts.add("Maxm",maxm);
    //Create ITensor
    std::vector<ITensor> mps(N);
    Index prevRight = Index("virtual1", Link);
    mps[0] = ITensor(Index("physical1", 2, Site), prevRight);
    randomize(mps[0]);
    for (int i = 2; i < N - 1; i++) {
        //Max bond dimension of 15
        Index left = Index("virtual2", 15, Link);
        //Index physical = Index("physical", 2, Site);
        mps[i] = ITensor(prevRight, sites(i), left);
        prevRight = left;
        randomize(mps[i]);
    }
    mps[N - 1] = ITensor(prevRight, Index("physical", 2, Site));
    randomize(mps[N-1]);

    //Canonical Form - initialize so that orthogonality center is site 1
    for (int i = N - 1; i > 1; i--) {
        ITensor V = mps[i]; //So V gets same indices
        ITensor D, U;
        svd(mps[i], U, D, V, opts);
        mps[i] = V;
        mps[i - 1] = mps[N-2]*U*D;
    }
    //TO DO: Implement TEBD

    for (int st = 0; st < steps; st++) {
        // Odd-even pairs
        for (int i = 0; i < N - 1; i+=2) {
        	ITensor hEven = -h*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 1))
        	- J*0.5*ITensor(sites.op("S+", i))*ITensor(sites.op("Id", i + 1))
        	- J*0.5*ITensor(sites.op("S-", i))*ITensor(sites.op("Id", i + 1));
        	auto expTemp = expHermitian(hEven, -T);
        	//the mps should already have its orthocenter at i
        	auto p1 = mps[i];
        	auto p2 =mps[i + 1];
        	auto p = p1*p2;
        	p = expTemp*p;
        	//p /= norm(p);
        	p = p.noprime();
        	auto U = mps[i];
        	ITensor D, V;
        	svd(p, U, D, V, opts);
            mps[i] = U;
            mps[i + 1] = D*V;
            //TO DO: How to keep it normalized?
            //Set the orthocenter to be the next i value
            if (i < N - 2) {
                U = mps[i + 1];
                svd(mps[i+1], U, D, V, opts);
                mps[i + 1] = U;
                mps[i + 2] = D*V*mps[i + 2];
            }
        }
        //If N is odd, then the orthocenter is N-1, and if N is even, it is N-2 (both of these cases contain N - 1)
        ITensor hlast = -J*0.5*ITensor(sites.op("S+", N)) - J*0.5*ITensor(sites.op("S-", N));
        mps[N - 1] = (mps[N-1]*expHermitian(hlast, -T)).noprime();

        //Even - odd pairs
        int start = (N % 2 == 0) ? N-1 : N-2;
    }
}



}
