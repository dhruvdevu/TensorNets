#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;

ITensor makeOtherIndicesIdentity1(std::vector<ITensor> I, int j) {
    ITensor ret = ITensor(1.0);
    for(int i = 0; i < I.size(); i++) {
        if (i != j) {
            ret = ret*I[i];
        }
    }
    return ret;
}

ITensor makeOtherIndicesIdentity2(std::vector<ITensor> I, int j) {
    ITensor ret = ITensor(1.0);
    for(int i =0; i < I.size(); i++) {
        if (i!=j && i!=j+1) {
            ret = ret*I[i];
        }
    }
    return ret;
}

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

    //DMRG
    auto hMPO = AutoMPO(sites);

    for (int j = 1; j < N; j++) {
    		hMPO += -J*0.5,"S+",j;
    		hMPO += -J*0.5,"S-",j;
    		hMPO += -h,"Sz",j,"Sz",j+1;
    }

    hMPO += -J*0.5, "S+", N;
    hMPO += -J*0.5, "S-", N;
    //HMPO is to calculate the energy conveniently
    auto HMPO = MPO(hMPO);
    auto dmrg_psi = MPS(sites)
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,40,100,200,200;
    sweeps.cutoff() = 1E-8;

    auto dmrg_energy = dmrg(dmrg_psi,HMPO,sweeps,{"Quiet",true});
    //                  ^ psi passed by reference,
    //                    can measure properties afterwards

    printfln("Ground state energy by dmrg = %.20f",dmrg_energy);
    //Create ITensor
    std::vector<ITensor> mps(N);
    Index prevRight = Index("virtual1", Link);
    mps[0] = ITensor(sites(1), prevRight);
    randomize(mps[0]);
    for (int i = 1; i < N - 1; i++) {
        //Max bond dimension of 15
        Index left = Index("virtual2", 15, Link);
        //Index physical = Index("physical", 2, Site);
        mps[i] = ITensor(prevRight, sites(i + 1), left);
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
        // Odd-even pairs (odd site first index)
        for (int i = 0; i < N - 1; i+=2) {
        	ITensor hEven = -h*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 2))
        	- J*0.5*ITensor(sites.op("S+", i + 1))*ITensor(sites.op("Id", i + 2))
        	- J*0.5*ITensor(sites.op("S-", i + 1))*ITensor(sites.op("Id", i + 2));
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
        int start = (N % 2 == 0) ? N-2 : N-1;
        for (int j = start; i > 1; i-=2) {
        	ITensor hOdd = -h*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 1))
        	- J*0.5*ITensor(sites.op("S+", i))*ITensor(sites.op("Id", i + 1))
        	- J*0.5*ITensor(sites.op("S-", i))*ITensor(sites.op("Id", i + 1));
        	auto expTemp = expHermitian(hOdd, -T);
        	//the mps should already have its orthocenter at i
        	auto p1 = mps[i - 1];
        	auto p2 =mps[i];
        	auto p = p1*p2;
        	p = expTemp*p;
        	//p /= norm(p);
        	p = p.noprime();
        	auto V = mps[i];
        	ITensor D, U;
        	svd(p, U, D, V, opts);
            mps[i] = V;
            mps[i - 1] = U*D;
            //TO DO: How to keep it normalized?
            //Set the orthocenter to be the next i value
            V = mps[i - 1];
            svd(mps[i - 1], U, D, V, opts);
            mps[i - 2] = mps[i - 1]*U*D;
            mps[i - 1] = V;
        }
        //orthogonality center is now the first index (index0, site1)
    }

//Calculate energy
ITensor psi = ITensor(1.0);
for (int i = 0; i < N; i++) {
    psi = psi*mps[i];
}
std::vector<ITensor> Id(N);
for(int i = 0; i < N; i++) {
    Id[i] = ITensor(sites.op("Id", i + 1));
}
ITensor H;
for (int i = 0; i < N; i++) {
    H += -h*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 2))*makeOtherIndicesIdentity2(Id, i)
    - J*0.5*ITensor(sites.op("S+", i + 1))*makeOtherIndicesIdentity1(Id, i)
    - J*0.5*ITensor(sites.op("S-", i + 1))*makeOtherIndicesIdentity1(Id, i);
}
H += -J*0.5*ITensor(sites.op("S+", N))*makeOtherIndicesIdentity1(N - 1) - J*0.5*ITensor(sites.op("S-", N))*makeOtherIndicesIdentity1(N - 1);


Real energy = (dag(prime(psi))*H*psi).real();
printfln("\nGround state energy by TEBD= %.10f",energy);
}
