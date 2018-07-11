#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;
//Only the orthogonality center needs to be normalized
//Whichever term gets the D

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
    for(int i = 0; i < I.size(); i++) {
        if (i!=j && i!=j+1) {
            ret = ret*I[i];
        }
    }
    return ret;
}

void normalizeCheck(std::vector<ITensor> mps, int N) {
    int b = 0;
    Real comp = norm(ITensor(1.0));
    for (int i = 0; i < N; i++) {
        if (norm(mps[i]) != comp) {
            printfln("%f, %f", i, norm(mps[i]));
            b++;
        }
    }
    printfln("\nNormalize check=%.10f", b);
}

int main() {
    int N = 10;
    int steps = 1000;
    int maxm = 15;
    int cutoff = 0;
    //Coupling constants
    float J = 1.0;
    float h = 1.0;
    Real T = 0.01;
    SiteSet sites = SpinHalf(N);


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
    auto dmrg_psi = MPS(sites);
    auto sweeps = Sweeps(5);
    sweeps.maxm() = 10,40,100,200,200;
    sweeps.cutoff() = 1E-8;

    auto dmrg_energy = dmrg(dmrg_psi,HMPO,sweeps,{"Quiet",true});
    //                  ^ psi passed by reference,
    //                    can measure properties afterwards

    printfln("Ground state energy by dmrg = %.20f",dmrg_energy);
    //Create ITensor
    std::vector<ITensor> mps(N);
    std::vector<Index> virtualIndices(N - 1);
    virtualIndices[0] = Index("virtual1", maxm, Link);
    mps[0] = ITensor(sites(1), virtualIndices[0]);
    randomize(mps[0]);
    mps[0] /= norm(mps[0]);
    //PrintData(mps[0]);
    for (int i = 1; i < N - 1; i++) {
        //Max bond dimension of 15
        virtualIndices[i] = Index("virtual2", maxm, Link);
        //Index physical = Index("physical", 2, Site);
        mps[i] = ITensor(virtualIndices[i - 1], sites(i + 1), virtualIndices[i]);
        randomize(mps[i]);
        mps[i] /= norm(mps[i]);
    }
    mps[N - 1] = ITensor(virtualIndices[N - 2], sites(N));
    randomize(mps[N-1]);
    mps[N-1] /= norm(mps[N-1]);
    normalizeCheck(mps, N);
    //PrintData(mps[N-1]);
    printfln("Orthogonalizing:");
    //Canonical Form - initialize so that orthogonality center is site 1
    ITensor V, D, U(commonIndex(mps[N-1], mps[N-2]));
    /**ITensor V = ITensor(sites(N));
    //PrintData(V);
    ITensor D;
    ITensor U = ITensor(virtualIndices[N - 2]); **/
    svd(mps[N - 1], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
    mps[N - 1] = V;
    mps[N - 2] = mps[N - 2]*U*D;
    mps[N - 2] /= norm(mps[N - 2]);
    printfln("Loop:");
    for (int i = N - 2; i > 1; i--) {
        //ITensor V = ITensor(sites(i + 1), virtualIndices[i]); //So V gets same indices
        ITensor D, V, U(commonIndex(mps[i], mps[i-1]));
        svd(mps[i], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
        mps[i] = V;
        mps[i - 1] = mps[i - 1]*U*D;
        mps[i - 1] /= norm(mps[i - 1]);
    }
    for (int j = 0; j < N; j++) {
        printfln("j=%d, rank=%d", j, rank(mps[j]));
    }
    printfln("TEBD:");

    for (int st = 0; st < steps; st++) {
        // Odd-even pairs (odd site first index)
        for (int i = 0; i < N - 1; i+=2) {
        	ITensor hEven = -h*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 2))
        	- J*0.5*ITensor(sites.op("S+", i + 1))*ITensor(sites.op("Id", i + 2))
        	- J*0.5*ITensor(sites.op("S-", i + 1))*ITensor(sites.op("Id", i + 2));
        	auto expTemp = expHermitian(hEven, -T);
        	//the mps should already have its orthocenter at i
        	auto p1 = mps[i];
        	auto p2 = mps[i + 1];
            //PrintData(p1);
            //PrintData(p2);
        	auto p = p1*p2;
            //PrintData(p);
        	p = expTemp*p;
        	//p /= norm(p);
        	p = p.noprime();
            ITensor U;
            if (i > 0) {
                U = ITensor(sites(i + 1),  commonIndex(mps[i], mps[i-1]));
            } else {
                U = ITensor(sites(i + 1));//virtualIndices[i]);
            }
        	ITensor D, V;
        	svd(p, U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
            mps[i] = U;
            mps[i] /= norm(mps[i]);
            mps[i + 1] = D*V;
            mps[i+1] /= norm(mps[i+1]);
            //Set the orthocenter to be the next i value
            if (i < N - 2) {
                ITensor D, V;
                U = ITensor(sites(i + 2), commonIndex(mps[i+1], mps[i]));//virtualIndices[i + 1]);
                svd(mps[i+1], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
                mps[i + 1] = U;
                mps[i + 1] /= norm(mps[i + 1]);
                mps[i + 2] = D*V*mps[i + 2];
                mps[i + 2] /= norm(mps[i + 2]);
            }
        }
        //If N is odd, then the orthocenter is N-1, and if N is even, it is N-2 (both of these cases contain N - 1)
        ITensor hlast = -J*0.5*ITensor(sites.op("S+", N)) - J*0.5*ITensor(sites.op("S-", N));
        mps[N - 1] = (mps[N-1]*expHermitian(hlast, -T)).noprime();
        mps[N - 1] /= norm(mps[N -1]);
        //Even - odd pairs
        int start = (N % 2 == 0) ? N-2 : N-1;
        for (int i = start; i > 1; i-=2) {
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
        	ITensor D, V, U(sites(i), commonIndex(mps[i-1], mps[i-2]));//, commonIndex(mps[i], mps[i-1]));//virtualIndices[i - 1]);
        	svd(p, U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
            mps[i] = V;
            //mps[i] /= norm(mps[i]);
            mps[i - 1] = U*D;
            mps[i - 1] /= norm(mps[i - 1]);
            //TO DO: How to keep it normalized?
            //Set the orthocenter to be the next i value
            U = ITensor(commonIndex(mps[i-1], mps[i - 2]));//virtualIndices[i - 2]);
            svd(mps[i - 1], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
            mps[i - 2] = mps[i - 2]*U*D;
            mps[i - 2] /= norm(mps[i - 2]);
            mps[i - 1] = V;
            //mps[i - 1] /= norm(mps[i - 1]);
        }
        //orthogonality center is now the first index (index0, site1)
    }

//Check normalized
normalizeCheck(mps, N);
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
printfln("Constructing H");for (int i = 0; i < N-1; i++) {
    H += -h*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 2))*makeOtherIndicesIdentity2(Id, i)
    - J*0.5*ITensor(sites.op("S+", i + 1))*makeOtherIndicesIdentity1(Id, i)
    - J*0.5*ITensor(sites.op("S-", i + 1))*makeOtherIndicesIdentity1(Id, i);
}
H += -J*0.5*ITensor(sites.op("S+", N))*makeOtherIndicesIdentity1(Id, N - 1) - J*0.5*ITensor(sites.op("S-", N))*makeOtherIndicesIdentity1(Id, N - 1);
printfln("Calculating energy");

for (int i = 0; i < N; i++) {

}
Real energy = (dag(prime(psi))*H*psi).real();
printfln("\nGround state energy by TEBD= %.10f",energy);
}
