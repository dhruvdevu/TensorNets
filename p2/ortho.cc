#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;
//Only the orthogonality center needs to be normalized
//Whichever term gets the D
int N = 15;
int steps = 1000;
int maxm = 100;
Real cutoff = 0.0;
//Coupling constants
float J1 = -1.0;
float J2 = 0.0;
Real T = 0.01;
SiteSet sites = SpinHalf(N);

void normalizeCheck(std::vector<ITensor> mps, int N) {
    ITensor psi = ITensor(1.0);
    for (int i = 0; i < N; i++) {
        //Print(mps[i]);
        psi = psi*mps[i];
    }
    Real n = norm(psi);//(dag(prime(psi))*psi).real();
    printfln("\nNormalize check = %.10f",n);
}

Real getEnergy(std::vector<ITensor> mps) {
  //TODO(Dhruv): Write this method to efficiently calculate the energy.
  ITensor H, U, D, V;
  Real energy = 0.0;
  for (int i = 0; i < N - 2; i++) {
      H = J1*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 2))*ITensor(sites.op("Id", i + 3))
      + J1*ITensor(sites.op("Sx", i + 1))*ITensor(sites.op("Sx", i + 2))*ITensor(sites.op("Id", i + 3))
      + J1*ITensor(sites.op("Sy", i + 1))*ITensor(sites.op("Sy", i + 2))*ITensor(sites.op("Id", i + 3))
      +J2*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 3))*ITensor(sites.op("Id", i + 2))
      + J2*ITensor(sites.op("Sx", i + 1))*ITensor(sites.op("Sx", i + 3))*ITensor(sites.op("Id", i + 2))
      + J2*ITensor(sites.op("Sy", i + 1))*ITensor(sites.op("Sy", i + 3))*ITensor(sites.op("Id", i + 2));
      ITensor psi = prime(prime(prime(mps[i]*mps[i+1]*mps[i+2], sites(i+1)), sites(i+2)), sites(i + 3));
      Real e = (dag(psi)*H*mps[i]*mps[i+1]*mps[i+3]).real();
      energy += e;
      if (i == 0) {
          U = ITensor(sites(i + 1));
      } else {
          U = ITensor(sites(i + 1), commonIndex(mps[i], mps[i-1]));
      }
      svd(mps[i]*mps[i+1], U, D, V);
      mps[i] = U;
      psi =  D*V;
      U = ITensor(sites(i+2), commonIndex(mps[i], psi));
      svd(psi, U, D, V);
      mps[i+1] = U;
      mps[i+2] = D*V;
  }
  H = J1*ITensor(sites.op("Sz", N))*ITensor(sites.op("Sz", N + 1))
  + J1*ITensor(sites.op("Sx", N))*ITensor(sites.op("Sx", N + 1))
  + J1*ITensor(sites.op("Sy", N))*ITensor(sites.op("Sy", N + 1));
  energy += (dag(prime(prime(mps[N-1]*mps[N-3], sites(N)), sites(N+1)))*H*mps[N-1]*mps[N-2]).real();
  return energy;
}

int main() {
    //DMRG
    auto hMPO = AutoMPO(sites);

    for (int j = 1; j < N - 1; j++) {
    		hMPO += -J1,"Sx",j,"Sx",j+1,"Id",j+2;
            hMPO += -J1,"Sy",j,"Sy",j+1,"Id",j+2;
            hMPO += -J1,"Sz",j,"Sz",j+1,"Id",j+2;
            hMPO += -J2,"Sx",j,"Sx",j+2,"Id",j+1;
            hMPO += -J2,"Sy",j,"Sy",j+2,"Id",j+1;
            hMPO += -J2,"Sz",j,"Sz",j+2,"Id",j+1;

    }

    hMPO += -J1,"Sx",N-1,"Sx",N;
    hMPO += -J1,"Sy",N-1,"Sy",N;
    hMPO += -J1,"Sz",N-1,"Sz",N;
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

    printfln("Orthogonalizing:");
    //Canonical Form - initialize so that orthogonality center is site 1
    ITensor V, D, U(commonIndex(mps[N-1], mps[N-2]));
    svd(mps[N - 1], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
    mps[N - 1] = V;
    mps[N - 2] = mps[N - 2]*U*D;
    mps[N - 2] /= norm(mps[N - 2]);
    printfln("Loop:");
    for (int i = N - 2; i >= 1; i--) {
        //ITensor V = ITensor(sites(i + 1), virtualIndices[i]); //So V gets same indices
        ITensor D, V, U(commonIndex(mps[i], mps[i-1]));
        svd(mps[i], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
        mps[i] = V;
        mps[i - 1] = mps[i - 1]*U*D;
        mps[i - 1] /= norm(mps[i - 1]);
    }
    printfln("TEBD:");
    float progress = 0.0;
    int barWidth = 70;


    for (int st = 0; st < steps; st++) {
        if (st % 100 == 0) {
            Print(st);
        }
        // Forward
        for (int i = 0; i < N - 2; i+=3) {
        	ITensor h = J1*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 2))*ITensor(sites.op("Id", i + 3))
            + J1*ITensor(sites.op("Sx", i + 1))*ITensor(sites.op("Sx", i + 2))*ITensor(sites.op("Id", i + 3))
            + J1*ITensor(sites.op("Sy", i + 1))*ITensor(sites.op("Sy", i + 2))*ITensor(sites.op("Id", i + 3))
            + J2*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Sz", i + 3))*ITensor(sites.op("Id", i + 2))
            + J2*ITensor(sites.op("Sx", i + 1))*ITensor(sites.op("Sx", i + 3))*ITensor(sites.op("Id", i + 2))
            + J2*ITensor(sites.op("Sy", i + 1))*ITensor(sites.op("Sy", i + 3))*ITensor(sites.op("Id", i + 2))
            + J1*ITensor(sites.op("Sz", i + 2))*ITensor(sites.op("Sz", i + 3))*ITensor(sites.op("Id", i))
            + J1*ITensor(sites.op("Sx", i + 2))*ITensor(sites.op("Sx", i + 3))*ITensor(sites.op("Id", i))
            + J1*ITensor(sites.op("Sy", i + 2))*ITensor(sites.op("Sy", i + 3))*ITensor(sites.op("Id", i));
        	auto expTemp = expHermitian(h, -T);
        	//the mps should already have its orthocenter at i
        	auto p1 = mps[i];
        	auto p2 = mps[i + 1];
            auto p3 = mps[i + 2];
        	auto p = p1*p2*p3;
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
            // mps[i] /= norm(mps[i]);
            ITensor temp = D*V;
            U = ITensor(sites(i + 2), commonIndex(mps[i], temp));
            svd(temp, U, D, V);
            mps[i + 1] = U;
            mps[i + 2] = D*V;
            mps[i + 2] /= norm(mps[i + 2]);
            //Set the orthocenter to be the next i value
            if (i < N - 3) {
                ITensor D, V;
                U = ITensor(sites(i + 3), commonIndex(mps[i+2], mps[i + 1]));
                svd(mps[i+2], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
                mps[i + 2] = U;
                // mps[i + 1] /= norm(mps[i + 1]);
                mps[i + 3] = D*V*mps[i + 3];
                mps[i + 3] /= norm(mps[i + 3]);
            }
        }
        //Case 1: no tail terms
        int start = N - 1;
        if (N % 3 == 0) {
            ITensor U = ITensor(commonIndex(mps[N-1], mps[N-2]));
            ITensor D, V;
            svd(mps[N-1], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
            mps[N-2] = mps[N-2]*U*D;
            mps[N-1] = V;
            mps[N-2] /= norm(mps[N-2]);
        } else if (N % 3 == 2) {
            ITensor U = ITensor(commonIndex(mps[N-2], mps[N-3]), sites(N-1));
            ITensor D, V;
            svd(mps[N-2], U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
            mps[N-2] = U;
            mps[N-1] = D*V*mps[N-1];
            mps[N-1] /= norm(mps[N-1]);

            ITensor h = J1*ITensor(sites.op("Sz", N - 1))*ITensor(sites.op("Sz", N))*ITensor(sites.op("Id", N - 2))
            + J1*ITensor(sites.op("Sx", N - 1))*ITensor(sites.op("Sx", N))*ITensor(sites.op("Id", N - 2))
            + J1*ITensor(sites.op("Sy", N - 1))*ITensor(sites.op("Sy", N))*ITensor(sites.op("Id", N - 2))
            + J2*ITensor(sites.op("Sz", N - 2))*ITensor(sites.op("Sz", N))*ITensor(sites.op("Id", N - 3))
            + J2*ITensor(sites.op("Sx", N - 2))*ITensor(sites.op("Sx", N))*ITensor(sites.op("Id", N - 3))
            + J2*ITensor(sites.op("Sy", N - 2))*ITensor(sites.op("Sy", N))*ITensor(sites.op("Id", N - 3));

            auto expTemp = expHermitian(h, -T);
            ITensor p = mps[N-1]*mps[N-2]*mps[N-3];
            p = expTemp*p;
            p = p.noprime();


        }
        
        //Backward
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


printfln("Calculating energy");
printfln("\nGround state energy by TEBD = %.25f", getEnergy(mps));
printfln("Ground state energy by dmrg = %.25f",dmrg_energy);

}
