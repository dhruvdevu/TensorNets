#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;
//Only the orthogonality center needs to be normalized
//Whichever term gets the D
int N = 19;
int steps = 100;
int maxm = 20;
Real cutoff = 0.0;
//Coupling constants
float J1 = -1.0;
float J2 = 0.0;
Real T = 0.1;
SiteSet sites = SpinHalf(N);

void normalizeCheck(MPS psi) {
    Real n = norm(psi);//(dag(prime(psi))*psi).real();
    printfln("\nNormalize check = %.10f",n);
}

Real getEnergy(MPS psi) {
  //TODO(Dhruv): Write this method to efficiently calculate the energy.
  ITensor H, U, D, V;
  Real energy = 0.0;
  for (int i = 1; i < N - 1; i++) {
      H = J1*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Id", i + 2))
      + J1*ITensor(sites.op("Sx", i))*ITensor(sites.op("Sx", i + 1))*ITensor(sites.op("Id", i + 2))
      + J1*ITensor(sites.op("Sy", i))*ITensor(sites.op("Sy", i + 1))*ITensor(sites.op("Id", i + 2))
      + J2*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 2))*ITensor(sites.op("Id", i + 1))
      + J2*ITensor(sites.op("Sx", i))*ITensor(sites.op("Sx", i + 2))*ITensor(sites.op("Id", i + 1))
      + J2*ITensor(sites.op("Sy", i))*ITensor(sites.op("Sy", i + 2))*ITensor(sites.op("Id", i + 1));
      ITensor p = prime(prime(prime(psi.A(i)*psi.A(i+1)*psi.A(i+2), sites(i)), sites(i + 1)), sites(i + 2));
      Real e = (dag(p)*H*psi.A(i)*psi.A(i+1)*psi.A(i+2)).real();
      energy += e;
      psi.position(i + 1);
  }
  H = J1*ITensor(sites.op("Sz", N))*ITensor(sites.op("Sz", N - 1))
  + J1*ITensor(sites.op("Sx", N))*ITensor(sites.op("Sx", N - 1))
  + J1*ITensor(sites.op("Sy", N))*ITensor(sites.op("Sy", N - 1));
  energy += (dag(prime(prime(psi.A(N)*psi.A(N-1), sites(N)), sites(N-1)))*H*psi.A(N)*psi.A(N-1)).real();
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
    auto psi = MPS(sites);
    psi.position(1);
    printfln("TEBD:");
    for (int st = 0; st < steps; st++) {
        if (st % 100 == 0) {
            Print(st);
        }

        // Forward
        for (int i = 1; i < N - 1; i+=1) {
        	ITensor h = J1*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 1))*ITensor(sites.op("Id", i + 2))
            + J1*ITensor(sites.op("Sx", i))*ITensor(sites.op("Sx", i + 1))*ITensor(sites.op("Id", i + 2))
            + J1*ITensor(sites.op("Sy", i))*ITensor(sites.op("Sy", i + 1))*ITensor(sites.op("Id", i + 2))
            + J2*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 2))*ITensor(sites.op("Id", i + 1))
            + J2*ITensor(sites.op("Sx", i))*ITensor(sites.op("Sx", i + 2))*ITensor(sites.op("Id", i + 1))
            + J2*ITensor(sites.op("Sy", i))*ITensor(sites.op("Sy", i + 2))*ITensor(sites.op("Id", i + 1));
        	auto expTemp = expHermitian(h, -T);
        	//the mps should already have its orthocenter at i
        	auto p1 = psi.A(i);
        	auto p2 = psi.A(i + 1);
            auto p3 = psi.A(i + 2);
        	auto p = p1*p2*p3;
            //PrintData(p);
        	p = expTemp*p;
        	//p /= norm(p);
        	p = p.noprime();
            ITensor U;
            if (i > 1) {
                U = ITensor(sites(i),  commonIndex(p1, psi.A(i - 1)));
            } else {
                U = ITensor(sites(i));
            }
            //U = psi.A(i);
        	ITensor D, V;
            // Print("SVD");
        	svd(p, U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
            // Print("2");
            psi.setA(i, U);
            ITensor temp = D*V;
            U = ITensor(sites(i + 2), commonIndex(psi.A(i), temp));
            // Print("3");
            svd(temp, U, D, V, {"Cutoff", cutoff, "Maxm", maxm});
            psi.setA(i + 1, U);
            psi.setA(i + 2, D*V);
            normalize(psi);
        }
        //Last 2 terms
        ITensor h = J1*ITensor(sites.op("Sz", N - 1))*ITensor(sites.op("Sz", N))
        + J1*ITensor(sites.op("Sx", N - 1))*ITensor(sites.op("Sx", N))
        + J1*ITensor(sites.op("Sy", N - 1))*ITensor(sites.op("Sy", N));
        auto expTemp = expHermitian(h, -T);
        //the mps should already have its orthocenter at i
        auto p1 = psi.A(N);
        auto p2 = psi.A(N - 1);
        auto p = expTemp*p1*p2;

        //auto p = p1*p2;
        p = p.noprime();
        ITensor U = ITensor(sites(N - 1), commonIndex(p, psi.A(N - 2)));//p1;//ITensor(sites(N - 1), commonIndex(mps[N - 2], mps[N - 3]));
        ITensor D, V;
        svd(p, U, D, V, {"Cutoff", cutoff, "Maxm", maxm});

        psi.setA(N - 1, U);
        psi.setA(N, D*V);

        normalize(psi);

        //Take back orthocenter

        psi.position(1);
        //orthogonality center is now the first index (index0, site1)
    }

//Check normalized

Print("normalize");
normalizeCheck(psi);
//Calculate energy


printfln("Calculating energy");
printfln("\nGround state energy by TEBD = %.25f", getEnergy(psi));
printfln("Ground state energy by dmrg = %.25f",dmrg_energy);

}
