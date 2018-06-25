#include "itensor/all.h"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace itensor;

void orthocenter(ITensor* psi, int pos) {
	
}

int main() {
//Number of sites
int N = 10;
//Coupling constants
float J = 1.0;
float h = 1.0;
//Constructing the Heisenberg Hamiltonian MPO
auto hMPO = AutoMPO(sites);

for (int j = 1; j < N; j++) {
		hMPO += -J*0.5,"S+",j;
		hMPO += -J*0.5,"S-",j;
		hMPO += -h,"Sz",j,"Sz",j+1;
}

hMPO += -J*0.5, "S+", N;
hMPO += -J*0.5, "S-", N;
auto sites = SpinHalf(N);
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


//Create MPS
std::vector <Index> inds(N);
for (int i = 1; i <= N; i++) {
	inds[i] = Index("name", 2, Site);
}

ITensor psi = ITensor(inds);

int steps = 50;
std::ofstream file;
file.open("10N50steps1e-3T.out");
for(int i = 0; i < steps; i++) {
//T : amount of imaginary time to evolve our state by
Real T = 0.001;

// Even-odd pairs
for (int i = 2; i < N; i+=2) {
	ITensor hEven = -h*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 1))
	- J*0.5*ITensor(sites.op("S+", i))*ITensor(sites.op("Id", i + 1))
	- J*0.5*ITensor(sites.op("S-", i))*ITensor(sites.op("Id", i + 1));
	auto expTemp = expHermitian(hEven, -T);
	//psi.position(i);
	orthocenter(psi, i)
	psi.norm();
	auto p1 = psi.A(i);
	auto p2 = psi.A(i + 1);
	auto p = p1*p2;
	p = expTemp*p;
	p /= norm(p);
	p = p.noprime();
	auto U = psi.A(i);
	ITensor D, V;
	svd(p, U, D, V, {"Cutoff", 0});
	//svd(p, U, D, V,{"Cutoff",1E-10});
	psi.setA(i, U);
	psi.setA(i + 1, D*V);
}

//psi.position(N);
orthocenter(psi, N);
psi.norm();
ITensor hlast = -J*0.5*ITensor(sites.op("S+", N)) - J*0.5*ITensor(sites.op("S-", N));
psi.setA(N, (psi.A(N)*expHermitian(hlast, -T)).noprime());
for (int i = 1; i < N; i += 2) {
		ITensor hOdd = -h*ITensor(sites.op("Sz", i))*ITensor(sites.op("Sz", i + 1))
		- J*0.5*ITensor(sites.op("S+", i))*ITensor(sites.op("Id", i + 1))
		- J*0.5*ITensor(sites.op("S-", i))*ITensor(sites.op("Id", i + 1));

	auto expTemp = expHermitian(hOdd, -T);
	//psi.position(i);
	//orthogonalize by hand:
	orthocenter(psi, i);
	psi.norm();
	auto p1 = psi.A(i);
	auto p2 = psi.A(i + 1);
	auto p = p1*p2;
	p = expTemp*p;
	p /= norm(p);
	p = p.noprime();
	auto U = psi.A(i);
	ITensor D, V;
	svd(p, U, D, V, {"Cutoff", 0});
	//svd(p, U, D, V, { "Cutoff",1E-10 });
	psi.setA(i, U);
	psi.setA(i + 1, D*V);
}


//Do I need an edge case for the last index in the chain?



psi.orthogonalize({"Cutoff", 0, "Maxm", 10});
psi.norm();

Real energy = overlap(psi, HMPO, psi);

if (file) {
	file << std::fixed << std::setprecision(10) << energy;
	file << "\n";
}
//printfln("Current Energy = %.20f",energy);
//printfln("Max bond dimernsion = %.20f", maxM(psi));

}

file.close();
Real energy = overlap(psi, HMPO, psi);

printfln("Ground state energy = %.20f",energy);
printfln("Max bond dimernsion = %.20f", maxM(psi));
return 0;
}
