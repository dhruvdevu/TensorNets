#include "itensor/all.h"

using namespace itensor;
int main() {
int N = 20;
auto sites = SpinHalf(N);
float J = 1.0;
float h = 1.0;

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
auto psi = MPS(sites);

int steps = 1000;

for(int i = 0; i < steps; i++) {
//T : amount of imaginary time to evolve our state by
float T = 0.1;

// Even-odd pairs
for (int i = 2; i < N; i+=2) {
    ITensor hEven = -h*sites.op("Sz", i)*sites.op("Sz", i + 1)
		- J*sites.op("S+", i)*sites.op("Id", i + 1)
		- J*sites.op("S-", i)*sites.op("Id", i + 1);
    auto expTemp = expHermitian(hEven, -T);
    psi.position(i);
    auto p1 = psi.A(i);
    auto p2 = psi.A(i + 1);
	auto p = p1*p2;
	p = expTemp*p;
	p /= norm(p);
	p = p.noprime();
	auto U = psi.A(i);
	ITensor D, V;
	//svd(p, U, D, V);
	svd(p, U, D, V,{"Cutoff",1E-10});
	psi.setA(i, U);
	psi.setA(i + 1, D*V);
	psi.norm();
}


for (int i = 1; i < N; i += 2) {
    ITensor hOdd = -h*sites.op("Sz", i)*sites.op("Sz", i + 1)
		- J*sites.op("S+", i)*sites.op("Id", i + 1)
		- J*sites.op("S-", i)*sites.op("Id", i + 1);
    auto expTemp = expHermitian(hOdd, -T);
    psi.position(i);
    auto p1 = psi.A(i);
    auto p2 = psi.A(i + 1);
	auto p = p1*p2;
	p = expTemp*p;
	p /= norm(p);
	p = p.noprime();
	auto U = psi.A(i);
	ITensor D, V;
	//svd(p, U, D, V);
	svd(p, U, D, V, { "Cutoff",1E-10 });
	psi.setA(i, U);
	psi.setA(i + 1, D*V);
	psi.norm();
}


//Do I need an edge case for the last index in the chain?

/*
psi.position(N);
ITensor hlast = -J*sites.op("S+", N) - J*sites.op("S-", N);
psi.setA(N, (psi.A(N)*hlast).noprime());
*/
psi.orthogonalize();

}
psi.norm();

Real energy = overlap(psi, HMPO, psi);

printfln("Ground state energy = %.20f",energy);

return 0;
}
