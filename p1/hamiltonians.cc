#include "itensor/all.h"
#include "string"
#include "vector"
using namespace itensor;


ITensor
makeOtherIndices1(std::vector<ITensor> I, int j) {
    ITensor ret = ITensor(1.0);
    for(int i = 0; i < I.size(); i++) {
        if (i != j) {
            ret = ret*I[i];
        }
    }
    return ret;
}

ITensor
makeOtherIndices2(std::vector<ITensor> I, int j) {
    ITensor ret = ITensor(1.0);
    for(int i =0; i < I.size(); i++) {
        if (i!=j && i!=j+1) {
            ret = ret*I[i];
        }
    }
    return ret;
}

ITensor
makeSz(Index const& s)
    {
    auto Sz = ITensor(s,prime(s));
    Sz.set(s(1),prime(s)(1), 0.5);
    Sz.set(s(2),prime(s)(2),-0.5);
    return Sz;
    }

ITensor
makeSx(Index const& s)
    {
    auto Sx = ITensor(s,prime(s));
    Sx.set(s(1),prime(s)(2), 0.5);
    Sx.set(s(2),prime(s)(1), 0.5);
    return Sx;
    }

ITensor
makeId(Index const& s)
    {
    auto Id = ITensor(s,prime(s));
    Id.set(s(1),prime(s)(1),1);
    Id.set(s(2),prime(s)(2),1);
    return Id;
    }

int main() {
auto s1 = Index("s1",2, Site);
auto s2 = Index("s2", 2, Site);
auto s3 = Index("s3", 2, Site);
Real J = 1;
Real h = 1;
Real beta = 15;
ITensor Sx1 = makeSx(s1);
ITensor Sx2 = makeSx(s2);
ITensor Sx3 = makeSx(s3);
ITensor Sz1 = makeSz(s1);
ITensor Sz2 = makeSz(s2);
ITensor Sz3 = makeSz(s3);
ITensor Id1 = makeId(s1);
ITensor Id2 = makeId(s2);
ITensor Id3 = makeId(s3);

ITensor H = -J*(Sz1*Sz2*Id3 + Id1*Sz2*Sz3) - h*(Sx1*Id2*Id3 + Id1*Sx2*Id3 + Id1*Id2*Sx3);
PrintData(H);
ITensor U;
ITensor D;
diagHermitian(H, U, D);
PrintData(D);

auto expH = expHermitian(H, -beta);
auto init_state = random(ITensor(s1, s2, s3));
init_state /= norm(init_state);
Real initEn = (dag(prime(init_state))*H*init_state).real();
printfln("\nInitial Energy expectation = %.10f", initEn);
auto fin_state = expH*init_state;
fin_state.noprime();
fin_state /= norm(fin_state);
Real finEn = (dag(prime(fin_state))*H*fin_state).real();
printfln("At beta = %.1f, energy = %.10f", beta, finEn);

// for an arbitrary number of sites:

size_t n = 10;
std::vector<Index> s(n);
std::vector<ITensor> Sx(n);
std::vector<ITensor> Sz(n);
std::vector<ITensor> Id(n);
for(int i = 0; i < n; i++) {
    s[i] = Index("name", 2, Site);
    printfln("Created index %.1f", i);
    Sx[i] = makeSx(s[i]);
    Sz[i] = makeSz(s[i]);
    Id[i] = makeId(s[i]);
}
printfln("\nFinished init\n");
ITensor nH;
for(int j = 0; j < n - 1; j++) {
    nH += -J*(Sz[j]*Sz[j+1]*makeOtherIndices2(Id, j)) - h*Sx[j]*makeOtherIndices1(Id, j);
}
nH += -h*Sx[n-1]*makeOtherIndices1(Id, n - 1);
printfln("\nConstructed Hamiltonian\n");
diagHermitian(nH, U, D);
PrintData(D);

auto psi = ITensor(s);
randomize(psi);
psi /= norm(psi);
expH = expHermitian(nH, -beta);
initEn = (dag(prime(psi))*nH*psi).real();
printfln("\nInitial Energy expectation = %.10f",initEn);
auto phi = expH*psi;
phi /= norm(phi);
phi.noprime();
finEn = (dag(prime(phi))*nH*phi).real();
printfln("\nAt beta = %.1f, energy = %.10f", beta, finEn); 
return 0;
}



