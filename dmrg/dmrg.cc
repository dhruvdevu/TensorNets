#include "itensor/all.h"

using namespace itensor;
int main() {
int N = 10;
auto sites = SpinHalf(N);
float J = 1.0;
float h = 1.0;
auto ampo = AutoMPO(sites);
for (int j = 1; j < N; j++) {
    ampo += -J*0.5,"S+",j;
    ampo += -J*0.5,"S-",j;
    ampo += -h,"Sz",j,"Sz",j+1;
}
ampo += -J*0.5,"S+",N;
ampo += -J*0.5,"S-",N;
auto H = MPO(ampo);
auto psi = MPS(sites);
auto sweeps = Sweeps(5);
sweeps.maxm() = 10,40,100,200,200;
sweeps.cutoff() = 0;

auto energy = dmrg(psi,H,sweeps,{"Quiet",false});
//                  ^ psi passed by reference,
//                    can measure properties afterwards

printfln("Ground state energy = %.20f",energy);

return 0;
}
