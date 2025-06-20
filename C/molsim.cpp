#include <string>
#include <cmath>
#include <string.h>
#include "molsim.hpp"

using namespace std;
using namespace molsim;

bool molsim::getlines(istream &fin,istringstream &sline) {
    static string line;
    sline.clear();
    if (getline(fin,line)) {
        sline.str(line);
        return true;
    } else
        return false;
}

inline double molsim::round_d(double x) {
    return floor(x+0.5);
}

void molsim::unfold(double **r,double **ro,int n,SimulationBox &box) {
    int j;
    double d;
    for (int i=0;i<n;i++)
        for (j=0;j<3;j++) {
            d=r[i][j]-ro[i][j];
            r[i][j]=ro[i][j]+d-box.boxl[j]*round_d(d/box.boxl[j]);
        }
}

int molsim::MoleculeType::setNumAtomType(int n) {
    nAtomType=n;
    atomType=new int[nAtomType];
    return 0;
}

int molsim::MoleculeType::setName(char mn[]) {
    strcpy(name,mn);
    return 0;
}

int molsim::AtomType::setName(char an[]) {
    strcpy(name,an);
    return 0;
}
