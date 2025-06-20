#include <sstream>

namespace molsim {
    static const int MASTER=0;
    class msError;
    class SimulationBox;
    class Molecule;
    class MoleculeType;
    class AtomType;
    class Stretching;
    class Bending;
    class Torsion;
    class U2;
    class U3;
    class U4;
    bool getlines(std::istream&,std::istringstream&);
    inline double round_d(double x);
    void unfold(double**,double**,int,SimulationBox&);
}

class molsim::msError {
    int id;
public:
    msError(int errorId=-1): id(errorId) {}
    ~msError() {}
    int getId() const {return id;}
    char *getMessage() const {return NULL;}
};

class molsim::SimulationBox {
public:
    bool pbc[3];
    int symm,dim; //symm: 1 - cubic, 2 - orthorhombic, 3 - tricilinic; dim: 2 or 3
    double vol,boxl[3],ang[3],height[3],hmat[3][3],hmati[3][3];
    SimulationBox() {}
    ~SimulationBox() {}
};

class molsim::Molecule {
public:
    int inBox,type;
    double com[3];
    double **atom;
    Molecule(): atom(NULL) {}
    ~Molecule() {
        if (atom) delete atom;
    }
};

class molsim::MoleculeType {
    static const int LEN_STR_MOLECULE_NAME=8;
public:
    int type,nAtomType,*atomType;
    char name[LEN_STR_MOLECULE_NAME];
    MoleculeType(): atomType(NULL) {}
    ~MoleculeType() {
        if (atomType) delete atomType;
    }
    int setNumAtomType(int);
    int setName(char mn[]);
};

class molsim::AtomType {
    static const int LEN_STR_ATOM_NAME=8,LEN_STR_ATOM_DESC=32;
public:
    bool hasCharge,hasDisper,Unonbond;
    int atomID;
    double mass,charge;
    char name[LEN_STR_ATOM_NAME],desc[LEN_STR_ATOM_DESC];
    AtomType() {}
    ~AtomType() {}
    int setName(char an[]);
};

class molsim::Stretching {

public:
    Stretching() {}
    ~Stretching() {}
};

class molsim::Bending {

public:
    Bending() {}
    ~Bending() {}
};

class molsim::Torsion {

public:
    Torsion() {}
    ~Torsion() {}
};

class molsim::U2 {

public:
    U2() {}
    ~U2() {}
};

class molsim::U3 {

public:
    U3() {}
    ~U3() {}
};

class molsim::U4 {

public:
    U4() {}
    ~U4() {}
};
