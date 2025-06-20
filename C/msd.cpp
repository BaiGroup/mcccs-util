/** Calculate time-dependent properties
//calculates block averages using the method of Flyvbjerg and Petersen JCP 91 (1989) pg 461
    @todo
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include "molsim.hpp"

using namespace std;
using namespace molsim;

class PostProcess {
    SimulationBox box;
    AtomType *atomtyp;
    MoleculeType *moltyp;
    Molecule *mol,*mol_unfolded;
    bool initialized;
    int nBlk,nElem,nAtmType,nMolecType,nAtom,nMol,nDataEntry,iFrame,**blkLen,**ref,***count;
    double dt,****avg,****data;
#ifdef __MICROPORE__
    int n_meso,n_total;
#endif

    int cIndex(int i) {
        return (nAtom+i);
    }

    int cIndex(int i,int j) {
        return ( nMolecType + (nMolecType+(nMolecType-i+1))*i/2 + j-i );
    }

    void output(char fname[],int* count[],double **avg[]) {
        ofstream fou(fname,ofstream::out);
        int iBlk=0,n;
        double msd;
        fou << "# time[ps] MSD[Ang^2] x y z\n";
        do {
            n=pow(double(nElem),iBlk);
            if (n>iFrame || iBlk>=nBlk) break;

            for (int iElem=0;iElem<blkLen[iBlk][nAtom];iElem++) {
                if (count[iBlk][iElem]>0.0) {
                    fou << dt*(iElem+1)*n;
                    msd=0.0;
                    for (int j=0;j<3;j++) {
                        avg[iBlk][iElem][j]/=double(count[iBlk][iElem]);
                        msd+=avg[iBlk][iElem][j];
                    }
                    fou << " " << msd << " " << avg[iBlk][iElem][0] << " "
                        << avg[iBlk][iElem][1] << " " << avg[iBlk][iElem][2] << "\n";
                }
            }

        } while (++iBlk);
        fou.close();

    }

public:
#ifdef __LAMMPS__
    map<int,int> atomMap;
#else
    map<string,int> atomMap;
#endif
    map<string,int> moleculeMap;

    PostProcess(int nBlock,int nElement,double deltaT,int nAtomType,int nMoleculeType,int nA): nBlk(nBlock),nElem(nElement),dt(deltaT),nAtmType(nAtomType),nMolecType(nMoleculeType),nAtom(nA) {
        initialized=false;
        nDataEntry=nMolecType*(nMolecType+3)/2;
        iFrame=0;
        box.dim=3;
        box.symm=2;
        for (int i=0;i<3;i++) box.pbc[i]=1;
        atomtyp=new AtomType[nAtmType];
        moltyp=new MoleculeType[nMolecType];
        mol=new Molecule[nAtom+nMolecType];
        mol_unfolded=new Molecule[nAtom];
        blkLen=new int*[nBlk];
        ref=new int*[nBlk];
        data=new double***[nBlk];
        for (int i=0;i<nBlk;i++) {
            blkLen[i]=new int[nAtom+1];
            ref[i]=new int[nAtom+1];
            for (int j=0;j<=nAtom;j++) {
                blkLen[i][j]=0;
                ref[i][j]=-1;
            }
            data[i]=new double**[nElem];
            for (int j=0;j<nElem;j++) {
                data[i][j]=new double*[nAtom+nMolecType];
                for (int k=0;k<nAtom+nMolecType;k++)
                    data[i][j][k]=new double[3];
            }
        }

        count=new int**[nMolecType*(nMolecType+1)];
        avg=new double***[nMolecType*(nMolecType+1)];
        for (int k=0;k<nDataEntry;k++) {
            count[k]=new int*[nBlk];
            avg[k]=new double**[nBlk];
            for (int i=0;i<nBlk;i++) {
                count[k][i]=new int[nElem];
                avg[k][i]=new double*[nElem];
                for (int j=0;j<nElem;j++) {
                    avg[k][i][j]=new double[3];
                    count[k][i][j]=0;
                    avg[k][i][j][0]=0.0;
                    avg[k][i][j][1]=0.0;
                    avg[k][i][j][2]=0.0;
                }
            }
        }

#ifdef __MICROPORE__
        n_meso=0;
        n_total=0;
#endif
    }

    ~PostProcess() {
        if (mol) delete [] mol;
        if (mol_unfolded) delete [] mol_unfolded;
        if (moltyp) delete [] moltyp;
        if (atomtyp) delete [] atomtyp;

        for (int i=0;i<nBlk;i++) {
            for (int j=0;j<nElem;j++) {
                for (int k=0;k<nAtom+nMolecType;k++) {
                    delete [] data[i][j][k];
                }
                delete [] data[i][j];
            }
            delete [] blkLen[i];
            delete [] ref[i];
            delete [] data[i];
        }

        for (int k=0;k<nDataEntry;k++) {
            for (int i=0;i<nBlk;i++) {
                for (int j=0;j<nElem;j++) {
                    delete [] avg[k][i][j];
                }
                delete [] count[k][i];
                delete [] avg[k][i];
            }
            delete [] count[k];
            delete [] avg[k];
        }

        delete [] blkLen;
        delete [] ref;
        delete [] count;
        delete [] avg;
        delete [] data;
    }

    void setAtomType(int iAtomType,char name[],double mass) {
        if (iAtomType>=nAtmType) throw msError(-1);
        atomtyp[iAtomType].setName(name);
        atomtyp[iAtomType].mass=mass;
    }

    void setMoleculeType(int iMoleculeType,char name[],vector<int> atoms) {
        if (iMoleculeType>=nMolecType) throw msError(-1);
        moltyp[iMoleculeType].setName(name);
        moltyp[iMoleculeType].setNumAtomType(atoms.size());
        for (int i=0;i<moltyp[iMoleculeType].nAtomType;i++)
            moltyp[iMoleculeType].atomType[i]=atoms[i];
    }

#ifdef __MICROPORE__
    bool isInMicropore(double com[]) {
        double r[3];
        for (int i=0; i<3; i++) {
            r[i] = com[i] - floor(com[i]/box.boxl[i]) * box.boxl[i];
        }
        if (((r[0]>20.022)&&(r[0]<40.044)) || ((r[1]>14.9)&&(r[1]<34.8))) {
            return true;
        } else {
            //cout << com[0] << ", " << com[1] << ", " << com[2] << "; folded: " << r[0] << ", " << r[1] << ", " << r[2];
            return false;
        }
    }
#endif

    void proc(const char *fname) {
#ifdef __LAMMPS__
#include "msd-pp-proc-lammps.cpp"
#else
#include "msd-pp-proc-gromacs.cpp"
#endif
    }

    void accumulateMsd(int iFrame) {
        int iBlk=0,iMolecType,jMolecType,n,idx,nMolPerType[nMolecType];
#ifdef __MICROPORE__
        bool isInside[nMol];
#endif

        // Calculate COM of all the molecules of the same type
        for (iMolecType=0;iMolecType<nMolecType;iMolecType++) {
            for (int j=0;j<3;j++) mol[cIndex(iMolecType)].com[j]=0.0;
            nMolPerType[iMolecType]=0;
        }

        for (int iMol=0;iMol<nMol;iMol++) {
#ifdef __MICROPORE__
            isInside[iMol] = isInMicropore(mol[iMol].com);
            n_total++;

            if (!isInside[iMol]) {
                n_meso++;
                //cout << " at frame " << iFrame << ", molecule " << iMol << endl;
                for (int i=0;i<nBlk;i++) {
                    blkLen[i][iMol]=0;
                    ref[i][iMol]=-1;
                }
                continue;
            }
#endif
            iMolecType=mol[iMol].type;
            for (int j=0;j<3;j++) {
                mol[cIndex(iMolecType)].com[j]+=mol[iMol].com[j];
            }
            nMolPerType[iMolecType]++;
        }
        return;

        for (iMolecType=0;iMolecType<nMolecType;iMolecType++) {
            for (int j=0;j<3;j++) mol[cIndex(iMolecType)].com[j]/=nMolPerType[iMolecType];
        }

        do {
            n=pow(double(nElem),iBlk);
            if (n>iFrame || iFrame%n) break;

            for (int iMol=0;iMol<nMol;iMol++) {
#ifdef __MICROPORE__
                if (!isInside[iMol]) continue;
#endif
                ref[iBlk][iMol]=(ref[iBlk][iMol]+1)%nElem;

                iMolecType=mol[iMol].type;
                for (int i=0;i<blkLen[iBlk][iMol];i++) {
                    if (i==0 && iBlk>0) continue;
                    count[iMolecType][iBlk][i]++;
                    for (int j=0;j<3;j++) {
                        avg[iMolecType][iBlk][i][j]+=pow(mol[iMol].com[j]-data[iBlk][(ref[iBlk][iMol]-i-1+nElem)%nElem][iMol][j],2);
                    }
                }

                for (int j=0;j<3;j++) {
                    data[iBlk][ref[iBlk][iMol]][iMol][j]=mol[iMol].com[j];
                }

                if (blkLen[iBlk][iMol]<nElem) blkLen[iBlk][iMol]++;
            }

            ref[iBlk][nAtom]=(ref[iBlk][nAtom]+1)%nElem;

            for (int iMolecType=0;iMolecType<nMolecType;iMolecType++) {
                for (int jMolecType=iMolecType;jMolecType<nMolecType;jMolecType++) {
                    idx=cIndex(iMolecType,jMolecType);
                    for (int i=0;i<blkLen[iBlk][nAtom];i++) {
                        if (i==0 && iBlk>0) continue;
                        count[idx][iBlk][i]++;
                        for (int j=0;j<3;j++) {
                            avg[idx][iBlk][i][j]+=(mol[cIndex(iMolecType)].com[j]-data[iBlk][(ref[iBlk][nAtom]-i-1+nElem)%nElem][cIndex(iMolecType)][j])*(mol[cIndex(jMolecType)].com[j]-data[iBlk][(ref[iBlk][nAtom]-i-1+nElem)%nElem][cIndex(jMolecType)][j]);
                        }
                    }
                }
            }

            for (int iMolecType=0;iMolecType<nMolecType;iMolecType++) {
                for (int j=0;j<3;j++) {
                    data[iBlk][ref[iBlk][nAtom]][cIndex(iMolecType)][j]=mol[cIndex(iMolecType)].com[j];
                }
            }

            if (blkLen[iBlk][nAtom]<nElem) blkLen[iBlk][nAtom]++;
        } while (++iBlk);
    }

    void output(char fNamePrefix[]){
        char fname[128];

#ifdef __MICROPORE__
        cout << float(n_meso)/n_total;
        return;
#endif
        for (int iMolecType=0;iMolecType<nMolecType;iMolecType++) {
            sprintf(fname,"%s-self%2.2d.xvg",fNamePrefix,iMolecType+1);
            output(fname,count[iMolecType],avg[iMolecType]);
        }

        for (int iMolecType=0;iMolecType<nMolecType;iMolecType++) {
            for (int jMolecType=iMolecType;jMolecType<nMolecType;jMolecType++) {
                sprintf(fname,"%s-corrected%2.2d-%2.2d.xvg",fNamePrefix,iMolecType+1,jMolecType+1);
                output(fname,count[cIndex(iMolecType,jMolecType)],avg[cIndex(iMolecType,jMolecType)]);
            }
        }
   }
};

int main(int argc,char *argv[]) {
    ifstream fin(argv[1]);
    istringstream sline;
    string fname;
    int nBlock,nElement,nAtomType,nMoleculeType,nAtom,type;
    double deltaT,mass;
    char name[8];
    vector<int> atoms;

    //cin.sync_with_stdio(false);
    //cout.sync_with_stdio(false);

    getlines(fin,sline); // number of blocks,number of block elements, time [ps] between two frames, number of atom types, number of molecule types
    getlines(fin,sline);
    sline >> nBlock >> nElement >> deltaT >> nAtomType >> nMoleculeType;
    nAtom=atoi(argv[3]);
    cout << "# MSD for " << nAtom << " atoms using " << nBlock << " blocks of length " << nElement << " each.\n# Time between frames: " << deltaT << " [ps]\n# Number of atom types: " << nAtomType << "\n# Number of molecule type: " << nMoleculeType << "\n# Atomic masses" << "\n";

    PostProcess pp(nBlock,nElement,deltaT,nAtomType,nMoleculeType,nAtom);

    getlines(fin,sline);
    for (int i=0;i<nAtomType;i++) {
        getlines(fin,sline);
        sline >> name >> mass; // atomic masses for each atom type
        pp.setAtomType(i,name,mass);
        pp.atomMap[name]=i;
        cout << "# " << name << " " << mass << "\n";
    }

    cout << "# Molecule composition";
    getlines(fin,sline);
    for (int i=0;i<nMoleculeType;i++) {
        getlines(fin,sline);
        sline >> name;
        pp.moleculeMap[name]=i;
        cout << "\n# " << name << ": ";
        atoms.clear();
        while (sline>>type) {
            atoms.push_back(type); // atom types for each molecule
            cout << type << " ";
        }
        pp.setMoleculeType(i,name,atoms);
    }

    cout << "\n# Movie files:";
    getlines(fin,sline);
    while (fin >> fname) {
        cout << " " << fname;
        pp.proc(fname.c_str());
    }
    fin.close();
    cout << "\n";
    pp.output(argv[2]);
    return 0;
}
