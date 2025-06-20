        ifstream fin(fname);
        istringstream sline;
        int nStep,atomID,molecID,atomType;
        double tmp,r[3];

        while (getlines(fin,sline)) { //ITEM: TIMESTEP
            getlines(fin,sline);
            sline >> nStep;
            getlines(fin,sline); //ITEM: NUMBER OF ATOMS
            getlines(fin,sline);
            sline >> nAtom;
            if (!initialized) {
                mol=new Molecule[nAtom];
                data=new double***[nBlk];
                for (int i=0;i<nBlk;i++) {
                    data[i]=new double**[nElem];
                    for (int j=0;j<nElem;j++) {
                        data[i][j]=new double*[nAtom];
                        for (int k=0;k<nAtom;k++)
                            data[i][j][k]=new double[3];
                    }
                }
            }
            getlines(fin,sline); //ITEM: BOX BOUNDS pp pp pp
            getlines(fin,sline);
            sline >> tmp >> box.boxl[0]; //x dimension, low bound, high bound
            getlines(fin,sline);
            sline >> tmp >> box.boxl[1]; //y
            getlines(fin,sline);
            sline >> tmp >> box.boxl[2]; //z
            getlines(fin,sline); //ITEM: ATOMS id type xu yu zu
            for (int i=0;i<nAtom;i++) {
                getlines(fin,sline);
                sline >> atomID >> atomType >> r[0] >> r[1] >> r[2];
                //sline >> atomID >> molecID >> atomType >> r[0] >> r[1] >> r[2];
                //cerr << atomID << " " << molecID << " " << atomType <<  " " << r[0] <<  " " << r[1] <<  " " << r[2] << endl;
                if (!initialized) atomMap[atomID]=i;
                mol[atomMap[atomID]].type=atomType;
                for (int j=0;j<3;j++) mol[atomMap[atomID]].com[j]=r[j];
            }
            if (!initialized) initialized=true;
            accumulateMsd(++iFrame);
        }
        fin.close();
