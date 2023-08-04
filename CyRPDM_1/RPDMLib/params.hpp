#ifndef PARAMS_H
#define PARAMS_H

class Params
{
    public:
        int Nx;
        int Nt;
        int Ndata;
        int n_save;
        double h;
        double dt;

        double phi_ext;
        double L;
        double T;

        double D_MV;
        double D_OV;

        double k1_0;
        double k2_0;
        double k3_0;
        double k4_0;
        double k5_0;

        double an;
        double e_f;
        double e_dl;
        double e_cdl;
        double d_dl;
        double d_cdl;

        double Temp;
        double Ff;
        double e_0;
        double R;

        double A_k;
        double A_D;
        double A_L;
        double A_C;
        double A_t;
        double A_phi;
};

#endif