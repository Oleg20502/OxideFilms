#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <sstream>
#include <filesystem>
#include <cmath>


using std::array, std::vector, std::size_t, std::string;
using std::pow, std::exp, std::stoi;


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


void gen_init_c(vector<double>& C_MV_0, vector<double>& C_OV_0,
                vector<double>& x, const int s1, const int s2, Params P)
{
    int Nx = x.size();
    for (int i = 0; i < Nx; ++i) {
        if (i < s1) {
            C_MV_0[i] = 150 / x[s1] * x[i] / P.A_C;
        }
        else {
            C_MV_0[i] = 150 / P.A_C;
        }
    }

    for (int i = 0; i < Nx; ++i) {
        if (i <= Nx - s2) {
            C_MV_0[i] = (210 - (210 - 185) / x[Nx - s2] * x[i]) / P.A_C;
        }
        else {
            C_MV_0[i] = 185 * (P.L - x[i]) / (P.L - x[Nx-s2]) / P.A_C;
        }
    }
}


int solve(Params P)
{
    const int Nx = P.Nx;
    const int Nt = P.Nt;
    const int Ndata = P.Ndata;
    const int n_save = P.n_save;
    const double h = P.h;
    const double dt = P.dt;


    const double D_MV = P.D_MV;
    const double D_OV = P.D_OV;

    const double k1_0 = P.k1_0;
    const double k2_0 = P.k2_0;
    const double k3_0 = P.k3_0;
    const double k4_0 = P.k4_0;
    const double k5_0 = P.k5_0;

    const double phi_ext = P.phi_ext;
    const double L = P.L;
    const double T = P.T;

    const double an = P.an;
    const double e_f = P.e_f;
    const double e_dl = P.e_dl;
    const double e_cdl = P.e_cdl;
    const double d_dl = P.d_dl;
    const double d_cdl = P.d_cdl;

    const double Temp = P.Temp;
    const double Ff = P.Ff;
    const double e_0 = P.e_0;
    const double R = P.R;

    const double A_k = P.A_k;
    const double A_D = P.A_D;
    const double A_L = P.A_L;
    const double A_C = P.A_C;
    const double A_t = P.A_t;
    const double A_phi = P.A_phi;

    const double A1 = -pow(Ff, 2) /e_f/e_0/R/Temp*(A_k*pow(A_L, 3)/A_D)*2;
    const double A2 = -e_dl/e_f/d_dl;
    const double A3 = e_f*(d_dl/e_dl + d_cdl/e_cdl);
    const double A4 = dt*0.5/h;
    const double A_MV = D_MV * dt/ pow(h, 2);
    const double A_OV = D_OV * dt/ pow(h, 2);


    vector<double> C_MV(Nx, 0.0);
    vector<double> C_OV(Nx, 0.0);
    vector<double> C_MV_2(Nx, 0.0);
    vector<double> C_OV_2(Nx, 0.0);
    vector<double> C_V(Nx, 0.0);

    vector<double> phi(Nx, 0.0);
    vector<double> E(Nx-1, 0.0);

    vector<double> A(Nx, 0.0);
    vector<double> B(Nx, 0.0);
    vector<double> C(Nx, 0.0);
    vector<double> d(Nx, 0.0);
    vector<double> Cm(Nx-1, 0.0);
    vector<double> dm(Nx-1, 0.0);


    vector<double> x(Nx, 0.0);
    for (int i = 0; i < Nx; ++i) {
        x[i] = i * h;
    }
    
    vector<double> t(Nt, 0.0);
    for (int i = 0; i < Nt; ++i) {
        t[i] = i * dt;
    }

    vector<double> t_data(P.Ndata, 0.0);
    vector<vector<double>> Data_C_MV(P.Ndata, vector<double>(P.Nx, 0.0));
    vector<vector<double>> Data_C_OV(P.Ndata, vector<double>(P.Nx, 0.0));
    vector<vector<double>> Data_phi(P.Ndata, vector<double>(P.Nx, 0.0));
    vector<vector<double>> Data_E(P.Ndata, vector<double>(P.Nx - 1, 0.0));
    vector<double> Data_k2(P.Ndata, 0.0);

    vector<double> C_MV_0(P.Nx, 0.0);
    vector<double> C_OV_0(P.Nx, 0.0);
    
    gen_init_c(C_MV_0, C_OV_0, x, round(0.1 * P.L / P.h), round(0.1 * P.L / P.h), P);


    for (int i = 0; i < Nx; ++i) {
        C_V[i] = C_OV[i] - C_MV[i];
    }
    B[0] = -1/h + A2;
    C[0] = 1/h;
    d[0] = A2*phi_ext;
    B[Nx-1] = A3/h + 1;
    A[Nx-1] = -A3/h;
    d[Nx-1] = 0;
    for (int i = 1; i < Nx-1; ++i) {
        B[i] = -2/pow(h, 2);
        C[i] = 1/pow(h, 2);
        A[i] = 1/pow(h, 2);
        d[i] = A1*C_V[i];
    }
    Cm[0] = C[0]/B[0];
    dm[0] = d[0]/B[0];

    for (int i = 1; i < Nx-1; ++i) {
        Cm[i] = C[i]/(B[i] - A[i]*Cm[i-1]);
        dm[i] = (d[i] - A[i]*dm[i-1])/(B[i] - A[i]*Cm[i-1]);
    }
    phi[Nx-1] = (d[Nx-1] - A[Nx-1]*dm[Nx-2])/(B[Nx-1] - A[Nx-1]*Cm[Nx-2]);
    for (int i = Nx-2; i >= 0; --i) {
        phi[i] = dm[i] - Cm[i] * phi[i+1];
    }
    for (int i = 0; i < Nx-1; ++i) {
        E[i] = -(phi[i+1] - phi[i])/h;
    }
    double fx_0 = phi[1] - phi[0];
    double fx_L = phi[Nx-1] - phi[Nx-2];
    double phi_mf = phi_ext - phi[0];
    double phi_fs = phi[Nx-1];
    double k1 = k1_0 * exp(an*phi_mf);
    double k2 = k2_0 * exp(an*phi_mf);
    double k3 = k3_0 * exp(an*phi_fs);
    double k4 = k4_0 * exp(an*phi_fs);


    for (int i = 0; i < Nx; ++i) {
        Data_C_MV[0][i] = C_MV[i];
        Data_C_OV[0][i] = C_OV[i];
        Data_phi[0][i] = phi[i];
        Data_E[0][i] = E[i];
    }
    Data_k2[0] = k2;

    double fx;
    int I_data = 1;
    for (int i = 1; i < Nt; ++i) {
        for (int j = 1; j < Nx-1; ++j) {
            fx = phi[j+1]-phi[j-1];
            C_MV_2[j] = C_MV[j] + A_MV*(C_MV[j+1] - 2*C_MV[j] + C_MV[j-1]) - 2*D_MV*A1*dt*C_V[j]*C_MV[j] - 0.5*A_MV*fx * (C_MV[j+1]-C_MV[j-1]);
            C_OV_2[j] = C_OV[j] + A_OV*(C_OV[j+1] - 2*C_OV[j] + C_OV[j-1]) + 2*D_OV*A1*dt*C_V[j]*C_OV[j] + 0.5*A_OV*fx * (C_OV[j+1]-C_OV[j-1]);
        }

        // граничные для MV
        C_MV_2[0] = C_MV_2[1] / (1 + k1*h/D_MV + 2*fx_0);
        C_MV_2[Nx-1] = (C_MV_2[Nx-2] + k3*h/D_MV)/(1 - 2*fx_L);
        
        // граничные для OV
        C_OV_2[0] = (C_OV_2[1] + k2*h/D_OV)/(1 - 2*fx_0);
        C_OV_2[Nx-1] = C_OV_2[Nx-2] / (1 + k4*h/D_OV + 2*fx_L);

        for (int m = 0; m < Nx; ++m) {
            C_MV[m] = C_MV_2[m];
            C_OV[m] = C_OV_2[m];
            C_V[m] = C_OV[m] - C_MV[m];
        }

        for (int m = 1; m < Nx-1; ++m) {
            d[m] = A1*C_V[m];
        }
        for (int p = 1; p < Nx-1; ++p) {
            Cm[p] = C[p]/(B[p] - A[p]*Cm[p-1]);
            dm[p] = (d[p] - A[p]*dm[p-1])/(B[p] - A[p]*Cm[p-1]);
        }
        phi[Nx-1] = (d[Nx-1] - A[Nx-1]*dm[Nx-2])/(B[Nx-1] - A[Nx-1]*Cm[Nx-2]);
        for (int p = Nx-2; p >= 0; --p) {
            phi[p] = dm[p] - Cm[p] * phi[p+1];
        }
        for (int i = 0; i < Nx-1; ++i) {
            E[i] = -(phi[i+1] - phi[i])/h;
        }
        fx_0 = phi[1] - phi[0];
        fx_L = phi[Nx-1] - phi[Nx-2];
        phi_mf = phi_ext - phi[0];
        phi_fs = phi[Nx-1];
        k1 = k1_0 * exp(an*phi_mf);
        k2 = k2_0 * exp(an*phi_mf);
        k3 = k3_0 * exp(an*phi_fs);
        k4 = k4_0 * exp(an*phi_fs);

        if (i % n_save == 0) {
            for (int p = 0; p < Nx; ++p) {
                Data_C_MV[I_data][p] = C_MV[p];
                Data_C_OV[I_data][p] = C_OV[p];
                Data_phi[I_data][p] = phi[p];
                Data_E[I_data][p] = E[p];
                t_data[I_data] = t[i];
            }
            Data_k2[I_data] = k2;
            std::cout << i << "  ";
            ++I_data;
        }
    }
    for (int p = 0; p < Nx; ++p) {
        Data_C_MV[I_data][p] = C_MV[p];
        Data_C_OV[I_data][p] = C_OV[p];
        Data_phi[I_data][p] = phi[p];
        Data_E[I_data][p] = E[p];
        t_data[I_data] = t[Nt-1];
    }
    Data_k2[I_data] = k2;
    std::cout << Nt-1 << '\n' << "Finished\n";

    return 0;
}


int save(vector<double> & t_data, vector<vector<double>> & Data_C_MV, 
        vector<vector<double>> & Data_C_OV, vector<vector<double>> & Data_phi, 
        vector<vector<double>> & Data_E, vector<double> & Data_k2, 
        Params P, string res_path, int dig = 10)
{
    std::ofstream input_t_data;
    input_t_data.open(res_path + "/t_data.txt", std::ios::out);
    input_t_data << std::fixed << std::showpoint << std::setprecision(dig);
    for (int i = 0; i < P.Ndata; ++i) {
        input_t_data << t_data[i] << '\n';
    }
    input_t_data.close();

    std::ofstream input_Data_C_MV;
    input_Data_C_MV.open(res_path + "/Data_C_MV.txt", std::ios::out);
    input_Data_C_MV << std::fixed << std::showpoint << std::setprecision(dig);
    for (int i = 0; i < P.Ndata; ++i) {
        for (int j = 0; j < P.Nx-1; ++j) {
            input_Data_C_MV << Data_C_MV[i][j] << ' ';
        }
        input_Data_C_MV << Data_C_MV[i][P.Nx-1] << '\n';
    }
    input_Data_C_MV.close();

    std::ofstream input_Data_C_OV;
    input_Data_C_OV.open(res_path + "/Data_C_OV.txt", std::ios::out);
    input_Data_C_OV << std::fixed << std::showpoint << std::setprecision(dig);
    for (int i = 0; i < P.Ndata; ++i) {
        for (int j = 0; j < P.Nx - 1; ++j) {
            input_Data_C_OV << Data_C_OV[i][j] << ' ';
        }
        input_Data_C_OV << Data_C_OV[i][P.Nx-1] << '\n';
    }
    input_Data_C_OV.close();

    std::ofstream input_Data_phi;
    input_Data_phi.open(res_path + "/Data_phi.txt", std::ios::out);
    input_Data_phi << std::fixed << std::showpoint << std::setprecision(dig);
    for (int i = 0; i < P.Ndata; ++i) {
        for (int j = 0; j < P.Nx - 1; ++j) {
            input_Data_phi << Data_phi[i][j] << ' ';
        }
        input_Data_phi << Data_phi[i][P.Nx-1] << '\n';
    }
    input_Data_phi.close();

    std::ofstream input_Data_E;
    input_Data_E.open(res_path + "/Data_E.txt", std::ios::out);
    input_Data_E << std::fixed << std::showpoint << std::setprecision(dig);
    for (int i = 0; i < P.Ndata; ++i) {
        for (int j = 0; j < P.Nx - 2; ++j) {
            input_Data_E << Data_E[i][j] << ' ';
        }
        input_Data_E << Data_E[i][P.Nx - 2] << '\n';
    }
    input_Data_E.close();

    std::ofstream input_Data_k2;
    input_Data_k2.open(res_path + "/Data_k2.txt", std::ios::out);
    input_Data_k2 << std::fixed << std::showpoint << std::setprecision(dig);
    for (int i = 0; i < P.Ndata; ++i) {
        input_Data_k2 << Data_k2[i] << '\n';
    }
    input_Data_k2.close();

    return 0;
}


int main(int argc, char* argv[])
{
    string res_path;
    int dig;
    int check = 0;

    Params P;
    string par_file = "C:/Users/OlegKashurin/home/Work/Films/OxideFilms/R-PDM/Release/params.txt";
    std::ifstream par_input;
    par_input.open(par_file, std::ios::in);
    if (!par_input) {
        std::cerr << "Unable to open file " << par_file << '\n';
        exit(1);
    }
    for (string line; getline(par_input, line); ) {
        std::stringstream input_line(line);
        string word1, word2;
        input_line >> word1 >> word2;
        if (word1 == "Nx") { P.Nx = stoi(word2); ++check; }
        else if (word1 == "n_save") { P.n_save = stoi(word2); ++check; }
        else if (word1 == "res_path") { res_path = word2; ++check; }
        else if (word1 == "phi_ext") { P.phi_ext = stod(word2); ++check; }
        else if (word1 == "L") { P.L = stod(word2); ++check; }
        else if (word1 == "T") { P.T = stod(word2); ++check; }
        else if (word1 == "D_MV") { P.D_MV = stod(word2); ++check; }
        else if (word1 == "D_OV") { P.D_OV = stod(word2); ++check; }
        else if (word1 == "k1_0") { P.k1_0 = stod(word2); ++check; }
        else if (word1 == "k2_0") { P.k2_0 = stod(word2); ++check; }
        else if (word1 == "k3_0") { P.k3_0 = stod(word2); ++check; }
        else if (word1 == "k4_0") { P.k4_0 = stod(word2); ++check; }
        else if (word1 == "k5_0") { P.k5_0 = stod(word2); ++check; }
        else if (word1 == "an") { P.an = stod(word2); ++check; }
        else if (word1 == "e_f") { P.e_f = stod(word2); ++check; }
        else if (word1 == "e_dl") { P.e_dl = stod(word2); ++check; }
        else if (word1 == "e_cdl") { P.e_cdl = stod(word2); ++check; }
        else if (word1 == "d_dl") { P.d_dl = stod(word2); ++check; }
        else if (word1 == "d_cdl") { P.d_cdl = stod(word2); ++check; }
        else if (word1 == "Temp") { P.Temp = stod(word2); ++check; }
        else if (word1 == "Ff") { P.Ff = stod(word2); ++check; }
        else if (word1 == "e_0") { P.e_0 = stod(word2); ++check; }
        else if (word1 == "R") { P.R = stod(word2); ++check; }
        else if (word1 == "A_k") { P.A_k = stod(word2); ++check; }
        else if (word1 == "A_D") { P.A_D = stod(word2); ++check; }
        else if (word1 == "A_L") { P.A_L = stod(word2); ++check; }
        else if (word1 == "dig") { dig = stoi(word2); ++check; }
        else {
            std::cerr << "Error: Argument " << word1 <<  " not recognized.\n";
            exit(-1);
        }       
    }
    par_input.close();

    if (check != 27) {
        std::cerr << "Error: Not all parameters were set.\n";
        exit(check);
    }
    
    P.A_C = P.A_k * P.A_L / P.A_D;
    P.A_t = pow(P.A_L, 2) / P.A_D;
    P.A_phi = P.R * P.Temp / P.Ff;
    P.phi_ext /= P.A_phi;

    P.h = P.L / (P.Nx-1);
    /*
    P.dt = 0.5 * pow(P.h, 2) / (P.D_MV + P.D_OV);
    P.Nt = t.size();*/

    P.Nt = ceil(2 * (P.D_MV + P.D_OV) * P.T / pow(P.h, 2)) + 1;
    P.dt = P.T / (P.Nt-1);
    P.Ndata = (P.Nt - 1) / P.n_save + 2;

    std::cout << "Nt " << P.Nt << '\n';
    std::cout << "Ndata " << P.Ndata << '\n';

    solve(P);

    
    return 0;
}