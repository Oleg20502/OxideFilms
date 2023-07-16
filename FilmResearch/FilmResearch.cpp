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
            C_MV_0[i] = 185 * (P.L - x[i]) / (P.L - x[Nx - s2]) / P.A_C;
        }
    }
}


int solve(vector<double>& x, vector<double>& t,
    vector<double>& C_MV_0, vector<double>& C_OV_0,
    vector<double>& t_data, vector<vector<double>>& Data_C_MV,
    vector<vector<double>>& Data_C_OV, vector<vector<double>>& Data_phi,
    vector<vector<double>>& Data_E, vector<double>& Data_k2,
    Params P)
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

    const double A1 = -pow(Ff, 2) / e_f / e_0 / R / Temp * (A_k * pow(A_L, 3) / A_D) * 2;
    const double A2 = -e_dl / e_f / d_dl;
    const double A3 = e_f * (d_dl / e_dl + d_cdl / e_cdl);
    const double A4 = dt * 0.5 / h;
    const double A_MV = D_MV * dt / pow(h, 2);
    const double A_OV = D_OV * dt / pow(h, 2);


    vector<double> C_MV = C_MV_0;
    vector<double> C_OV = C_OV_0;
    vector<double> C_MV_2(Nx, 0.0);
    vector<double> C_OV_2(Nx, 0.0);
    vector<double> C_V(Nx, 0.0);

    vector<double> phi(Nx, 0.0);
    vector<double> E(Nx - 1, 0.0);

    vector<double> A(Nx, 0.0);
    vector<double> B(Nx, 0.0);
    vector<double> C(Nx, 0.0);
    vector<double> d(Nx, 0.0);
    vector<double> Cm(Nx - 1, 0.0);
    vector<double> dm(Nx - 1, 0.0);

    for (int i = 0; i < Nx; ++i) {
        C_V[i] = C_OV[i] - C_MV[i];
    }
    B[0] = -1 / h + A2;
    C[0] = 1 / h;
    d[0] = A2 * phi_ext;
    B[Nx - 1] = A3 / h + 1;
    A[Nx - 1] = -A3 / h;
    d[Nx - 1] = 0;
    for (int i = 1; i < Nx - 1; ++i) {
        B[i] = -2 / pow(h, 2);
        C[i] = 1 / pow(h, 2);
        A[i] = 1 / pow(h, 2);
        d[i] = A1 * C_V[i];
    }
    Cm[0] = C[0] / B[0];
    dm[0] = d[0] / B[0];

    for (int i = 1; i < Nx - 1; ++i) {
        Cm[i] = C[i] / (B[i] - A[i] * Cm[i - 1]);
        dm[i] = (d[i] - A[i] * dm[i - 1]) / (B[i] - A[i] * Cm[i - 1]);
    }
    phi[Nx - 1] = (d[Nx - 1] - A[Nx - 1] * dm[Nx - 2]) / (B[Nx - 1] - A[Nx - 1] * Cm[Nx - 2]);
    for (int i = Nx - 2; i >= 0; --i) {
        phi[i] = dm[i] - Cm[i] * phi[i + 1];
    }
    for (int i = 0; i < Nx - 1; ++i) {
        E[i] = -(phi[i + 1] - phi[i]) / h;
    }
    double fx_0 = phi[1] - phi[0];
    double fx_L = phi[Nx - 1] - phi[Nx - 2];
    double phi_mf = phi_ext - phi[0];
    double phi_fs = phi[Nx - 1];
    double k1 = k1_0 * exp(an * phi_mf);
    double k2 = k2_0 * exp(an * phi_mf);
    double k3 = k3_0 * exp(an * phi_fs);
    double k4 = k4_0 * exp(an * phi_fs);


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
        for (int j = 1; j < Nx - 1; ++j) {
            fx = phi[j + 1] - phi[j - 1];
            C_MV_2[j] = C_MV[j] + A_MV * (C_MV[j + 1] - 2 * C_MV[j] + C_MV[j - 1]) - 2 * D_MV * A1 * dt * C_V[j] * C_MV[j] - 0.5 * A_MV * fx * (C_MV[j + 1] - C_MV[j - 1]);
            C_OV_2[j] = C_OV[j] + A_OV * (C_OV[j + 1] - 2 * C_OV[j] + C_OV[j - 1]) + 2 * D_OV * A1 * dt * C_V[j] * C_OV[j] + 0.5 * A_OV * fx * (C_OV[j + 1] - C_OV[j - 1]);
        }

        // граничные для MV
        C_MV_2[0] = C_MV_2[1] / (1 + k1 * h / D_MV + 2 * fx_0);
        C_MV_2[Nx - 1] = (C_MV_2[Nx - 2] + k3 * h / D_MV) / (1 - 2 * fx_L);

        // граничные для OV
        C_OV_2[0] = (C_OV_2[1] + k2 * h / D_OV) / (1 - 2 * fx_0);
        C_OV_2[Nx - 1] = C_OV_2[Nx - 2] / (1 + k4 * h / D_OV + 2 * fx_L);

        for (int m = 0; m < Nx; ++m) {
            C_MV[m] = C_MV_2[m];
            C_OV[m] = C_OV_2[m];
            C_V[m] = C_OV[m] - C_MV[m];
        }

        for (int m = 1; m < Nx - 1; ++m) {
            d[m] = A1 * C_V[m];
        }
        for (int p = 1; p < Nx - 1; ++p) {
            Cm[p] = C[p] / (B[p] - A[p] * Cm[p - 1]);
            dm[p] = (d[p] - A[p] * dm[p - 1]) / (B[p] - A[p] * Cm[p - 1]);
        }
        phi[Nx - 1] = (d[Nx - 1] - A[Nx - 1] * dm[Nx - 2]) / (B[Nx - 1] - A[Nx - 1] * Cm[Nx - 2]);
        for (int p = Nx - 2; p >= 0; --p) {
            phi[p] = dm[p] - Cm[p] * phi[p + 1];
        }
        for (int i = 0; i < Nx - 1; ++i) {
            E[i] = -(phi[i + 1] - phi[i]) / h;
        }
        fx_0 = phi[1] - phi[0];
        fx_L = phi[Nx - 1] - phi[Nx - 2];
        phi_mf = phi_ext - phi[0];
        phi_fs = phi[Nx - 1];
        k1 = k1_0 * exp(an * phi_mf);
        k2 = k2_0 * exp(an * phi_mf);
        k3 = k3_0 * exp(an * phi_fs);
        k4 = k4_0 * exp(an * phi_fs);

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
        t_data[I_data] = t[Nt - 1];
    }
    Data_k2[I_data] = k2;
    std::cout << Nt - 1 << '\n' << "Finished\n";

    return 0;
}


int calc_L_phi(double L_b, double L_e, double dL, Params P)
{

    return 0;
}


int main()
{
    int Nx, Ndata;
    double h;
    double L_b, L_e, dL;
    int n_dL;
    string res_path;
    int dig;
    int check = 0;

    Params Par;
    string par_file = "C:/Users/OlegKashurin/home/Work/Films/OxideFilms/FilmResearch/Release/params.txt";
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
        if (word1 == "L_b") { L_b = stod(word2); ++check; }
        else if (word1 == "L_e") { L_e = stod(word2); ++check; }
        else if (word1 == "dL") { dL = stod(word2); ++check; }  // может быть n_dL
        else if (word1 == "Nx") { Par.Nx = stoi(word2); ++check; }
        else if (word1 == "n_save") { Par.n_save = stoi(word2); ++check; }
        else if (word1 == "res_path") { res_path = word2; ++check; }
        else if (word1 == "phi_ext") { Par.phi_ext = stod(word2); ++check; }
        else if (word1 == "L") { Par.L = stod(word2); ++check; }
        else if (word1 == "T") { Par.T = stod(word2); ++check; }
        else if (word1 == "D_MV") { Par.D_MV = stod(word2); ++check; }
        else if (word1 == "D_OV") { Par.D_OV = stod(word2); ++check; }
        else if (word1 == "k1_0") { Par.k1_0 = stod(word2); ++check; }
        else if (word1 == "k2_0") { Par.k2_0 = stod(word2); ++check; }
        else if (word1 == "k3_0") { Par.k3_0 = stod(word2); ++check; }
        else if (word1 == "k4_0") { Par.k4_0 = stod(word2); ++check; }
        else if (word1 == "k5_0") { Par.k5_0 = stod(word2); ++check; }
        else if (word1 == "an") { Par.an = stod(word2); ++check; }
        else if (word1 == "e_f") { Par.e_f = stod(word2); ++check; }
        else if (word1 == "e_dl") { Par.e_dl = stod(word2); ++check; }
        else if (word1 == "e_cdl") { Par.e_cdl = stod(word2); ++check; }
        else if (word1 == "d_dl") { Par.d_dl = stod(word2); ++check; }
        else if (word1 == "d_cdl") { Par.d_cdl = stod(word2); ++check; }
        else if (word1 == "Temp") { Par.Temp = stod(word2); ++check; }
        else if (word1 == "Ff") { Par.Ff = stod(word2); ++check; }
        else if (word1 == "e_0") { Par.e_0 = stod(word2); ++check; }
        else if (word1 == "R") { Par.R = stod(word2); ++check; }
        else if (word1 == "A_k") { Par.A_k = stod(word2); ++check; }
        else if (word1 == "A_D") { Par.A_D = stod(word2); ++check; }
        else if (word1 == "A_L") { Par.A_L = stod(word2); ++check; }
        else if (word1 == "dig") { dig = stoi(word2); ++check; }
        else {
            std::cerr << "Error: Argument " << word1 << " not recognized.\n";
            exit(-1);
        }
    }
    par_input.close();

    if (check != 30) {
        std::cerr << "Error: Not all parameters were set.\n";
        exit(check);
    }

    Par.A_C = Par.A_k * Par.A_L / Par.A_D;
    Par.A_t = pow(Par.A_L, 2) / Par.A_D;
    Par.A_phi = Par.R * Par.Temp / Par.Ff;


    Par.phi_ext /= Par.A_phi; // v, что-то еще


    Nx = Par.Nx;
    Par.h = Par.L / (Nx - 1);
    vector<double> x(Nx, 0.0);
    for (int i = 0; i < Nx; ++i) {
        x[i] = i * Par.h;
    }
    /*
    Par.dt = 0.5 * pow(Par.h, 2) / (Par.D_MV + Par.D_OV);
    vector<double> t;
    int tmp = 0;
    while (Par.dt * tmp < Par.T + Par.dt) {
        t.push_back(Par.dt * tmp);
        ++tmp;
    }
    Par.Nt = t.size();*/

    Par.Nt = ceil(2 * (Par.D_MV + Par.D_OV) * Par.T / pow(Par.h, 2)) + 1;
    Par.dt = Par.T / (Par.Nt - 1);
    vector<double> t(Par.Nt, 0.0);
    for (int i = 0; i < Par.Nt; ++i) {
        t[i] = i * Par.dt;
    }

    Par.Ndata = (Par.Nt - 1) / Par.n_save + 2;
    Ndata = Par.Ndata;
    vector<double> t_data(Ndata, 0.0);
    vector<vector<double>> Data_C_MV(Ndata, vector<double>(Nx, 0.0));
    vector<vector<double>> Data_C_OV(Ndata, vector<double>(Nx, 0.0));
    vector<vector<double>> Data_phi(Ndata, vector<double>(Nx, 0.0));
    vector<vector<double>> Data_E(Ndata, vector<double>(Nx - 1, 0.0));
    vector<double> Data_k2(Ndata, 0.0);

    vector<double> C_MV_0(Nx, 0.0);
    vector<double> C_OV_0(Nx, 0.0);
    gen_init_c(C_MV_0, C_OV_0, x, round(0.1 * Par.L / Par.h), round(0.1 * Par.L / Par.h), Par);

    std::cout << "Nt " << Par.Nt << '\n';
    std::cout << "Ndata " << Ndata << '\n';
    std::cout << "phi_ext " << Par.phi_ext << '\n';

    calc_L_phi(L_b, L_e, dL, Par);




    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
