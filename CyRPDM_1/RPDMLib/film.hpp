#ifndef FILM_H
#define FILM_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>


using std::vector;
using std::string;
using std::pow;
using std::exp;


class Film
{
    private:
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

        double an1;
        double an2;
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

        double A0;
        double A1;
        double A2;
        double A3;
        double A4;
        double A_MV;
        double A_OV;

        vector<double> C_MV;
        vector<double> C_OV;
        vector<double> C_MV_2;
        vector<double> C_OV_2;
        vector<double> C_MV_0;
        vector<double> C_OV_0;
        vector<double> C_V;

        vector<double> phi;
        vector<double> E;

        vector<double> A;
        vector<double> B;
        vector<double> C;
        vector<double> d;
        vector<double> Cm;
        vector<double> dm;

        vector<double> x;
        vector<double> t;

        vector<double> t_data;
        vector<vector<double>> Data_C_MV;
        vector<vector<double>> Data_C_OV;
        vector<vector<double>> Data_phi;
        vector<vector<double>> Data_E;
        vector<double> Data_k2;


        int check_init;

    public:
        Film()
        {
            check_init = 0;
        }

        int init()
        {
            if (check_init < 27) {
                std::cerr << "Error: Not all parameters were set.\n";
                return check_init;
            }
            A_C = A_k * A_L / A_D;
            A_t = pow(A_L, 2) / A_D;

            h = L / (Nx-1);
            Nt = ceil(2 * (D_MV + D_OV) * T / pow(h, 2)) + 1;
            dt = T / (Nt-1);
            Ndata = (Nt - 1) / n_save + 2;

            A0 = A_phi * Ff / (R * Temp);
            A1 = -Ff/e_f/e_0*(A_k*pow(A_L, 3)/A_D)*2 / A_phi;
            A2 = -e_dl/e_f/d_dl;
            A3 = e_f*(d_dl/e_dl + d_cdl/e_cdl);
            A4 = dt*0.5/h;
            A_MV = D_MV * dt/ pow(h, 2);
            A_OV = D_OV * dt/ pow(h, 2);

            C_MV.resize(Nx);
            C_OV.resize(Nx);
            C_MV_2.resize(Nx);
            C_OV_2.resize(Nx);
            C_V.resize(Nx);

            phi.resize(Nx);
            E.resize(Nx-1);

            A.resize(Nx);
            B.resize(Nx);
            C.resize(Nx);
            d.resize(Nx);
            Cm.resize(Nx-1);
            dm.resize(Nx-1);

            x.resize(Nx);
            for (int i = 0; i < Nx; ++i) {
                x[i] = i * h;
            }
            
            t.resize(Nt);
            for (int i = 0; i < Nt; ++i) {
                t[i] = i * dt;
            }

            t_data.resize(Ndata);
            Data_k2.resize(Ndata);
            Data_C_MV.resize(Ndata);
            Data_C_OV.resize(Ndata);
            Data_phi.resize(Ndata);
            Data_E.resize(Ndata);
            for (int i = 0; i < Ndata; ++i) {
                Data_C_MV[i].resize(Nx);
                Data_C_OV[i].resize(Nx);
                Data_phi[i].resize(Nx);
                Data_E[i].resize(Nx-1);
            }

            /* C_MV_0 = vector<double> (Nx, 0.0);
            C_OV_0 = vector<double> (Nx, 0.0); */
            
            init_C(C_MV, C_OV, x, round(0.1 * L / h), round(0.1 * L / h), A_C);
            
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

            std::cout << "Nx " << Nx << '\n';
            std::cout << "Nt " << Nt << '\n';
            std::cout << "Ndata " << Ndata << '\n';

            return check_init;
        }

        int solve()
        {         
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
            double k1 = k1_0 * exp(an1*A0*phi_mf);
            double k2 = k2_0 * exp(an1*A0*phi_mf);
            double k3 = k3_0 * exp(an2*A0*phi_fs);
            double k4 = k4_0 * exp(an2*A0*phi_fs);

            for (int i = 0; i < Nx; ++i) {
                Data_C_MV[0][i] = C_MV[i];
                Data_C_OV[0][i] = C_OV[i];
                Data_phi[0][i] = phi[i];
            }
            for (int i = 0; i < Nx-1; ++i) {
                Data_E[0][i] = E[i];
            }
            Data_k2[0] = k2;

            double fx;
            int I_data = 1;
            for (int i = 1; i < Nt; ++i) {
                for (int j = 1; j < Nx-1; ++j) {
                    fx = (phi[j+1]-phi[j-1]);
                    C_MV_2[j] = C_MV[j] + A_MV*(C_MV[j+1]-2*C_MV[j]+C_MV[j-1]) - 2*D_MV*dt*d[j]*C_MV[j] - 0.5*A_MV*fx*A0*(C_MV[j+1]-C_MV[j-1]);
                    C_OV_2[j] = C_OV[j] + A_OV*(C_OV[j+1]-2*C_OV[j]+C_OV[j-1]) + 2*D_OV*dt*d[j]*C_OV[j] + 0.5*A_OV*fx*A0*(C_OV[j+1]-C_OV[j-1]);
                }

                // граничные для MV
                C_MV_2[0] = C_MV_2[1] / (1 + k1*h/D_MV + 2*fx_0*A0);
                C_MV_2[Nx-1] = (C_MV_2[Nx-2] + k3*h/D_MV)/(1 - 2*fx_L*A0);
                
                // граничные для OV
                C_OV_2[0] = (C_OV_2[1] + k2*h/D_OV)/(1 - 2*fx_0*A0);
                C_OV_2[Nx-1] = C_OV_2[Nx-2] / (1 + k4*h/D_OV + 2*fx_L*A0);

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
                k1 = k1_0 * exp(an1*A0*phi_mf);
                k2 = k2_0 * exp(an1*A0*phi_mf);
                k3 = k3_0 * exp(an2*A0*phi_fs);
                k4 = k4_0 * exp(an2*A0*phi_fs);

                if (i % n_save == 0) {
                    for (int p = 0; p < Nx; ++p) {
                        Data_C_MV[I_data][p] = C_MV[p];
                        Data_C_OV[I_data][p] = C_OV[p];
                        Data_phi[I_data][p] = phi[p];
                        t_data[I_data] = t[i];
                    }
                    for (int p = 0; p < Nx-1; ++p) {
                        Data_E[I_data][p] = E[p];
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
                t_data[I_data] = t[Nt-1];
            }
            for (int p = 0; p < Nx-1; ++p) {
                Data_E[I_data][p] = E[p];
            }
            Data_k2[I_data] = k2;
            std::cout << Nt-1 << '\n' << "Finished\n";

            return 0;
        }

        static void init_C(vector<double>& C_MV_0, vector<double>& C_OV_0,
                           vector<double>& x, int s1, int s2, double A_C)
        {
            int Nx = x.size();
            for (int i = 0; i < Nx; ++i) {
                if (i < s1) {
                    C_MV_0[i] = 150 / x[s1] * x[i] / A_C;
                }
                else {
                    C_MV_0[i] = 150 / A_C;
                }
            }

            for (int i = 0; i < Nx; ++i) {
                if (i <= Nx - s2) {
                    C_OV_0[i] = (210 - (210 - 185) / x[Nx - s2] * x[i]) / A_C;
                }
                else {
                    C_OV_0[i] = 185 * (x[Nx-1] - x[i]) / (x[Nx-1] - x[Nx-s2]) / A_C;
                }
            }
        }

        void save_x(string res_path, int dig = 10)
        {
            std::ofstream input_x_data;
            input_x_data.open(res_path + "/x.txt", std::ios::out);
            input_x_data << std::fixed << std::showpoint << std::setprecision(dig);
            for (int i = 0; i < Nx; ++i) {
                input_x_data << x[i] << '\n';
            }
            input_x_data.close();
        }

        void save_t(string res_path, int dig = 10)
        {
            std::ofstream input_t_data;
            input_t_data.open(res_path + "/t_data.txt", std::ios::out);
            input_t_data << std::fixed << std::showpoint << std::setprecision(dig);
            for (int i = 0; i < Ndata; ++i) {
                input_t_data << t_data[i] << '\n';
            }
            input_t_data.close();
        }

        void save_C_MV(string res_path, int dig = 10)
        {
            std::ofstream input_Data_C_MV;
            input_Data_C_MV.open(res_path + "/Data_C_MV.txt", std::ios::out);
            input_Data_C_MV << std::fixed << std::showpoint << std::setprecision(dig);
            for (int i = 0; i < Ndata; ++i) {
                for (int j = 0; j < Nx-1; ++j) {
                    input_Data_C_MV << Data_C_MV[i][j] << ' ';
                }
                input_Data_C_MV << Data_C_MV[i][Nx-1] << '\n';
            }
            input_Data_C_MV.close();
        }

        void save_C_OV(string res_path, int dig = 10)
        {
            std::ofstream input_Data_C_OV;
            input_Data_C_OV.open(res_path + "/Data_C_OV.txt", std::ios::out);
            input_Data_C_OV << std::fixed << std::showpoint << std::setprecision(dig);
            for (int i = 0; i < Ndata; ++i) {
                for (int j = 0; j < Nx - 1; ++j) {
                    input_Data_C_OV << Data_C_OV[i][j] << ' ';
                }
                input_Data_C_OV << Data_C_OV[i][Nx-1] << '\n';
            }
            input_Data_C_OV.close();
        }

        void save_phi(string res_path, int dig = 10)
        {
            std::ofstream input_Data_phi;
            input_Data_phi.open(res_path + "/Data_phi.txt", std::ios::out);
            input_Data_phi << std::fixed << std::showpoint << std::setprecision(dig);
            for (int i = 0; i < Ndata; ++i) {
                for (int j = 0; j < Nx - 1; ++j) {
                    input_Data_phi << Data_phi[i][j] << ' ';
                }
                input_Data_phi << Data_phi[i][Nx-1] << '\n';
            }
            input_Data_phi.close();
        }

        void save_E(string res_path, int dig = 10)
        {
            std::ofstream input_Data_E;
            input_Data_E.open(res_path + "/Data_E.txt", std::ios::out);
            input_Data_E << std::fixed << std::showpoint << std::setprecision(dig);
            for (int i = 0; i < Ndata; ++i) {
                for (int j = 0; j < Nx - 2; ++j) {
                    input_Data_E << Data_E[i][j] << ' ';
                }
                input_Data_E << Data_E[i][Nx - 2] << '\n';
            }
            input_Data_E.close();
        }

        void save_k2(string res_path, int dig = 10)
        {
            std::ofstream input_Data_k2;
            input_Data_k2.open(res_path + "/Data_k2.txt", std::ios::out);
            input_Data_k2 << std::fixed << std::showpoint << std::setprecision(dig);
            for (int i = 0; i < Ndata; ++i) {
                input_Data_k2 << Data_k2[i] << '\n';
            }
            input_Data_k2.close();
        }

        ~Film() {}

        int get_Nx() { return Nx; }
        void set_Nx(int a) { Nx = a; ++check_init; }

        int get_Nt() { return Nt; }
        void set_Nt(int a) { Nt = a; }

        int get_Ndata() { return Ndata; }
        void set_Ndata(int a) { Ndata = a; }

        int get_n_save() { return n_save; }
        void set_n_save(int a) { n_save = a; ++check_init; }

        double get_h() { return h; }
        void set_h(double a) { h = a; }

        double get_dt() { return dt; }
        void set_dt(double a) { dt = a; }

        double get_phi_ext() { return phi_ext; }
        void set_phi_ext(double a) { phi_ext = a; ++check_init; }

        double get_L() { return L; }
        void set_L(double a) { L = a; ++check_init; }

        double get_T() { return T; }
        void set_T(double a) { T = a; ++check_init; }

        double get_D_MV() { return D_MV; }
        void set_D_MV(double a) { D_MV = a; ++check_init; }

        double get_D_OV() { return D_OV; }
        void set_D_OV(double a) { D_OV = a; ++check_init; }

        double get_k1_0() { return k1_0; }
        void set_k1_0(double a) { k1_0 = a; ++check_init; }

        double get_k2_0() { return k2_0; }
        void set_k2_0(double a) { k2_0 = a; ++check_init; }

        double get_k3_0() { return k3_0; }
        void set_k3_0(double a) { k3_0 = a; ++check_init; }

        double get_k4_0() { return k4_0; }
        void set_k4_0(double a) { k4_0 = a; ++check_init; }

        double get_k5_0() { return k5_0; }
        void set_k5_0(double a) { k5_0 = a; ++check_init; }

        double get_an1() { return an1; }
        void set_an1(double a) { an1 = a; ++check_init; }

        double get_an2() { return an2; }
        void set_an2(double a) { an2 = a; ++check_init; }

        double get_e_f() { return e_f; }
        void set_e_f(double a) { e_f = a; ++check_init; }

        double get_e_dl() { return e_dl; }
        void set_e_dl(double a) { e_dl = a; ++check_init; }

        double get_e_cdl() { return e_cdl; }
        void set_e_cdl(double a) { e_cdl = a; ++check_init; }

        double get_d_dl() { return d_dl; }
        void set_d_dl(double a) { d_dl = a; ++check_init; }

        double get_d_cdl() { return d_cdl; }
        void set_d_cdl(double a) { d_cdl = a; ++check_init; }

        double get_Temp() { return Temp; }
        void set_Temp(double a) { Temp = a; ++check_init; }

        double get_Ff() { return Ff; }
        void set_Ff(double a) { Ff = a; ++check_init; }

        double get_e_0() { return e_0; }
        void set_e_0(double a) { e_0 = a; ++check_init; }

        double get_R() { return R; }
        void set_R(double a) { R = a; ++check_init; }

        double get_A_k() { return A_k; }
        void set_A_k(double a) { A_k = a; ++check_init; }

        double get_A_D() { return A_D; }
        void set_A_D(double a) { A_D = a; ++check_init; }

        double get_A_L() { return A_L; }
        void set_A_L(double a) { A_L = a; ++check_init; }

        double get_A_phi() { return A_phi; }
        void set_A_phi(double a) { A_phi = a; ++check_init; }

        double get_A_C() { return A_C; }

        double get_A_t() { return A_t; }

        vector<double> get_C_MV() { return C_MV; }
        void set_C_MV(vector<double> & a) { C_MV = a; }

        vector<double> get_C_OV() { return C_OV; }
        void set_C_OV(vector<double> & a) { C_OV = a; }

};

#endif // FILM_H