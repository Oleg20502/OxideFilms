
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "film.hpp":
    cdef cppclass CppFilm "Film":
        CppFilm() except +
        int init(vector[double] C_MV_0, vector[double] C_OV_0)
        int init(int mode)
        int solve()

        void save_x(string res_path, int dig)
        void save_t(string res_path, int dig)
        void save_C_MV(string res_path, int dig)
        void save_C_OV(string res_path, int dig)
        void save_phi(string res_path, int dig)
        void save_E(string res_path, int dig)
        void save_k2(string res_path, int dig)

        int get_Nx()
        void set_Nx(int a)

        int get_Nt()
        void set_Nt(int a)

        int get_Ndata()
        void set_Ndata(int a)

        int get_n_save()
        void set_n_save(int a)

        double get_h()
        void set_h(double a)

        double get_dt()
        void set_dt(double a)

        double get_phi_ext()
        void set_phi_ext(double a)

        double get_L()
        void set_L(double a)

        double get_T()
        void set_T(double a)

        double get_D_MV()
        void set_D_MV(double a)

        double get_D_OV()
        void set_D_OV(double a)

        double get_k1_0()
        void set_k1_0(double a)

        double get_k2_0()
        void set_k2_0(double a)

        double get_k3_0()
        void set_k3_0(double a)

        double get_k4_0()
        void set_k4_0(double a)

        double get_k5_0()
        void set_k5_0(double a)

        double get_an1()
        void set_an1(double a)

        double get_an2()
        void set_an2(double a)

        double get_e_f()
        void set_e_f(double a)

        double get_e_dl()
        void set_e_dl(double a)

        double get_e_cdl()
        void set_e_cdl(double a)

        double get_d_dl()
        void set_d_dl(double a)

        double get_d_cdl()
        void set_d_cdl(double a)

        double get_Temp()
        void set_Temp(double a)

        double get_Ff()
        void set_Ff(double a)

        double get_e_0()
        void set_e_0(double a)

        double get_R()
        void set_R(double a)

        double get_A_k()
        void set_A_k(double a)

        double get_A_D()
        void set_A_D(double a)

        double get_A_L()
        void set_A_L(double a)

        double get_A_phi()
        void set_A_phi(double a)

        double get_A_C()

        double get_A_t()

        vector[double]& get_C_MV()
        # void set_C_MV(vector[double] a)

        vector[double]& get_C_OV()
        # void set_C_OV(vector[double] a)

        vector[double]& get_x()

        vector[double]& get_t_data()

        vector[double]& get_Data_k2()

        vector[vector[double]]& get_Data_C_MV()

        vector[vector[double]]& get_Data_C_OV()

        vector[vector[double]]& get_Data_phi()

        vector[vector[double]]& get_Data_E()
