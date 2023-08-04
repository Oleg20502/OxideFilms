# distutils: language = c++
# cython: language_level=3

from libcpp.vector cimport vector
cimport numpy as cnp
import numpy as np
from cfilm cimport CppFilm

cdef class Film:
    cdef CppFilm* thisptr

    def __cinit__(self):
        self.thisptr = new CppFilm()

    def __dealloc__(self):
        del self.thisptr

    def init(self):
        return self.thisptr.init()

    def solve(self):
        return self.thisptr.solve()
    
    def save_x(self, res_path, dig = 10):
        self.thisptr.save_x(res_path.encode(), dig)

    def save_t(self, res_path, dig = 10):
        self.thisptr.save_t(res_path.encode(), dig)
    
    def save_C_MV(self, res_path, dig = 10):
        self.thisptr.save_C_MV(res_path.encode(), dig)

    def save_C_OV(self, res_path, dig = 10):
        self.thisptr.save_C_OV(res_path.encode(), dig)
    
    def save_phi(self, res_path, dig = 10):
        self.thisptr.save_phi(res_path.encode(), dig)

    def save_E(self, res_path, dig = 10):
        self.thisptr.save_E(res_path.encode(), dig)
    
    def save_k2(self, res_path, dig = 10):
        self.thisptr.save_k2(res_path.encode(), dig)
    

    @property
    def Nx(self):
        return self.thisptr.get_Nx()
    
    @Nx.setter
    def Nx(self, a):
        self.thisptr.set_Nx(a)


    @property
    def Nt(self):
        return self.thisptr.get_Nt()
    
    @Nt.setter
    def Nt(self, a):
        self.thisptr.set_Nt(a)
    

    @property
    def Ndata(self):
        return self.thisptr.get_Ndata()
    
    @Ndata.setter
    def Ndata(self, a):
        self.thisptr.set_Ndata(a)
    

    @property
    def n_save(self):
        return self.thisptr.get_n_save()
    
    @n_save.setter
    def n_save(self, a):
        self.thisptr.set_n_save(a)
    

    @property
    def h(self):
        return self.thisptr.get_h()
    
    @h.setter
    def h(self, a):
        self.thisptr.set_h(a)
    

    @property
    def dt(self):
        return self.thisptr.get_dt()
    
    @dt.setter
    def dt(self, a):
        self.thisptr.set_dt(a)
    

    @property
    def phi_ext(self):
        return self.thisptr.get_phi_ext()
    
    @phi_ext.setter
    def phi_ext(self, a):
        self.thisptr.set_phi_ext(a)
    

    @property
    def L(self):
        return self.thisptr.get_L()
    
    @L.setter
    def L(self, a):
        self.thisptr.set_L(a)
    

    @property
    def T(self):
        return self.thisptr.get_T()
    
    @T.setter
    def T(self, a):
        self.thisptr.set_T(a)
    

    @property
    def D_MV(self):
        return self.thisptr.get_D_MV()
    
    @D_MV.setter
    def D_MV(self, a):
        self.thisptr.set_D_MV(a)
    

    @property
    def D_OV(self):
        return self.thisptr.get_D_OV()
    
    @D_OV.setter
    def D_OV(self, a):
        self.thisptr.set_D_OV(a)
    

    @property
    def k1_0(self):
        return self.thisptr.get_k1_0()
    
    @k1_0.setter
    def k1_0(self, a):
        self.thisptr.set_k1_0(a)
    

    @property
    def k2_0(self):
        return self.thisptr.get_k2_0()
    
    @k2_0.setter
    def k2_0(self, a):
        self.thisptr.set_k2_0(a)
    

    @property
    def k3_0(self):
        return self.thisptr.get_k3_0()
    
    @k3_0.setter
    def k3_0(self, a):
        self.thisptr.set_k3_0(a)
    

    @property
    def k4_0(self):
        return self.thisptr.get_k4_0()
    
    @k4_0.setter
    def k4_0(self, a):
        self.thisptr.set_k4_0(a)
    

    @property
    def k5_0(self):
        return self.thisptr.get_k5_0()
    
    @k5_0.setter
    def k5_0(self, a):
        self.thisptr.set_k5_0(a)
    

    @property
    def an1(self):
        return self.thisptr.get_an1()
    
    @an1.setter
    def an1(self, a):
        self.thisptr.set_an1(a)


    @property
    def an2(self):
        return self.thisptr.get_an2()
    
    @an2.setter
    def an2(self, a):
        self.thisptr.set_an2(a)
    

    @property
    def e_f(self):
        return self.thisptr.get_e_f()
    
    @e_f.setter
    def e_f(self, a):
        self.thisptr.set_e_f(a)
    

    @property
    def e_dl(self):
        return self.thisptr.get_e_dl()
    
    @e_dl.setter
    def e_dl(self, a):
        self.thisptr.set_e_dl(a)
    

    @property
    def e_cdl(self):
        return self.thisptr.get_e_cdl()
    
    @e_cdl.setter
    def e_cdl(self, a):
        self.thisptr.set_e_cdl(a)
    

    @property
    def d_dl(self):
        return self.thisptr.get_d_dl()
    
    @d_dl.setter
    def d_dl(self, a):
        self.thisptr.set_d_dl(a)
    

    @property
    def d_cdl(self):
        return self.thisptr.get_d_cdl()
    
    @d_cdl.setter
    def d_cdl(self, a):
        self.thisptr.set_d_cdl(a)
    

    @property
    def Temp(self):
        return self.thisptr.get_Temp()
    
    @Temp.setter
    def Temp(self, a):
        self.thisptr.set_Temp(a)
    

    @property
    def Ff(self):
        return self.thisptr.get_Ff()
    
    @Ff.setter
    def Ff(self, a):
        self.thisptr.set_Ff(a)
    

    @property
    def e_0(self):
        return self.thisptr.get_e_0()
    
    @e_0.setter
    def e_0(self, a):
        self.thisptr.set_e_0(a)
    

    @property
    def R(self):
        return self.thisptr.get_R()
    
    @R.setter
    def R(self, a):
        self.thisptr.set_R(a)
    

    @property
    def A_k(self):
        return self.thisptr.get_A_k()
    
    @A_k.setter
    def A_k(self, a):
        self.thisptr.set_A_k(a)
    

    @property
    def A_D(self):
        return self.thisptr.get_A_D()
    
    @A_D.setter
    def A_D(self, a):
        self.thisptr.set_A_D(a)
    

    @property
    def A_L(self):
        return self.thisptr.get_A_L()
    
    @A_L.setter
    def A_L(self, a):
        self.thisptr.set_A_L(a)


    @property
    def A_phi(self):
        return self.thisptr.get_A_phi()

    @A_phi.setter
    def A_phi(self, a):
        self.thisptr.set_A_phi(a)

    
    @property
    def A_C(self):
        return self.thisptr.get_A_C()

    @property
    def A_t(self):
        return self.thisptr.get_A_t()

    
    @property
    def C_MV(self):
        cdef int Nx = self.thisptr.get_Nx()
        return np.asarray(<cnp.float_t[:Nx]> self.thisptr.get_C_MV().data())
    
    @C_MV.setter
    def C_MV(self, a):
        Nx = a.size()
        if Nx != self.Nx:
            print(f'WARNING: Size {Nx} of a is not equal to size {self.Nx} of the class')
        cdef vector[double] arr
        cdef int i
        for i in range(Nx):
            arr.push_back(a[i])
        self.thisptr.set_C_MV(arr)


    
    
        
    
    


    
