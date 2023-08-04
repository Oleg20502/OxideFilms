# distutils: language = c++
# cython: language_level=3

from cparams cimport CppParams

cdef class Params:
    cdef CppParams* thisptr

    def __cinit__(self):
        self.thisptr = new CppParams()

    def __dealloc__(self):
        del self.thisptr

    @property
    def Nx(self):
        return self.thisptr.Nx

    @Nx.setter
    def Nx(self, a):
        self.thisptr.Nx = a


    @property
    def Nt(self):
        return self.thisptr.Nt

    @Nt.setter
    def Nt(self, a):
        self.thisptr.Nt = a


    @property
    def Ndata(self):
        return self.thisptr.Ndata

    @Ndata.setter
    def Ndata(self, a):
        self.thisptr.Ndata = a


    @property
    def h(self):
        return self.thisptr.h

    @h.setter
    def h(self, a):
        self.thisptr.h = a