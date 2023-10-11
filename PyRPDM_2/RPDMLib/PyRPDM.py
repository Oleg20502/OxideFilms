from scipy.optimize import newton_krylov
import numpy as np


class Params:
    def __init__(self):
        self.Nd = None
        self.Ff = None
        self.eps_0 = None
        self.R = None
        self.Temp = None

        self.an1 = None
        self.an2 = None
        self.eps_f = None
        self.eps_dl = None
        self.eps_cdl = None
        self.d_dl = None
        self.d_cdl = None

        # MV first, OV second
        self.z = None
        self.D = None

        self.k1_0 = None
        self.k2_0 = None
        self.k3_0 = None
        self.k4_0 = None
        self.k5 = None

        self.A_L = None
        self.A_J = None
        self.A_D = None
        self.A_C = None
        self.A_phi = None

        self.L = None
        self.T = None
        self.phi_ext = None

        self.x = None
        self.t = None
        self.n_save = None
    



class Film:
    def __init__(self, Params):
        
        self.Nd = Params.Nd

        self.Ff = Params.Ff #
        self.eps_0 = Params.eps_0
        self.R = Params.R #
        self.Temp = Params.Temp

        self.an1 = Params.an1
        self.an2 = Params.an2
        self.eps_f = Params.eps_f
        self.eps_dl = Params.eps_dl
        self.eps_cdl = Params.eps_cdl
        self.d_dl = Params.d_dl
        self.d_cdl = Params.d_cdl

        # MV first, OV second
        self.z = Params.z
        self.D = Params.D

        self.k1_0 = Params.k1_0
        self.k2_0 = Params.k2_0
        self.k3_0 = Params.k3_0
        self.k4_0 = Params.k4_0
        self.k5 = Params.k5 #

        self.A_L = Params.A_L
        self.A_J = Params.A_J #
        self.A_D = Params.A_D
        self.A_C = Params.A_C #
        self.A_phi = Params.A_phi #

        self.L = Params.L #
        self.T = Params.T #
        self.phi_ext = Params.phi_ext

        self.x = Params.x #
        self.t = Params.t
        self.n_save = Params.n_save


        self.i_time = 0
        self.Nx = len(self.x)
        self.Nt = len(self.t)
        self.dx = self.x[1] - self.x[0]
        #self.dt = self.t[1] - self.t[0]   #
        self.N_data = (self.Nt-2) // self.n_save + 2

        self.Data_C = np.zeros((self.N_data, self.Nd * self.Nx))
        self.Data_phi = np.zeros((self.N_data, self.Nx))
        self.Data_k2 = np.zeros(self.N_data)
        self.Data_t = np.zeros(self.N_data)


        self.C = np.zeros(self.Nd * self.Nx)
        for nd in range(self.Nd):
            for i in range(self.Nx):
                self.C[i+self.Nx*nd] = Params.C_0[i+self.Nx*nd]


        self.k1 = 0.0  #
        self.k2 = 0.0  #
        self.k3 = 0.0  #
        self.k4 = 0.0  #

        self.phi = np.zeros(self.Nx)


        self.A1 = -Params.Ff**2/self.eps_f/self.eps_0/Params.R / \
                   Params.Temp*(self.A_J*self.A_L**3/self.A_D)
        self.A2 = -self.eps_dl/self.eps_f/self.d_dl
        self.A3 = self.eps_f*(self.d_dl/self.eps_dl + self.d_cdl/self.eps_cdl)

        self.Am = np.zeros((self.Nx, self.Nx))
        self.Am[0, 0] = -1/self.dx + self.A2
        self.Am[0, 1] = 1/self.dx
        self.Am[self.Nx-1, self.Nx-1] = self.A3/self.dx + 1
        self.Am[self.Nx-1, self.Nx-2] = -self.A3/self.dx
        self.dm = np.zeros(self.Nx)
        self.dm[self.Nx-1] = 0
        self.dm[0] = self.A2 * self.phi_ext
        for i in range(1, self.Nx-1):
            self.Am[i, i] = -2/self.dx**2
            self.Am[i, i+1] = 1/self.dx**2
            self.Am[i, i-1] = 1/self.dx**2
        


    def poisson(self, C):
        for i in range(1, self.Nx-1):
            self.dm[i] = 0
            for nd in range(self.Nd):
                self.dm[i] += self.A1 * self.z[nd] * C[i+self.Nx*nd]
        return np.linalg.solve(self.Am, self.dm.T).T
    

    def update_k(self, phi_mf, phi_fs):
        k1 = self.k1_0 * np.exp(self.an1 * phi_mf)
        k2 = self.k2_0 * np.exp(self.an1 * phi_mf)
        k3 = self.k3_0 * np.exp(self.an2 * phi_fs)
        k4 = self.k4_0 * np.exp(self.an2 * phi_fs)
        return k1, k2, k3, k4


    def equation(self, C):
        dx = self.dx 
        dt = self.t[self.i_time+1] - self.t[self.i_time] # 

        phi = self.poisson(C)
        phi_mf = self.phi_ext - phi[0]
        phi_fs = phi[-1]
        k1, k2, k3, k4 = self.update_k(phi_mf, phi_fs)

        Eq = np.zeros(self.Nd * self.Nx)
        for nd in range(self.Nd):
            for i in range(1, self.Nx-1):
                k = i + nd * self.Nx

                dfi = phi[i+1] - phi[i]
                J1 = -self.D[nd] * (C[k+1] * self.B(-self.z[nd]*dfi) - C[k] * self.B(self.z[nd]*dfi)) / dx

                dfi = phi[i] - phi[i-1]
                J2 = -self.D[nd] * (C[k] * self.B(-self.z[nd]*dfi) - C[k-1] * self.B(self.z[nd]*dfi)) / dx
            
                Eq[k] = C[k] - self.C[k] + dt * (J1 - J2) / dx

        # === Boundary conditions ===
        dfi = phi[1] - phi[0]
        J0 = -self.D[0] * (C[1] * self.B(-self.z[0]*dfi) - C[0] * self.B(self.z[0]*dfi)) / dx
        Eq[0] = J0 + k1 * C[0]
        J0 = -self.D[1] * (C[1+self.Nx] * self.B(-self.z[1]*dfi) - C[self.Nx] * self.B(self.z[1]*dfi)) / dx
        Eq[self.Nx] = J0 - k2

        dfi = phi[self.Nx-1] - phi[self.Nx-2]
        J0 = -self.D[0] * (C[self.Nx-1] * self.B(-self.z[0]*dfi) - C[self.Nx-2] * self.B(self.z[0]*dfi)) / dx
        Eq[self.Nx-1] = J0 + k3
        J0 = -self.D[1] * (C[2*self.Nx-1] * self.B(-self.z[1]*dfi) - C[2*self.Nx-2] * self.B(self.z[1]*dfi)) / dx
        Eq[2*self.Nx-1] = J0 - k4 * C[2*self.Nx-1]
        # ============================
      
        return Eq
        

    def update_data(self, i_data):
        self.phi = self.poisson(self.C)
        phi_mf = self.phi_ext - self.phi[0]
        phi_fs = self.phi[-1]
        self.k1, self.k2, self.k3, self.k4 = self.update_k(phi_mf, phi_fs)

        self.Data_C[i_data] = self.C
        self.Data_phi[i_data] = self.phi
        self.Data_k2[i_data] = self.k2
        self.Data_t[i_data] = self.t[self.i_time]


    def solve(self):
        i_data = 0
        for it in range(self.Nt-1):
            self.i_time = it
            self.C = newton_krylov(self.equation, self.C)

            if (it+1) % self.n_save == 0:
                print(it+1)
                self.update_data(i_data)
                i_data += 1
        
        self.update_data(i_data)
        print('Finished!')
    
    

    def B(self, x):
        if x != 0.0:
            return x/(np.exp(x) - 1)
        else:
            return 1
                            


    










        
