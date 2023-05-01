#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 22:03:37 2018

@author: jv437
"""

# coding: utf-8
# Do not import unicode_literals or hamiltonian spec breaks in Python 2.
import numpy as np
import qutip
import pickle


phi_stepsize = 0.00025


N_max_charge = 150
N_max_b = 110
N_max_charge = 80
N_max_b = 20
N_max_a = 4

# p_snail = {'Ej' : 17.1,
#            'phi_ext': 0.0,
#            'Ec' : 0.095,
#            'NoJ' : 3,
#            'alpha' : 0.1,
#            'Ng' : 0.0,
#            "N_max_charge":N_max_charge, 
#            "N_max_b":N_max_b,
#            "quantum_modes":1
#            }

# p_transmon = {'Ej' : 17.1,
#            'Ec' : 0.095,
#            'Ng' : 0.0,
#            "N_max_charge":N_max_charge, 
#            "N_max_b":N_max_b,
#            "quantum_modes":1
#            }

# p_snail = {'Ej' : 17.1,
#            'phi_ext': 0.33,
#            'Ec' : 0.095,
#            'NoJ' : 3,
#            'alpha' : 0.1,
#            'Ng' : 0.0,
#            "N_max_charge":N_max_charge, 
#            "N_max_b":N_max_b,
#            "quantum_modes":1
#            }
class QubitPump:
    def __init__(self):
        pass
    
    def export_paras_string(self):
        return str(self.paras).replace(',', ',\n')
    
    def serialize_paras(self):
        return pickle.dumps(self.paras)
        
    @classmethod
    def fromparas(cls, paras):
        paras = pickle.load(paras)
        return cls(paras)
        
    


    def phi_basis_converter(self, phi_value):
        N_max_charge = self.N_max_charge
        phi_ket = np.array([0.0+0.0j]*(2*N_max_charge + 1))
        for k in range(0, 2*N_max_charge+1):
            phi_ket += np.asarray(np.exp(1j*(k - N_max_charge) * phi_value) * \
                                  qutip.fock(2*N_max_charge+1, k).full().flatten())
        return phi_ket
            
    def _change_basis_matrix_phi(self, \
                    phi = np.arange(-np.pi, np.pi, 2*np.pi*phi_stepsize)):
        if 'phi' not in self.__dict__.keys() or not np.array_equal(self.phi, phi):
                phi_ket_list = []
                for i in range(0, len(phi)):
                    phi_ket_list.append(self.phi_basis_converter(phi[i]))
                self.change_basis_matrix_phi = qutip.Qobj(np.column_stack([x for x in \
                                            phi_ket_list])).dag()
                self.phi =phi
        return self.change_basis_matrix_phi
    
    def qubit_report(self):
        evals = self.evals_ch
        ge = evals[1] - evals[0]
        ef = evals[2] - evals[1]
        print("qubit frequency: %.4f GHz\nanharmonicity: %.4f MHz"%(ge, (ef-ge)*1e3))
        return ge, (ef-ge)
class TransmonPump(QubitPump):
    """
    transmon class takes in parameters. Calculates various operators associated 
    to it, evals, evecs
    """
    def __init__(self, static_paras = {}):
        self.qubit_type = 'transmon with classical cavity and quantum qubit'
        self.quantum_modes = 1.5
        if static_paras == None:
            static_paras = {}

        for key in p_transmon:
            print(key)
            if key in static_paras.keys():
                setattr(self, key, static_paras[key])
            else:
                setattr(self, key, p_transmon[key])
            
        
        self.paras = {
         'qubit_type': self.qubit_type,
         'quantum_modes' : self.quantum_modes,
         'class_name': self.__class__.__name__,}
        for key in p_transmon:
            self.paras[key] = getattr(self,key)
            
        self.export_paras_string()
        self.serialize_paras()
        self.N_constructor()
        self.phi_constructor()
        self.H_charge_constructor()
        self.evecbasis_converter()
        self.H_eig_constructor()
        self.hamiltonian_static_builder()
        

    def N_constructor(self):
        N_max_charge = self.N_max_charge
        self.N_small = qutip.Qobj(\
                         np.diag(range(-N_max_charge, N_max_charge + 1)) \
                         ) #  notes 
        
    def phi_constructor(self):
        N_max_charge = self.N_max_charge

        self.cos_phi = qutip.Qobj(\
                                 np.diag(0.5 * np.ones(2 * N_max_charge), k = 1)\
                                 + \
                                 np.diag(0.5 * np.ones(2 * N_max_charge), k = -1) \
                                 )  #  cos (\varphi) 
        self.sin_phi = qutip.Qobj(\
                                 np.diag(-0.5j * np.ones(2 * N_max_charge), k = 1)\
                                 + \
                                 np.diag(0.5j * np.ones(2 * N_max_charge), k = -1) \
                                 )  #  sin (phi)
        self.Id_ch = qutip.qeye(2 * N_max_charge + 1)  #  identity
        
    def H_charge_constructor(self):     
        Ec = self.Ec
        Ej = self.Ej
        Ng = self.Ng
        N_small = self.N_small
        Id_ch = self.Id_ch
        cos_phi = self.cos_phi
        self.H_charge = qutip.Qobj(\
                         4 * Ec * (N_small  -Ng*Id_ch )**2 - \
                         Ej * cos_phi
                         ) 
        self.evals_ch, self.evecs_ch = self.H_charge.eigenstates()
        
    def H_eig_constructor(self):
        Ec = self.Ec
        Ng = self.Ng
        N_small_eigv_basis = self.N_small_eigv_basis
        Id_eigv_basis = self.Id_eigv_basis
        
        self.H_0_eig_kinetic = (\
                 4.0 * Ec \
                 * (N_small_eigv_basis - Ng * Id_eigv_basis )**2 \
                 )                                      
    
    def hamiltonian_static_builder(self):
        Ej = self.Ej
        cos_phi_eigv_basis = self.cos_phi_eigv_basis
        sin_phi_eigv_basis = self.sin_phi_eigv_basis

        self.h_0 = self.H_0_eig_kinetic
        self.h_1 = -1.0 * Ej* cos_phi_eigv_basis
        self.h_1_coef = 'cos(oscillating_prefactor * sin(omega_p * t))'
        self.h_2 = Ej * sin_phi_eigv_basis
        self.h_2_coef = 'sin(oscillating_prefactor * sin(omega_p * t))'
        # Qutip hamiltonian object
        
        self.h = [self.h_0, [self.h_1, self.h_1_coef], [self.h_2, self.h_2_coef]]



    def evecbasis_converter(self):
        N_small = self.N_small
        evecs_ch = self.evecs_ch 
        cos_phi = self.cos_phi
        sin_phi = self.sin_phi
        N_max_b = self.N_max_b

        # Matrix to pass from charge states basis to qubit eigenstate basis. 
        cbm = qutip.Qobj(\
                                     np.column_stack(x.full() \
                                                     for x in evecs_ch) \
                                     )
        # Changing basis to qubit eigenstates 
        self.N_small_eigv_basis = cbm.dag() * N_small * cbm
        self.cos_phi_eigv_basis = cbm.dag() * cos_phi * cbm
        self.sin_phi_eigv_basis = cbm.dag() * sin_phi * cbm
        self.Id_eigv_basis = qutip.qeye(N_max_b)
        self.N_small_eigv_basis = qutip.Qobj(self.N_small_eigv_basis[:N_max_b, :N_max_b])
        self.cos_phi_eigv_basis = qutip.Qobj(self.cos_phi_eigv_basis[:N_max_b, :N_max_b])
        self.sin_phi_eigv_basis = qutip.Qobj(self.sin_phi_eigv_basis[:N_max_b, :N_max_b])
        self.Id_eigv_basis = qutip.qeye(N_max_b)
        self.change_basis_matrix = qutip.Qobj(\
                                     np.column_stack(x.full() \
                                                     for x in evecs_ch[:N_max_b]) \
                                     )        

    def compute_epsilon_p(self, omega_p, nbar):
        omega_q = np.sqrt(8*self.Ec*self.Ej)
        
        epsilon_p = np.sqrt(nbar) * (omega_p**2-omega_q**2)/(2*omega_p)\
                    *2**(9.0/4)*(self.Ec/self.Ej)**0.25
        return epsilon_p
    
    def hamiltonian_builder(self, omega_p, nbar):
        
        epsilon_p = self.compute_epsilon_p(omega_p, nbar)

        oscillating_prefactor = epsilon_p / omega_p

        args = {'oscillating_prefactor': oscillating_prefactor, \
                'omega_p': omega_p} 
        return self.h, args
          


class SnailPump(QubitPump):
    """
    snail class takes in parameters. Calculates various operators associated 
    to it, evals, evecs
    """
    def __init__(self, static_paras = {}):
        self.qubit_type = 'snail with classical cavity and quantum qubit'
        
        self.quantum_modes = 1
        if static_paras == None:
            static_paras = {}

        for key in static_paras:
            print(key)
            if key in static_paras.keys():
                setattr(self, key, static_paras[key])
            else:
                setattr(self, key, static_paras[key])
            
        
        self.paras = {
         'qubit_type': self.qubit_type,
         'quantum_modes' : self.quantum_modes,
         'class_name': self.__class__.__name__,}
        for key in static_paras:
            self.paras[key] = getattr(self,key)
            
        self.alpha = static_paras['alpha']
        self.phi_ext = static_paras['phi_ext']
        self.export_paras_string()
        self.serialize_paras()
        self.N_constructor()
        self.phi_constructor()
        self.H_charge_constructor()
        self.evecbasis_converter()
        self.H_eig_constructor()
        self.hamiltonian_static_builder()
        
        

    def N_constructor(self):
        alpha = self.alpha
        NoJ = self.NoJ
        phi_ext = self.phi_ext
        N_max_charge = self.N_max_charge
        self.N_small = qutip.Qobj(\
                         np.diag(range(-N_max_charge, N_max_charge + 1))/ \
                         NoJ \
                         ) #  notes   
        self.N_large = qutip.Qobj(\
                         np.diag(range(-N_max_charge, N_max_charge + 1)) \
                         ) #  notes                   
    

    def phi_constructor(self):
        NoJ = self.NoJ
        alpha = self.alpha
        N_max_charge = self.N_max_charge

        self.cos_phi_overNoJ = qutip.Qobj(\
                                 np.diag(0.5 * np.ones(2 * N_max_charge), k = 1)\
                                 + \
                                 np.diag(0.5 * np.ones(2 * N_max_charge), k = -1) \
                                 )  #  cos (\varphi/M) 
        self.sin_phi_overNoJ = qutip.Qobj(\
                                 np.diag(-0.5j * np.ones(2 * N_max_charge), k = 1)\
                                 + \
                                 np.diag(0.5j * np.ones(2 * N_max_charge), k = -1) \
                                 )  #  sin (\varphi/M)
        self.cos_phi = qutip.Qobj(\
                         np.diag(0.5 * np.ones(2 * N_max_charge -(NoJ-1)), \
                                 k = NoJ) + \
                         np.diag(0.5 * np.ones(2 * N_max_charge -(NoJ-1)), \
                                 k = -NoJ) \
                                 )  #  cos(\varphi) 
        self.sin_phi = qutip.Qobj(\
                         np.diag(-0.5j * np.ones(2 * N_max_charge - (NoJ - 1)), \
                                 k = NoJ) + \
                         np.diag(0.5j * np.ones(2 * N_max_charge - (NoJ - 1)), \
                                 k = - NoJ) \
                                 )   #  sin (\varphi)
        self.Id_ch = qutip.qeye(2 * N_max_charge + 1)  #  identity

    
    def H_charge_constructor(self):     
        Ec = self.Ec
        Ej = self.Ej # Ej of one big junction
        NoJ = self.NoJ
        alpha = self.alpha
        phi_ext = self.phi_ext
        Ng = self.Ng
        N_small = self.N_small
        Id_ch = self.Id_ch
        cos_phi = self.cos_phi
        cos_phi_overNoJ = self.cos_phi_overNoJ
        sin_phi_overNoJ = self.sin_phi_overNoJ
        self.H_charge = qutip.Qobj(\
                         4 * Ec * (N_small  -Ng*Id_ch )**2 - \
                         alpha * Ej * cos_phi - \
                         NoJ  * Ej * cos_phi_overNoJ * np.cos(phi_ext/NoJ) - \
                         NoJ * Ej * sin_phi_overNoJ * np.sin(phi_ext/NoJ)
                         ) 
        self.evals_ch, self.evecs_ch = self.H_charge.eigenstates()

    def H_eig_constructor(self):
        Ec = self.Ec
        Ng = self.Ng
        alpha = self.alpha
        N_small_eigv_basis = self.N_small_eigv_basis
        Id_eigv_basis = self.Id_eigv_basis
        
        self.H_0_eig_kinetic = (\
                 4.0 * Ec \
                 * (N_small_eigv_basis - Ng * Id_eigv_basis )**2 \
                 ) 
        
    def hamiltonian_static_builder(self):
        Ej = self.Ej
        cos_phi_eigv_basis = self.cos_phi_eigv_basis
        sin_phi_eigv_basis = self.sin_phi_eigv_basis
        cos_phi_overNoJ_eigv_basis = self.cos_phi_overNoJ_eigv_basis
        sin_phi_overNoJ_eigv_basis = self.sin_phi_overNoJ_eigv_basis
        NoJ = self.NoJ
        alpha = self.alpha
        phi_ext = self.phi_ext

        self.h_0 = self.H_0_eig_kinetic
        self.h_1 = -1.0 * alpha * Ej* cos_phi_eigv_basis
        self.h_1_coef = 'cos(oscillating_prefactor * sin(omega_p * t))'
        self.h_2 = alpha * Ej * sin_phi_eigv_basis
        self.h_2_coef = 'sin(oscillating_prefactor * sin(omega_p * t))'
        self.h_3 = - 1.0 * NoJ * Ej * cos_phi_overNoJ_eigv_basis * np.cos(phi_ext/NoJ)
        self.h_3_coef = 'cos(1.0/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        self.h_4 = 1.0 * NoJ * Ej * sin_phi_overNoJ_eigv_basis * np.cos(phi_ext/NoJ)
        self.h_4_coef = 'sin(1.0/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        self.h_5 = - 1.0 * NoJ * Ej * sin_phi_overNoJ_eigv_basis * np.sin(phi_ext/NoJ)
        self.h_5_coef = 'cos(1.0/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        self.h_6 = - 1.0 * NoJ * Ej * cos_phi_overNoJ_eigv_basis * np.sin(phi_ext/NoJ)
        self.h_6_coef = 'sin(1.0/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        # Qutip hamiltonian object
        
        self.h = [self.h_0, [self.h_1, self.h_1_coef], [self.h_2, self.h_2_coef],\
                  [self.h_3, self.h_3_coef], [self.h_4, self.h_4_coef], \
                      [self.h_5, self.h_5_coef], [self.h_6, self.h_6_coef]]
        

    def evecbasis_converter(self):
        N_small = self.N_small
        N_large = self.N_large
        evecs_ch = self.evecs_ch 
        cos_phi = self.cos_phi
        sin_phi = self.sin_phi
        cos_phi_overNoJ = self.cos_phi_overNoJ
        sin_phi_overNoJ = self.sin_phi_overNoJ
        N_max_b = self.N_max_b

        # Matrix to pass from charge states basis to qubit eigenstate basis. 
        cbm = qutip.Qobj(\
                                     np.column_stack(x.full() \
                                                     for x in evecs_ch) \
                                     )
        # Changing basis to qubit eigenstates 
        self.N_small_eigv_basis = cbm.dag() * N_small * cbm
        self.cos_phi_eigv_basis = cbm.dag() * cos_phi * cbm
        self.cos_phi_overNoJ_eigv_basis = cbm.dag() * cos_phi_overNoJ * cbm
        self.sin_phi_eigv_basis = cbm.dag() * sin_phi * cbm
        self.sin_phi_overNoJ_eigv_basis = cbm.dag() * sin_phi_overNoJ * cbm
        self.N_large_eigv_basis = cbm.dag() * N_large* cbm
        self.Id_eigv_basis = qutip.qeye(N_max_b)
        self.N_small_eigv_basis = qutip.Qobj(self.N_small_eigv_basis[:N_max_b, :N_max_b])
        self.N_large_eigv_basis = qutip.Qobj(self.N_large_eigv_basis[:N_max_b, :N_max_b])
        self.cos_phi_eigv_basis = qutip.Qobj(self.cos_phi_eigv_basis[:N_max_b, :N_max_b])
        self.sin_phi_eigv_basis = qutip.Qobj(self.sin_phi_eigv_basis[:N_max_b, :N_max_b])
        self.cos_phi_overNoJ_eigv_basis = qutip.Qobj(self.cos_phi_overNoJ_eigv_basis[:N_max_b, :N_max_b])
        self.sin_phi_overNoJ_eigv_basis = qutip.Qobj(self.sin_phi_overNoJ_eigv_basis[:N_max_b, :N_max_b])
        self.Id_eigv_basis = qutip.qeye(N_max_b)
        self.change_basis_matrix = qutip.Qobj(\
                                     np.column_stack(x.full() \
                                                     for x in evecs_ch[:N_max_b]) \
                                     )
    

    def compute_epsilon_p(self, omega_p, nbar):
        El = self.Ej/(self.NoJ*self.alpha)
        omega_q = np.sqrt(8*self.Ec*(self.Ej+El))
        
        delta = omega_p - omega_q
        alpha = (self.Ej + El/self.NoJ**2)/(self.Ej + El)*self.Ec
#        print (delta+alpha)/delta, delta, alpha, omega_q
        nbar = nbar*(delta+alpha)/delta
        
        epsilon_p = np.sqrt(nbar) * (omega_p**2-omega_q**2)/(2*omega_p)\
                    *2**(9.0/4)*(self.Ec/(self.Ej+El))**0.25
        return epsilon_p
    
    def hamiltonian_builder(self, omega_p, nbar):
        
        epsilon_p = self.compute_epsilon_p(omega_p, nbar)

        oscillating_prefactor = epsilon_p / omega_p

        args = {'oscillating_prefactor': oscillating_prefactor, \
                'omega_p': omega_p, 'NoJ': self.NoJ} 
        return self.h, args



class IstPump(QubitPump):
    """
    ist class takes in parameters. Calculates various operators associated 
    to it, evals, evecs
    """
    def __init__(self, static_paras = {}):
        self.qubit_type = 'ist with classical cavity and quantum qubit'
        self.quantum_modes = 1
        if static_paras == None:
            static_paras = {}

        for key in p_ist:
            print(key)
            if key in static_paras.keys():
                setattr(self, key, static_paras[key])
            else:
                setattr(self, key, p_ist[key])
            
        
        self.paras = {
         'qubit_type': self.qubit_type,
         'quantum_modes' : self.quantum_modes,
         'class_name': self.__class__.__name__,}
        for key in p_ist:
            self.paras[key] = getattr(self,key)
            
            
        self.export_paras_string()
        self.serialize_paras()
        self.N_constructor()
        self.phi_constructor()
        self.H_charge_constructor()
        self.evecbasis_converter()
        self.H_eig_constructor()
        self.hamiltonian_static_builder()
        

    def N_constructor(self):
        NoJ = self.NoJ
        N_max_charge = self.N_max_charge
        self.N_small = qutip.Qobj(\
                         np.diag(range(-N_max_charge, N_max_charge + 1))/ \
                         NoJ \
                         ) #  notes   
        self.N_large = qutip.Qobj(\
                         np.diag(range(-N_max_charge, N_max_charge + 1)) \
                         ) #  notes                   
    

    def phi_constructor(self):
        NoJ = self.NoJ
        N_max_charge = self.N_max_charge

        self.cos_phi_overNoJ = qutip.Qobj(\
                                 np.diag(0.5 * np.ones(2 * N_max_charge), k = 1)\
                                 + \
                                 np.diag(0.5 * np.ones(2 * N_max_charge), k = -1) \
                                 )  #  cos (\varphi/M) 
        self.sin_phi_overNoJ = qutip.Qobj(\
                                 np.diag(-0.5j * np.ones(2 * N_max_charge), k = 1)\
                                 + \
                                 np.diag(0.5j * np.ones(2 * N_max_charge), k = -1) \
                                 )  #  sin (\varphi/M)
        self.cos_phi = qutip.Qobj(\
                         np.diag(0.5 * np.ones(2 * N_max_charge -(NoJ-1)), \
                                 k = NoJ) + \
                         np.diag(0.5 * np.ones(2 * N_max_charge -(NoJ-1)), \
                                 k = -NoJ) \
                                 )  #  cos(\varphi) 
        self.sin_phi = qutip.Qobj(\
                         np.diag(-0.5j * np.ones(2 * N_max_charge - (NoJ - 1)), \
                                 k = NoJ) + \
                         np.diag(0.5j * np.ones(2 * N_max_charge - (NoJ - 1)), \
                                 k = - NoJ) \
                                 )   #  sin (\varphi)
        self.Id_ch = qutip.qeye(2 * N_max_charge + 1)  #  identity

    
    def H_charge_constructor(self):     
        Ec = self.Ec
        Ej = self.Ej
        NoJ = self.NoJ
        beta = self.beta
        Ng = self.Ng
        N_small = self.N_small
        Id_ch = self.Id_ch
        cos_phi = self.cos_phi
        cos_phi_overNoJ = self.cos_phi_overNoJ
        self.H_charge = qutip.Qobj(\
                         4 * Ec * (N_small  -Ng*Id_ch )**2 - \
                         Ej * cos_phi - \
                         NoJ * beta * Ej * cos_phi_overNoJ 
                         ) 
        self.evals_ch, self.evecs_ch = self.H_charge.eigenstates()


    def H_eig_constructor(self):
        Ec = self.Ec
        Ng = self.Ng
        N_small_eigv_basis = self.N_small_eigv_basis
        Id_eigv_basis = self.Id_eigv_basis
        
        self.H_0_eig_kinetic = (\
                 4.0 * Ec \
                 * (N_small_eigv_basis - Ng * Id_eigv_basis )**2 \
                 ) 
        
    def hamiltonian_static_builder(self):
        Ej = self.Ej
        cos_phi_eigv_basis = self.cos_phi_eigv_basis
        sin_phi_eigv_basis = self.sin_phi_eigv_basis
        cos_phi_overNoJ_eigv_basis = self.cos_phi_overNoJ_eigv_basis
        sin_phi_overNoJ_eigv_basis = self.sin_phi_overNoJ_eigv_basis
        NoJ = self.NoJ
        beta = self.beta

        self.h_0 = self.H_0_eig_kinetic
        self.h_1 = -1.0 * Ej* cos_phi_eigv_basis
        self.h_1_coef = 'cos(oscillating_prefactor * sin(omega_p * t))'
        self.h_2 = Ej * sin_phi_eigv_basis
        self.h_2_coef = 'sin(oscillating_prefactor * sin(omega_p * t))'
        self.h_3 = - 1.0 * NoJ * beta * Ej * cos_phi_overNoJ_eigv_basis 
        self.h_3_coef = 'cos(1.0/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        self.h_4 = 1.0 * NoJ * beta * Ej * sin_phi_overNoJ_eigv_basis
        self.h_4_coef = 'sin(1.0/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        # Qutip hamiltonian object
        
        self.h = [self.h_0, [self.h_1, self.h_1_coef], [self.h_2, self.h_2_coef],\
                  [self.h_3, self.h_3_coef], [self.h_4, self.h_4_coef]]
        

    def evecbasis_converter(self):
        N_small = self.N_small
        N_large = self.N_large
        evecs_ch = self.evecs_ch 
        cos_phi = self.cos_phi
        sin_phi = self.sin_phi
        cos_phi_overNoJ = self.cos_phi_overNoJ
        sin_phi_overNoJ = self.sin_phi_overNoJ
        N_max_b = self.N_max_b

        # Matrix to pass from charge states basis to qubit eigenstate basis. 
        cbm = qutip.Qobj(\
                                     np.column_stack(x.full() \
                                                     for x in evecs_ch) \
                                     )
        # Changing basis to qubit eigenstates 
        self.N_small_eigv_basis = cbm.dag() * N_small * cbm
        self.cos_phi_eigv_basis = cbm.dag() * cos_phi * cbm
        self.cos_phi_overNoJ_eigv_basis = cbm.dag() * cos_phi_overNoJ * cbm
        self.sin_phi_eigv_basis = cbm.dag() * sin_phi * cbm
        self.sin_phi_overNoJ_eigv_basis = cbm.dag() * sin_phi_overNoJ * cbm
        self.N_large_eigv_basis = cbm.dag() * N_large* cbm
        self.Id_eigv_basis = qutip.qeye(N_max_b)
        self.N_small_eigv_basis = qutip.Qobj(self.N_small_eigv_basis[:N_max_b, :N_max_b])
        self.N_large_eigv_basis = qutip.Qobj(self.N_large_eigv_basis[:N_max_b, :N_max_b])
        self.cos_phi_eigv_basis = qutip.Qobj(self.cos_phi_eigv_basis[:N_max_b, :N_max_b])
        self.sin_phi_eigv_basis = qutip.Qobj(self.sin_phi_eigv_basis[:N_max_b, :N_max_b])
        self.cos_phi_overNoJ_eigv_basis = qutip.Qobj(self.cos_phi_overNoJ_eigv_basis[:N_max_b, :N_max_b])
        self.sin_phi_overNoJ_eigv_basis = qutip.Qobj(self.sin_phi_overNoJ_eigv_basis[:N_max_b, :N_max_b])
        self.Id_eigv_basis = qutip.qeye(N_max_b)
        self.change_basis_matrix = qutip.Qobj(\
                                     np.column_stack(x.full() \
                                                     for x in evecs_ch[:N_max_b]) \
                                     )
    

    def compute_epsilon_p(self, omega_p, nbar):
        El = self.Ej*self.beta/self.NoJ
        omega_q = np.sqrt(8*self.Ec*(self.Ej+El))
        
        delta = omega_p - omega_q
        alpha = (self.Ej + El/self.NoJ**2)/(self.Ej + El)*self.Ec
#        print (delta+alpha)/delta, delta, alpha, omega_q
        nbar = nbar*(delta+alpha)/delta
        
        epsilon_p = np.sqrt(nbar) * (omega_p**2-omega_q**2)/(2*omega_p)\
                    *2**(9.0/4)*(self.Ec/(self.Ej+El))**0.25
        return epsilon_p
    
    def hamiltonian_builder(self, omega_p, nbar):
        
        epsilon_p = self.compute_epsilon_p(omega_p, nbar)

        oscillating_prefactor = epsilon_p / omega_p

        args = {'oscillating_prefactor': oscillating_prefactor, \
                'omega_p': omega_p, 'NoJ': self.NoJ} 
        return self.h, args
    
class IstQuantumCavityPump(IstPump):
    """
    ist class takes in parameters. Calculates various operators associated 
    to it, evals, evecs
    """
    def __init__(self, static_paras = {}):
        self.qubit_type = 'ist with quantum cavity and quantum qubit'
        self.quantum_modes = 2.0 #cavity and qubit
        if static_paras == None:
            static_paras = {}

        for key in p_cavity_ist:
            print(key)
            if key in static_paras.keys():
                setattr(self, key, static_paras[key])
            else:
                setattr(self, key, p_cavity_ist[key])
            
        
        self.paras = {
         'qubit_type': self.qubit_type,
         'quantum_modes' : self.quantum_modes,
         'class_name': self.__class__.__name__,}
        for key in p_cavity_ist:
            self.paras[key] = getattr(self,key)
            
        self.export_paras_string()
        self.serialize_paras()
        self.N_constructor()
        self.phi_constructor()
        self.H_charge_constructor()
        self.evecbasis_converter()
        self.H_eig_constructor()
        self.hamiltonian_static_builder()
        
    
    def hamiltonian_static_builder(self):
        N_small_eigv_basis = self.N_small_eigv_basis
        N_large_eigv_basis = self.N_large_eigv_basis
        cos_phi_eigv_basis = self.cos_phi_eigv_basis
        sin_phi_eigv_basis = self.sin_phi_eigv_basis
        cos_phi_overNoJ_eigv_basis = self.cos_phi_overNoJ_eigv_basis
        sin_phi_overNoJ_eigv_basis = self.sin_phi_overNoJ_eigv_basis
        Id_eigv_basis = self.Id_eigv_basis
        NoJ = self.NoJ
        beta = self.beta
        N_a = self.N_max_a
        omega_a = self.omega_a
        Eg = self.Eg
        Ej = self.Ej
        Ec = self.Ec
        Ng = self.Ng
        
        a = qutip.destroy(N_a)
        X_a = (a + a.dag())/ np.sqrt(2)  
        P_a = (a - a.dag()) / 1.0j / np.sqrt(2)     
        Id_a = qutip.qeye(N_a)
        
        self.N_a_tensor = qutip.tensor(a.dag() * a, Id_eigv_basis)
        self.P_a_tensor = qutip.tensor(P_a, Id_eigv_basis)
        self.X_a_tensor = qutip.tensor(X_a, Id_eigv_basis)
        self.N_small_tensor = qutip.tensor(Id_a, N_small_eigv_basis)
        self.N_large_tensor = qutip.tensor(Id_a, N_large_eigv_basis)
        self.cos_phi_tensor = qutip.tensor(Id_a, cos_phi_eigv_basis)
        self.cos_phi_overNoJ_tensor=qutip.tensor(Id_a, cos_phi_overNoJ_eigv_basis)
        self.sin_phi_tensor = qutip.tensor(Id_a, sin_phi_eigv_basis)
        self.sin_phi_overNoJ_tensor=qutip.tensor(Id_a, sin_phi_overNoJ_eigv_basis)
        self.Id_tensor = qutip.tensor(Id_a, Id_eigv_basis)
        self.h_cav = omega_a * self.N_a_tensor
        self.h_coupling = Eg * (self.N_small_tensor - Ng * self.Id_tensor) * self.P_a_tensor * np.sqrt(2)
        self.h_qub_kinetic = 4 * Ec * (self.N_small_tensor - Ng * self.Id_tensor)**2
        self.h_0 = self.h_cav + self.h_qub_kinetic + self.h_coupling
        self.h_1 = -1.0 * Ej* self.cos_phi_tensor
        self.h_1_coef = 'cos(oscillating_prefactor * sin(omega_p * t))'
        self.h_2 = Ej * self.sin_phi_tensor
        self.h_2_coef = 'sin(oscillating_prefactor * sin(omega_p * t))'
        self.h_3 = - 1.0 * NoJ * beta * Ej * self.cos_phi_overNoJ_tensor
        self.h_3_coef = 'cos(1/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        self.h_4 = 1.0 * NoJ * beta * Ej * self.sin_phi_overNoJ_tensor
        self.h_4_coef = 'sin(1/ NoJ * oscillating_prefactor  * sin(omega_p * t))'
        # Qutip hamiltonian object
        
        self.h = [self.h_0, [self.h_1, self.h_1_coef], [self.h_2, self.h_2_coef],\
                  [self.h_3, self.h_3_coef], [self.h_4, self.h_4_coef]]
    
    def evecbasis2chargebasis(self, state):
        return self.change_basis_matrix * state
    
    def evecbasis2phasebasis(self, state, phi = np.arange(np.pi, np.pi, 2*np.pi*0.0025)):
        return self._change_basis_matrix_phi(phi) * self.change_basis_matrix * state
    
    def quantum_number(self, state):
        phi = np.arange(-np.pi, np.pi, 2*np.pi*0.0025)
        phase_state = self.evecbasis2phasebasis(state, phi).full().flatten()
    
        ps = np.real(phase_state)
        avg = np.max(np.abs(ps))/6
        nodess = 0
        sign_track = 0
        for i in range(len(ps)):
            if abs(ps[i]) < avg:
                ps[i] = 0
            if ps[i] * ps[i - 1] < 0:
                nodess +=1 
            if ps[i - 1] == 0 and ps[i]!=0 and ps[i]*sign_track < 0:
    #            print i, array[i], array[i - 1], sign_track
                nodess +=1
            if ps[i] != 0: sign_track = ps[i]

    
        array = np.gradient(np.abs(phase_state))
        avg = np.max(np.abs(array))/6
        nodes = 0
        sign_track = 0
        for i, a in enumerate(array):
            if abs(a) < avg:
                array[i] = 0
                
            if array[i] * array[i - 1] < 0:
                nodes +=1
            if array[i - 1] == 0 and array[i]!=0 and a*sign_track < 0:
                nodes +=1
            if array[i] != 0: sign_track = array[i]
        nodes = (nodes - 1)/2
        if abs(nodess - nodes) > 8: return -1
        return nodes


       
