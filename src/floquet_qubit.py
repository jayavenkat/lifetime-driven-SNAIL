#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 10:02:40 2019

@author: xiaoxuisaac
"""

import numpy as np
import qutip
import pickle
import time
import os
import floquet_v2
#from ist_new import IstSemiclassicalCavity
#from ist_new import IstQuantumCavity 
from ist_new import IstPump, TransmonPump, SnailPump
qubit_list = [IstPump, TransmonPump, SnailPump]


data_dir = './floquet_data/'
PERIOD_MESH =501

qubit_dict = {}
for q in qubit_list:
    qubit_dict[q.__name__] = q

class Floquet1dQubit():
    def __init__(self, static_paras = {}, Qubit_Type = TransmonPump):
        '''
        TODO: need to implement change basis matrix to phi
        TODO: need to construct Wigner function representation
        TODO: time evolution stuffs
        TODO: saving figures?
        '''
        self.propagator_map = {}
        self.floquet_modes_map = {}
        self.quasienergies_map = {}
        self.quantum_number_map = {}
        if 'class_name' in static_paras.keys():
            Qubit_Type = qubit_dict[static_paras['class_name']]
        
        self.static_qubit = Qubit_Type(static_paras)
        self.static_paras = self.static_qubit.paras
        self.change_basis_matrix = self.static_qubit.change_basis_matrix
    
    @classmethod
    def fromdir(cls, directory, dummy = False, fast = True):
        if directory[-1]!='/': directory += '/'
        paras = pickle.load( open( directory + 'paras.dat', "rb" ) )
        obj = cls(paras)
        obj.directory = directory
        if not dummy:
            obj.quasienergies_map = obj.load_data_map('quasienergies_map')
            try:
                obj.quantum_number_map = obj.load_data_map('quantum_number_map')
            except:
                pass
            if not fast:
                obj.propagator_map = obj.load_data_map('propagator_map')
                obj.floquet_modes_map = obj.load_data_map('floquet_modes_map')
        return obj
    
    @classmethod
    def fromdir_sub(cls, directory, nbar_min, nbar_max, nbar_step, fast = False):
        if directory[1]!='/': directory += '/'
        paras = pickle.load( open( directory + 'paras.dat', "rb" ) )
        obj = cls(paras)
        obj.directory = directory
        obj.quasienergies_map = obj.load_data_map_sub('quasienergies_map', nbar_min, nbar_max, nbar_step)
        try:
            obj.quantum_number_map = obj.load_data_map_sub('quantum_number_map', nbar_min, nbar_max, nbar_step)
        except:
            pass
        if not fast:
            obj.propagator_map = obj.load_data_map_sub('propagator_map', nbar_min, nbar_max, nbar_step)
            obj.floquet_modes_map = obj.load_data_map_sub('floquet_modes_map', nbar_min, nbar_max, nbar_step)
        
        return obj
        
    
    def export_data(self, name = None):
        if name == None:
            name = time.strftime("%d%b%Y-%H%M%S", time.localtime())
        directory = data_dir+name + '/'
        if 'directory' in self.__dict__.keys(): directory = self.directory
        self.directory = directory
        
        #create directory if not exists
        if not os.path.exists(directory):
            os.makedirs(directory)

        #export parameters in strings
        with open(directory + 'info.txt', 'w+') as f:
            if os.stat(directory + 'info.txt').st_size==0:
                f.write(self.static_qubit.export_paras_string())
                
        #export serealized parameters
        pickle.dump(self.static_qubit.paras, open(directory+'paras.dat', "wb" ))
        self.export_data_map(self.propagator_map, 'propagator_map')
        self.export_data_map(self.floquet_modes_map, 'floquet_modes_map')
        self.export_data_map(self.quasienergies_map, 'quasienergies_map')
        self.export_data_map(self.quantum_number_map, 'quantum_number_map')
        
        
    def export_data_map(self, dict_map, name):
        directory = self.directory + name +'/'
        if not os.path.exists(directory): os.makedirs(directory)   
        for key in dict_map.keys():
            index = pickle.loads(key)
            file_name = self._get_file_name(index['omega_p'], index['nbar'])
            qutip.qsave(dict_map[key], directory + file_name)
        
    def load_data_map(self, name):
        directory = self.directory + name +'/'
        data_map = {}
        for f in os.listdir(directory):
            file_name =  os.path.splitext(f)[0]
            print (file_name)
            omega_p = float(file_name.split(' ')[1])
            nbar = float(file_name.split(' ')[3])
            index = self._get_index(omega_p, nbar)
            data_map[index] = qutip.qload(directory+file_name)
        return data_map
    
 
    def load_data_map_sub(self, name, nbar_min, nbar_max, nbar_step):
        directory = self.directory + name +'/'
        data_map = {}
        flag = 0
        for f in os.listdir(directory):
            file_name =  os.path.splitext(f)[0]
            omega_p = float(file_name.split(' ')[1])
            nbar = float(file_name.split(' ')[3])
            if nbar not in np.arange(nbar_min, nbar_max+nbar_step, nbar_step):
                continue
            elif nbar in np.arange(nbar_min, nbar_max+nbar_step, nbar_step):
                flag = 1
                index = self._get_index(omega_p, nbar)
                data_map[index] = qutip.qload(directory+file_name)
        if flag == 0:
            print ("the nbar mesh entered hasn't been swept")
        return data_map
                                        
    def _get_index(self, omega_p, nbar):
        nbar = float(nbar)
        omega_p = float(omega_p)
        return pickle.dumps({"omega_p":omega_p, "nbar":nbar})
    
    def _get_file_name(self, omega_p, nbar):
        return 'omega_p ' + str(float(omega_p)) +' nbar '+str(float(nbar)) 
    
    def _get_hamiltonian(self, omega_p, nbar):

        return self.static_qubit.hamiltonian_builder(omega_p, nbar)

    def _load_single_file(self, map_name, omega_p, nbar):
        try: 
            directory = self.directory + map_name +'/'
            file_name = self._get_file_name(omega_p, nbar) 
            data_map = getattr(self, map_name)
            index = self._get_index(omega_p, nbar)
            data_map[index] = qutip.qload(directory+file_name)
#            setattr(self, map_name, data_map)
            return True
        except:
            return False
    
    def _get_propagator(self, omega_p, nbar):
        index = self._get_index(omega_p, nbar)
        H, args = self._get_hamiltonian(omega_p, nbar)
        T = 2 * np.pi / omega_p
        U = qutip.propagator(H, np.linspace(0, T, PERIOD_MESH), [], args)[-1]
        self.propagator_map[index] = U
        return U
    
    def propagator(self, omega_p, nbar):
        index = self._get_index(omega_p, nbar)
        if index not in self.propagator_map.keys() \
            and not self._load_single_file('propagator_map', omega_p, nbar):
            return self._get_propagator(omega_p, nbar)
        return self.propagator_map[index]


    def _get_floquet_mode_energy(self, omega_p, nbar):
        index = self._get_index(omega_p, nbar)
        H, args = self._get_hamiltonian(omega_p, nbar)
        U = self.propagator(omega_p, nbar)
        T = 2 * np.pi / omega_p
        f_modes, quasi_energies = floquet_v2.floquet_modes(H, T, args, sort = True, U=U)
        self.floquet_modes_map[index] = f_modes
        self.quasienergies_map[index] = quasi_energies
        
    
    def floquet_modes(self,omega_p, nbar):
        index = self._get_index(omega_p, nbar)
        if index not in self.floquet_modes_map.keys()\
            and not self._load_single_file('floquet_modes_map', omega_p, nbar):
            self._get_floquet_mode_energy(omega_p, nbar)
        return self.floquet_modes_map[index]

    def quasienergies(self, omega_p, nbar):
        index = self._get_index(omega_p, nbar)
        if index not in self.quasienergies_map.keys():
            self._get_floquet_mode_energy(omega_p, nbar)
        return self.quasienergies_map[index]

    def floquet_mode_charge_basis(self,floquet_mode):
        return self.change_basis_matrix * floquet_mode
    

    def floquet_mode_phase_basis(self, charge_state, \
                phi = np.arange(-2*np.pi, 2*np.pi, 4*np.pi*0.0025)):
        change_basis_matrix_phi = self.static_qubit._change_basis_matrix_phi(phi)
        return change_basis_matrix_phi * charge_state

    def _count_nodes(self, state):
        charge_state = self.floquet_mode_charge_basis(state)
        phi = np.arange(-np.pi, np.pi, 2*np.pi*0.0025)
        phase_state = self.floquet_mode_phase_basis(charge_state, phi).full().flatten()
    
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
#        avg = np.max(np.abs(array))/8
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
    
    def count_nodes(self, omega_p, nbar, state_n):
        index = self._get_index(omega_p, nbar)
#        state = self.floquet_modes(omega_p, nbar)[state_n]
#        return self._count_nodes(state)
        if index not in self.quantum_number_map.keys():
            self.quantum_number_map[index] = {}
        if state_n not in self.quantum_number_map[index].keys():
            state = self.floquet_modes(omega_p, nbar)[state_n]
            qn = self._count_nodes(state)
            self.quantum_number_map[index][state_n] = qn
        return self.quantum_number_map[index][state_n]
    

    def quantum_number(self, omega_p, nbar, state_n):
        index = self._get_index(omega_p, nbar)
        if index not in self.quantum_number_map.keys():
            self.quantum_number_map[index] = {}
        if state_n not in self.quantum_number_map[index].keys():
            state = self.floquet_modes(omega_p, nbar)[state_n]
            qn = self.static_qubit.quantum_number(state)
            self.quantum_number_map[index][state_n] = qn
        return self.quantum_number_map[index][state_n]
    
        
            
            
        
        