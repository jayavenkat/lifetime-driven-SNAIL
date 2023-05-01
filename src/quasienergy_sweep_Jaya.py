#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 20:25:47 2019

@author: jv437
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pickle
from floquet_qubit import Floquet1dQubit
from wcs_colors import ccycle_floquet as ccycle
from plot_settings import *
#%%
PERIOD_MESH =501
N_max_charge = 80
N_max_b = 20

if 1:
    windows = [
               np.linspace(0, 1, 1)
              ]
    omega_p = 12
    
    static_paras = {'qubit_type': 'snail with classical cavity and quantum qubit',
                    'quantum_modes': 1,
                    'class_name': 'SnailPump',
                    'Ej': 3500*1.5*1.4,
                    'Ec': 0.0014,
                    'NoJ' : 3,
                    'alpha' : 0.1,
                    'Ng' : 0.0,
                    'phi_ext': 0.2,
                    'N_max_charge': 300,
                    'N_max_b': 20}
    
    # static_paras = {'Ej' : 90,
    #            'phi_ext': 0.33*3,
    #            'Ec' : 0.060,
    #            'NoJ' : 3,
    #            'alpha' : 0.1,
    #            'Ng' : 0.0,
    #            "N_max_charge":N_max_charge, 
    #            "N_max_b":N_max_b,
    #            "quantum_modes":1
    #            }
    #floquet_qubit = Floquet1dQubit.static_paras
    
    floquet_qubit = Floquet1dQubit(static_paras)
    print(floquet_qubit.static_qubit.qubit_report())
    
    floquet_qubit = Floquet1dQubit(static_paras)
#%%
    #%%
    '''
    Setting up sweep
    '''
    if 1:
        windows = [np.linspace(0.0, 0.05, 51)]
        windows = [np.linspace(0.05, 0.5, 51)]
        windows = [np.linspace(0.0, 1.0, 11)]
        static_paras = {'qubit_type': 'snail with classical cavity and quantum qubit',
                        'quantum_modes': 1,
                        'class_name': 'SnailPump',
                        'Ej': 320,
                        'Ec': 0.033,
                        'NoJ' : 3,
                        'alpha' : 0.1,
                        'Ng' : 0.0,
                        'phi_ext': 0.3,
                        'N_max_charge': 75,
                        'N_max_b': 20}
        static_paras = {'qubit_type': 'snail with classical cavity and quantum qubit',
                        'quantum_modes': 1,
                        'class_name': 'SnailPump',
                        'Ej': 7350,
                        'Ec': 0.0014,
                        'NoJ' : 3,
                        'alpha' : 0.1,
                        'Ng' : 0.0,
                        'phi_ext': 0.2,
                        'N_max_charge': 300,
                        'N_max_b': 50}
        floquet_qubit = Floquet1dQubit(static_paras)
        print(floquet_qubit.static_qubit.qubit_report())
        count = 0
        for window in windows:
            for nbar in window:
                floquet_qubit.quasienergies(omega_p, nbar)
                floquet_qubit.export_data('Snail_KC_v4')
                count += 1

        floquet_qubit.export_data('Snail_KC_v4')

#%%
'''
checking states and wave functions
'''
if 1:
    nbar = 0 #0.008#60
    fq = floquet_qubit
    qes = fq.quasienergies(omega_p, nbar)
    # #    omega_p = 8.9
    floquet_mode_index = 15#24
    # floquet_qubit = fq
    state = floquet_qubit.floquet_modes(omega_p, nbar)[floquet_mode_index]
    charge_state = floquet_qubit.floquet_mode_charge_basis(state)
    
    for i in range(N_max_b):
        if (fq.count_nodes(omega_p, nbar, i) == 0 or fq.count_nodes(omega_p, nbar, i)== 1):# or fq.count_nodes(omega_p, nbar, i) == 2):
            print(i, fq.count_nodes(omega_p, nbar, i), fq.quasienergies(omega_p, nbar)[i])
    phi = np.arange(-1*np.pi, 1*np.pi, 2*np.pi*0.00025)
    
    phase_state = floquet_qubit.floquet_mode_phase_basis(charge_state, phi)
    gradient = np.gradient(np.abs(phase_state.full().flatten()))
    #    print n_nodes(state)
    
    fig, ax = plt.subplots()
    #    ax.plot(phi/2/np.pi, np.abs(phase_state.full().flatten()))
    #    ax.plot(phi/2/np.pi, gradient)
    ax.plot(phi/2/np.pi, np.abs(phase_state.full().flatten()))
    ax.plot(phi/2/np.pi, 1- np.cos(phi))
    ax.axvline(x=0)
    ax.set_xlabel(r'$ \varphi/2 \pi$')
    ax.set_ylabel(r'$|\langle \varphi | \Phi_g (0) \rangle |^2$')
    title1 =  r'$| {{{}}} \rangle, $'.format((int(fq.count_nodes(omega_p, nbar, floquet_mode_index))))
    title2 = r'$nbar = {{{}}}, $'.format((nbar))
    title3 = r'$ N_g = {{{}}} $'.format((int(fq.static_paras['Ng'])))
    title4 = r'$, qe = {{{}}} $'.format((fq.quasienergies(omega_p, nbar)[floquet_mode_index]))
    plt.suptitle(title1 + title2 + title3+ title4)
#%%
if 1:
    fq = Floquet1dQubit.fromdir("./floquet_data/Snail_KC_v4/")
#%%
#%%
'''
Assigning g
'''
if 1: # calibarate nbar to starkshift
    qe_m = fq.quasienergies_map
    ge = []
    nbars = []
    gs = []
    es = []
    g_dict = {}
    e_dict = {}
    for index in qe_m:
        info = pickle.loads(index)
        nbar = info['nbar']
        omega_p = info['omega_p']
        g = None
        e = None
        
        if 1  and omega_p == 12:
            for i, qe in enumerate(qe_m[index]):
                nn = fq.count_nodes(omega_p, nbar, i)
                if nn == 0 or nn == 1 and fq.quasienergies(omega_p, nbar)[i] < 0: 
                    print(nbar, i, qe)
                    g = qe
                if nn == 1 or nn == 0 and fq.quasienergies(omega_p, nbar)[i] > 0.0: 
                    print(nbar, i, qe)
                    e = qe
                        
            if g != None and e!= None and (e-g):
                nbars.append(nbar)
                ge.append(e-g)
                gs.append(g)
                es.append(e)
                
                g_dict[nbar] = g
                e_dict[nbar] = e
    fq.export_data()
    nbars_sort, gs_sort = (list(t) for t in zip(*sorted(zip(nbars, gs))))
    nbars2_sort, es_sort = (list(t) for t in zip(*sorted(zip(nbars, es))))
#%%
fig, ax = plt.subplots()

ax.plot(np.array(nbars_sort), (np.array(ge)), '.', color = 'black')
# ax.set_ylim(4.5, 5.1)
ax.set_xlabel('nbar')
ax.set_ylabel(r'$\Delta_{\rm{ac}}/2\pi (\rm{GHz})$')
# ax.set_ylim(6.015, 6.022)
set_size(4, 2, ax)    
#%%
if 1:
    fig, ax = plt.subplots()
    qe_m = fq.quasienergies_map
    labels = []
    windows = [[[0.0, 1.0], [[-omega_p, omega_p]]]]    
    
    nbars = []
#    ges = []
    for index in qe_m:
        info = pickle.loads(index)
        nbar = info['nbar']
        omega_p = info['omega_p']
        
#        ss = ge_fit(nbar) - ge_fit(0)
#        ss = ss / -0.0375
#        x = ss
        if 1 :
            
            for window in windows:
                # if nbar >= window[0][0] and nbar < window[0][1] :
                    nbars.append(nbar)
                    for i, qe in enumerate(qe_m[index]):
                        if 1: #i%10 == 0:
                            hit = False
                            qe = (qe)% (omega_p/2)
    #                        qe = (qe - g_dict[nbar])% (omega_p)
    #                        qe = (qe - g_fit(nbar)-omega_p/4.0) % omega_p 
    #                        qe = qe - (qe > omega_p/2)*omega_p 
                            for freq in window[1]:
                                
                                if 1: #qe>= freq[0] and qe < freq[1]:
                                    hit = True
                                    break
                            if hit:                        
                                ax.scatter(nbar, qe , c = 'grey', s = 1.2, alpha = 0.5)
    #                            print(nbar, 'scattered')
            #                    nn = n_nodes(fq.floquet_modes(omega_p, nbar)[i])
                                nn = fq.count_nodes(omega_p, nbar, i)
    #                            if nn < 12 and nn >= 0:
                                
                                if nn in [0, 1]:#, 4, 5, 6]:#, 2, 3]:
                                    if nn not in labels:
                                        labels.append(nn)
                                    else:
                                        ax.scatter(nbar, qe , s = 8, label = r'$| {{{}}} \rangle$'.format((int(nn))), color = states[int(nn)])
            
    # xlim = ax.get_xbound()
    # ylim = ax.get_ybound()
    
    labels.sort()
    # for nn in labels:
    #     ax.scatter(-100, 0, c = [ccycle[int(nn)]], s = 30, label = '|' + str(nn)+'>')
#coun
    # ax.set_xlim(xlim)
    # ax.set_ylim(2, 2.5)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.04,1.04), loc="upper left")
#    ax.legend(bbox_to_anchor=(1.04,1.04), loc="upper left")
    ax.set_xlabel(r'$\bar{n}$')
    ax.set_ylabel(r'$\epsilon_n$')
    title = 'quasienergies vs nbar for Ng = ' + str(fq.static_paras['Ng'])
    plt.suptitle(title)
#    ax.set_xlim(0.0, 0.15)