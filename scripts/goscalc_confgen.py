import numpy as np
from matplotlib import pyplot as plt
import mendeleev as md
from panta_rhei.gui.eels_atlas_data import eels_atlas
import json

edgedict   =   {'K' :[1,0,'K1',1], 'K1' :[1,0,'K1' ,1],
                'L1':[2,0,'L1',1], 'L2,3':[2,1,'L23',1], 'L2':[2,1,'L23',1/3], 'L3':[2,1,'L23',2/3],
                'M1':[3,0,'M1',1], 'M2,3':[3,1,'M23',1], 'M2':[3,1,'M23',1/3], 'M3':[3,1,'M23',2/3], 'M4,5':[3,2,'M45',1], 'M4':[3,2,'M45',2/5], 'M5':[3,2,'M45',3/5],
                'N1':[4,0,'N1',1], 'N2,3':[4,1,'N23',1], 'N2':[4,1,'N23',1/3], 'N3':[4,1,'N23',2/3], 'N4,5':[4,2,'N45',1], 'N4':[4,2,'N45',2/5], 'N5':[4,2,'N45',3/5],
                                   'N6,7':[4,3,'N67',1], 'N6':[4,3,'N67',3/7], 'N7':[4,3,'N67',2/3],
                'O1':[5,0,'O1',1], 'O2,3':[5,1,'O23',1], 'O2':[5,1,'O23',1/3], 'O3':[5,1,'O23',2/3], 'O4,5':[5,2,'O45',1], 'O4':[5,2,'O45',2/5], 'O5':[5,2,'O45',3/5]}


ldict = {'s':0, 'p':1, 'd':2, 'f':3}

def long_econf(Z):
    c0 = md.element(Z).econf
    
    if ']' in c0:
        c = c0
        cl = []
        while ']'in c:
            i=c.find(']')
            cl.append(c[i+1:])
            sym = c[1:i]
            c = md.element(sym).econf
        cl.append(c)
    else:
        cl = [c0]
    cl.reverse()
    
    conf = cl[0]
    
    for i in range(1, len(cl)):
        conf += cl[i]
    return conf

def split_subshells(econf):
    if ' ' in econf:
        subshells = []
        while ' ' in econf:
            i = econf.find(' ')
            ssh = econf[0:i]
            subshells.append(ssh)
            econf = econf[i+1:]
        subshells.append(econf)
    else:
        subshells = [econf]
    return subshells

def unpack_subshell(ssh):
    ssh_arr = np.zeros(4, dtype='int')
    ssh_arr[0] = ssh[0]
    ssh_arr[1] = ldict[ssh[1]]
    if len(ssh) == 4:
        ne = int(ssh[2:4])
    elif len(ssh) == 3:
        ne = int(ssh[2])
    else:
        ne = 1
    ssh_arr[2] = np.floor(ne/2.)
    ssh_arr[3] = np.ceil(ne/2.)
    return ssh_arr
    
def subshells_to_table(subshells):
    table = np.zeros([len(subshells), 4], dtype='int')
    i = 0
    for ssh in subshells:
        table[i, :] = unpack_subshell(ssh)
        i += 1
    return table

def electron_table(Z):
    econf = long_econf(Z)
    subshells = split_subshells(econf)
    table = subshells_to_table(subshells)
    return table

def order_table(table):
    ordered_table = table[table[:, 1].argsort()]  # sort by l
    ordered_table = ordered_table[ordered_table[:, 0].argsort(kind='mergesort')]  # sort by n
    return ordered_table


def create_wavegen(Z, outpath='wavegen.dat', dfa='LDA'):
    table = electron_table(Z)
    table = order_table(table)
#    string = ['     {}     {}     {}. {}.']
    string = '{:6d}{:6d}{:6d}. {:d}.\n'
    with open(outpath, 'w') as f:
        f.write(dfa+'\n')
        f.write(' {:.8f}\n'.format(Z))
        for i in range(table.shape[0]):
            f.write(string.format(table[i,0], table[i,1], table[i,2], table[i,3]))
        
def collect_names():
    names = []
    for Z in range(1,110):
        sym = md.element(Z).symbol
        element = eels_atlas[sym]
        for e in element.edges:
            if e.name not in names:
                names.append(e.name)
    return names

def create_config_json(edge):
    sym = edge.element
    Z = md.element(sym).atomic_number
    n, l, name, _ = edgedict[edge.name]
    dE = edge.energy
    e_start = 0.2
    e_steps = 96
    e_inc = 1.1
    settings = { "dft_filename": "waveup.dat",
                "output_dir_name": 'outputs/'+sym+'_'+name,
                "n_bound": n,
                "l_bound": l,
                "max_considered_lfree": 15,
                "energy_free_start": e_start,
                "energy_free_steps": e_steps,
                "energy_free_increase": e_inc,
                "max_kvalue_Ang": 70}
    with open('config.json', 'w') as f:
        json.dump(settings, f, indent = 4 )
    
