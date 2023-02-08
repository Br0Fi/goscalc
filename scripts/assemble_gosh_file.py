import numpy as np
from os import path
import os
import h5py as h5
import json

edgenames = ['K', 'K1', 'L1', 'L2,3', 'L23', 'L2', 'L3',
                        'M1', 'M2,3', 'M23', 'M3', 'M2', 'M4,5', 'M45', 'M4', 'M5', 
                        'N1', 'N2,3', 'N23', 'N3', 'N2', 'N4,5', 'N45', 'N4', 'N5', 'N6,7', 'N67', 'N6', 'N7', 
                        'O1', 'O2,3', 'O23', 'O3', 'O2', 'O4,5', 'O45', 'O4', 'O5']

edgedict   =   {'K' :[1,0,'K1',1], 'K1' :[1,0,'K1' ,1],
                'L1':[2,0,'L1',1], 'L23':[2,1,'L23',1], 'L2,3':[2,1,'L23',1], 'L2':[2,1,'L23',1/3], 'L3':[2,1,'L23',2/3],
                'M1':[3,0,'M1',1], 'M23':[3,1,'M23',1], 'M2,3':[3,1,'M23',1], 'M2':[3,1,'M23',1/3], 'M3':[3,1,'M23',2/3], 'M45':[3,2,'M45',1], 'M4,5':[3,2,'M45',1], 'M4':[3,2,'M45',2/5], 'M5':[3,2,'M45',3/5],
                'N1':[4,0,'N1',1], 'N23':[4,1,'N23',1], 'N2,3':[4,1,'N23',1], 'N2':[4,1,'N23',1/3], 'N3':[4,1,'N23',2/3], 'N45':[4,2,'N45',1], 'N4,5':[4,2,'N45',1], 'N4':[4,2,'N45',2/5], 'N5':[4,2,'N45',3/5],
                                   'N67':[4,3,'N67',1], 'N6,7':[4,3,'N67',1], 'N6':[4,3,'N67',3/7], 'N7':[4,3,'N67',4/7],
                'O1':[5,0,'O1',1], 'O23':[5,1,'O23',1], 'O2,3':[5,1,'O23',1], 'O2':[5,1,'O23',1/3], 'O3':[5,1,'O23',2/3], 'O45':[5,2,'O45',1], 'O4,5':[5,2,'O45',1], 'O4':[5,2,'O45',2/5], 'O5':[5,2,'O45',3/5]}

k0=2*np.pi/1.97e-12

def load_dft_GOS_dat(gos_folder):
    qaxis = np.loadtxt(path.join(gos_folder,"k.dat"))*1e10 # convert from inverse angstroms to inverse meters
    eaxis = np.loadtxt(path.join(gos_folder,"free_energies.dat"))
    gos = np.loadtxt(path.join(gos_folder,"gos.dat"))
    with open(path.join(gos_folder,"gos.dat")) as f:
        f.readline()
        l = f.readline()
        eonset = l.split('\t')[0]
        eonset = eonset[1:-2]
        eonset = float(eonset)-0.2
    return gos, qaxis, eaxis, eonset

def downsample_gos(gos_folder, nq=192, ne=None):
    gos, qaxis, eaxis, eonset = load_dft_GOS_dat(gos_folder)
    if ne == None :
        ne = len(eaxis)
    qaxis = qaxis[qaxis<= 0.2*k0]
    if nq == None :
        nq = len(qaxis)
    q_samples = np.linspace(0,len(qaxis)-1, nq).round().astype('int')
    e_samples = np.linspace(0,len(eaxis)-1, ne).round().astype('int')
    return gos[q_samples,:][:, e_samples], qaxis[q_samples], eaxis[e_samples], eonset

def nest_metadata(group, k, v):
    if type(v) == dict :
        gi = group.create_group(k)
        for ki, vi in v.items():
            nest_metadata(gi, ki, vi)
    else:
            if type(v) == int :
                group.attrs[k] = np.int32(v)
            elif type(v) == float :
                group.attrs[k] = np.float32(v)
            else : 
                group.attrs[k] = v

def downsample_all(folder_in, file_out):
    contents = os.listdir(folder_in)
    edges_in = [folder for folder in contents if path.isdir(os.path.join(folder_in, folder)) ]
    #edges_in = edges_in[0::15]
    elements = [ edge[:edge.find('_')] for edge in edges_in ]
    edges = [ edgedict[edge[1+edge.find('_'):] ][2] for edge in edges_in ]
    
    folders_in  = [path.join(folder_in, edge)  for edge in edges_in]

    with open(path.join(folder_in,'metadata.json')) as f:
        metadata = json.load(f)
    
    with h5.File(file_out, 'w') as h:

        #dt = np.dtype([ ('edge', 'S4'),('ref_edge', 'S4' ),('normalisation', 'f4'), ('n', 'i2'),( 'l', 'i2')] )
        #table = h.create_dataset('edges_table', 38, dtype=dt)
        #for i,n in enumerate(edgenames):
        #    ep = edgedict[n]
        #    table[i] = n, ep[2], ep[3], ep[0], ep[1]
        #print("added edges notation table")
        
        h.attrs['file_format'] = 'GOSH'
        h.attrs['file_format_version'] = '0.7.2'
        group_metadata = h.create_group('metadata')
        for k, v in metadata.items():
            nest_metadata(group_metadata, k, v)

        print("added metadata")

        for i, f in enumerate(folders_in):
            gos, qaxis, eaxis, eonset = downsample_gos(f)
            gos = gos[:,:, np.newaxis] # make it 3D as required
            edge_group_name = '{}/{}'.format(elements[i],edges[i])
            edge_group = h.create_group(edge_group_name)
            gos_data = edge_group.create_dataset('data',
                                        shape = gos.shape,
                                        data = gos,
                                        dtype = 'f',
                                        compression = 'gzip')
            q_data = edge_group.create_dataset('q',
                                        shape = qaxis.shape,
                                        data = qaxis.astype('f'),
                                        dtype = 'f',
                                        compression = 'gzip')
            e_data = edge_group.create_dataset('free_energies',
                                        shape = eaxis.shape,
                                        data = eaxis.astype('f'),
                                        dtype = 'f',
                                        compression = 'gzip')
            v_data = edge_group.create_dataset('variants',
                                        shape = [1],
                                        data = np.asarray(['default'], dtype='S128'),
                                        dtype = 'S128',
                                        compression = 'gzip')
            edge_metadata = edge_group.create_group('metadata')
            edge_metadata.attrs['onset'] = eonset
            print("downsampled {}, completed:{}/{}".format(edge_group_name, i, len(folders_in)))
