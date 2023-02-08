import mendeleev as md
from panta_rhei.gui.eels_atlas_data import eels_atlas
from goscalc_confgen import *
import subprocess

skiplist = ['L2','M2','N2','O2','M4','N4', 'O4', 'N6']

def compute_element(Z):
    create_wavegen(Z)
    subprocess.call('./wavegen_mod')
    sym = md.element(Z).symbol
    edges = eels_atlas[sym].edges
    if edges != []:
        for e in edges:
            if e.name not in skiplist :
                create_config_json(e)
                subprocess.call('./goscalc')
        

for Z in range(1,99):
    compute_element(Z)
