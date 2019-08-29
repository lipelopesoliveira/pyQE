#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 22:23:38 2017

@author: felipe
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.lines import Line2D
from PyQt5.QtWidgets import QFileDialog, QMainWindow, QApplication
import random
import matplotlib.style
import matplotlib as mpl
from cycler import cycler
from matplotlib import colors as mcolors

#mpl.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')
#mpl.rcParams['lines.markersize'] = np.sqrt(20)

#k_path_label  = [r'X', r'$\Gamma$', r'Y$\mid$L','$\Gamma$', r'Z$\mid$N',r'$\Gamma$', r'M$\mid\Gamma$', 'R']
#k_path_label  = [r'$\Gamma$', 'X','K', r'$\Gamma$', 'L', 'X', 'W','L']
k_path_label  = [r'$\Gamma$', 'X','P', 'N', r'$\Gamma$', 'Z', 'S']

rc('font', **{'size': 16})
rc('text', usetex=True)

class Build_Files(QMainWindow):
    
    def __init__(self, prefix):
        super().__init__()
        self.prefix = prefix
        self.path = QFileDialog.getExistingDirectory()
        
    def open_gv_file(self, prefix):
        
        ''' Imput=file.gjf | Function that opens a GausView file (.gjf) and returns the positions of the atoms, 
        dimensions of the cell, number of atoms, and a list of atoms.'''
        if os.path.exists(os.path.join(self.path, prefix + '.gjf')) == True:
            gv_file = os.path.join(self.path, prefix + '.gjf')
            
        else:
            print("Cannot find GaussView file! Please select the file: ")
            gv_file = QFileDialog.getOpenFileName(self, 'Open .gjf file', self.path, "GaussView Job File (*.gjf)")[0]

        #Open file ***.gjf from GausView with atoms coordinates and cell parameters
        arq = open(gv_file, 'r') 
        file_a = arq.readlines()
        arq.close

        #Split the file in two, pos with the atoms coordinates and cel with cell parameters
        file_b = file_a[5:-1]
        pos = file_b[:-3]

        #Get cell parameters
        cel = []
        for i in file_b:
            if 'Tv' in i: cel +=[i]
        
        #Create a list of the different types of atoms
        list_type = []
        list_atoms = []
    
        for i in range(len(pos)):
            
            pos[i] = pos[i].strip()
            list_atoms+=[[pos[i][0:2].rstrip(' ')]]
            
            if (pos[i][0] in list_type) == False:
                list_type += pos[i][0]
        #Process pos and cel to put in the format for Quantum Espresso read        
        for i in range(len(cel)):
            
            cel[i] = cel[i].replace('Tv         ', ' ')  
            cel[i] = cel[i].rstrip('\n')
                
        #Get the number of different atoms
        nat = len(pos)
        #Return in list format of atoms coordinates, cell parameters, number of different atoms and a list of different atoms    
        return pos, cel, nat, list_type, list_atoms
        
    
    def read_qe_in(self, prefix):

        if os.path.exists(os.path.join(self.path, self.prefix+'.scf.in')) == True:
            qe_file = os.path.join(self.path, self.prefix+'.scf.in')
            
        else:
            print("Cannot find input file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE SCF file (*.scf.in);; QE IN file (*.in) ;;All Files (*.*)")[0]
        
        temp = open(qe_file, 'r')
        file = temp.readlines()
        temp.close()
        
        data = {}
        
        n_atom = 0

        for i in range(len(file)):
            if 'nat' in file[i]:
                n_atom = file[i].rstrip('\n')
                n_atom = n_atom.rstrip(',')
                n_atom = n_atom.lstrip(' ')
                n_atom = n_atom.lstrip('nat = ')
                n_atom = int(n_atom)
        
        data['n_atom'] = n_atom
        cell = []
        
        for i in range(len(file)):
            if 'CELL_PARAMETERS' in file[i]:
                v1 = file[i+1][0:-1].split(' ')
                v1 = list(filter(lambda a: a != '', v1))
                v1 = [float(i) for i in v1]
                
                v2 = file[i+2][:-1].split(' ')
                v2 = list(filter(lambda a: a != '', v2))
                v2 = [float(i) for i in v2]
                
                v3 = file[i+3][:-1].split(' ')
                v3 = list(filter(lambda a: a != '', v3))
                v3 = [float(i) for i in v3]
                cell += [v1, v2, v3]
        
        data['cell_cart'] = cell

        atomic_pos = []

        for i in range(len(file)):
            if 'ATOMIC_POSITIONS' in file[i]:
                
                atoms = [x.rstrip('\n') for x in file[i+1:i+n_atom+1]]
                atoms = [x.split(' ') for x in atoms]
                atomic_pos = [list(filter(lambda a: a != '', x)) for x in atoms]
        
        data['atomic_pos'] = atomic_pos
        
        list_type = []
        list_atom = []
        
        for i in atomic_pos:
            list_atom += [[i[0]]]
            
            if (i[0] in list_type) == False:
                list_type += i[0]
                
        data['list_atom'] = list_atom
        data['list_type'] = list_type
        
        return data
    
    def read_vcrelax_out(self):
        
        if os.path.exists(os.path.join(self.path, self.prefix+'.vcrelax.out')) == True:
            qe_file = os.path.join(self.path, self.prefix+'.vcrelax.out')
            
        else:
            print("Cannot find vcrelax file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE VC-Relax file (*.vcrelax.out);; QE VC-Relax file (*.vc-relax.out) ;;All Files (*.*)")[0]
        
        temp = open(qe_file, 'r')
        file = temp.readlines()
        temp.close()
        
        data = {}
        
        n_atom = 0
        
        for i in file:
            if 'number of atoms/cell ' in i: n_atom = int(i[33:-1])
        data['n_atom'] = n_atom
        
        energy = np.array([])
        for i in file:
            if '!' in i: energy = np.append(energy ,[float(i[32:-3])])
        data['energy'] = energy
        
        force = np.array([])
        for i in file:
            if 'Total force' in i: force = np.append(force ,[float(i[21:35])])
        data['force'] = force
        
        total_stress = np.array([])
        for i in file:
            if ' total   stress' in i: total_stress = np.append(total_stress ,[float(i[71:-1])])
        data['total_stress'] = total_stress
        
        matrix_stress = []
        
        for i in range(len(file)):
            if ' total   stress' in file[i]:
                
                v1 = file[i+1][45:-1].split(' ')
                v1 = list(filter(lambda a: a != '', v1))
                v2 = file[i+2][45:-1].split(' ')
                v2 = list(filter(lambda a: a != '', v2))
                v3 = file[i+3][45:-1].split(' ')
                v3 = list(filter(lambda a: a != '', v3))
                matrix_stress += [[v1, v2, v3]]
        
        data['matrix_stress'] = np.array(matrix_stress)
        
        cell = []
        
        for i in range(len(file)):
            if 'CELL_PARAMETERS' in file[i]:
                v1 = file[i+1][0:-1].split(' ')
                v1 = list(filter(lambda a: a != '', v1))
                v1 = [float(i) for i in v1]
                
                v2 = file[i+2][:-1].split(' ')
                v2 = list(filter(lambda a: a != '', v2))
                v2 = [float(i) for i in v2]
                
                v3 = file[i+3][:-1].split(' ')
                v3 = list(filter(lambda a: a != '', v3))
                v3 = [float(i) for i in v3]
                cell += [[v1, v2, v3]]
        
        data['cell'] = cell
        
        cell_type = []
        
        for i in range(len(file)):
            if 'CELL_PARAMETERS' in file[i]: cell_type = file[i][16:-1]
        
        data['cell_type'] = cell_type
        
        atomic_pos_type = []
        
        for i in range(len(file)):
            if 'ATOMIC_POSITIONS' in file[i]: atomic_pos_type = file[i][16:-1]
        
        data['atomic_pos_type'] = atomic_pos_type
        
        atomic_pos = []
        for i in range(len(file)):
            if 'ATOMIC_POSITIONS' in file[i]:
                atomic_pos += [file[i+1:i+n_atom+1]]
        
        data['atomic_pos'] = atomic_pos
        
        opt_coord = []
        pos_begin = 0
        pos_end = 0
        
        for i in range(len(file)):
            
            if 'Begin final coordinates' in file[i]: pos_begin = i
            if 'End final coordinates' in file[i]: pos_end = i
        
        opt_coord = file[pos_begin+4:pos_end]
        
        data['opt_coord'] = opt_coord
        
        v1 = opt_coord[1].rstrip('\n').split(' ')
        v1 = list(filter(lambda a: a != '', v1))
        v1 = [float(i) for i in v1]
                
        v2 = opt_coord[2].rstrip('\n').split(' ')
        v2 = list(filter(lambda a: a != '', v2))
        v2 = [float(i) for i in v2]
                
        v3 = opt_coord[3].rstrip('\n').split(' ')
        v3 = list(filter(lambda a: a != '', v3))
        v3 = [float(i) for i in v3]
       
        data['opt_cell'] = [v1, v2, v3]
        
        return data

class Bulk_Modulus(QMainWindow):
    
    def __init__(self, prefix):
        super().__init__()
        self.prefix = prefix
        self.path = QFileDialog.getExistingDirectory()

    def get_bulk_data(self):
        
        files_list = os.listdir(self.path)
        file_names = []
        
        for i in files_list:
            if '.scf.out' in i:
                file_names += [i]

        cell_vol_list = []
        cell_energy_list = []
        n_atoms = 1
        for i in file_names:
            tmp = open(os.path.join(self.path, i)).readlines()
            for line in tmp:
                if 'unit-cell volume' in line:
                    cell_vol_list += [float(line[33:].rstrip('(a.u.)^3\n'))]
                elif '!' in line:
                    cell_energy_list += [float(line[33:].rstrip('Ry\n'))]
                elif 'number of atoms/cell' in line:
                    n_atoms = int(line[33:].rstrip('\n'))

        #Convert energy from [Ry] to [eV/atom] and volume from [a.u.^3] to [A^3/atom]
        cell_vol_list = np.array(cell_vol_list)*(0.529177**3)/n_atoms
        cell_energy_list = np.array(cell_energy_list)*13.605698066/n_atoms

        #Save data as a comma separetad value in the current working directory
        save = open(os.path.join(self.path, self.prefix + '.csv'), 'w')
        
        save.write('ev/atom;A^3/atom\n' )
        
        for i in range(len(cell_vol_list)):
            save.write(str(cell_vol_list[i]) + ';' + str(cell_energy_list[i])+'\n')

        save.close()

    def generate_bulk_data(self):

        dir_list = os.listdir(os.getcwd())
        dir_list = [f for f in os.listdir(os.getcwd()) if os.path.isdir(f)==True]
        #dir_list.remove('pyQE.py')
        for structure in dir_list:
            files_list = os.listdir(os.path.join(os.getcwd(), structure))
            
            file_names = []
        
            for i in files_list:
                if '.scf.out' in i:
                    file_names += [i]

            cell_vol_list = []
            cell_energy_list = []
            n_atoms = 1
            for i in file_names:
                tmp = open(os.path.join(os.path.join(os.getcwd(), structure), i)).readlines()
                for line in tmp:
                    if 'unit-cell volume' in line:
                        cell_vol_list += [float(line[33:].rstrip('(a.u.)^3\n'))]
                    elif '!' in line:
                        cell_energy_list += [float(line[33:].rstrip('Ry\n'))]
                    elif 'number of atoms/cell' in line:
                        n_atoms = int(line[33:].rstrip('\n'))

            #Calculate Bulk modulus and convert to GPa - 2*d^2E/dV^2*V0
            fit = np.polyfit(cell_vol_list,cell_energy_list, 2)
            bulk_m = fit[0]*2*1.47e13*1e-9*(-fit[1]/(2*fit[0]))
            print('Bulk Modulus of', structure, bulk_m, 'GPa')

            #Convert volume from [a.u.^3] to [A^3/atom]
            cell_vol_list = np.array(cell_vol_list)*(0.529177**3)/n_atoms
            #Convert energy from [Ry] to [eV/atom]
            cell_energy_list = np.array(cell_energy_list)*13.605698066/n_atoms

            #Save data as a comma separetad value in the current working directory
            save = open(os.path.join(os.getcwd(), structure + '.csv'), 'w')

            for i in range(len(cell_vol_list)):
                save.write(str(cell_vol_list[i]) + ';' + str(cell_energy_list[i])+'\n')

            save.close()

    def plot_bulk(self):

        files_list = [f for f in os.listdir(self.path) if f.endswith('.csv')]
        data = {}

        for name in files_list:
            tmp = np.transpose(np.loadtxt(os.path.join(self.path, name), delimiter=';'))
            data.update({name[:-4]: tmp})

        markers = ['v', 'o', '^', 's', 'P', '*', 'd', 'h', '+', 'x', '8' , '<', '>']

        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)

        # Sort colors by hue, saturation, value and name.
        by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                        for name, color in colors.items())
        sorted_names = [name for hsv, name in by_hsv]
        sorted_names = [i for i in sorted_names if len(i)>3]

        
        for key in data:
            cor = random.choice(sorted_names)
            plt.plot(data[key][0],data[key][1] + 11.84567612*13.605698066, 'o', ms=6, lw=1.5, alpha=1, mfc=cor, color='black', label=key)
            sorted_names.remove(cor)

        plt.ylabel('Energy (eV/atom)', fontsize=18) 
        plt.xlabel(r'Volume (\AA \textsuperscript{3}/atom)', fontsize=18)   
        plt.grid()
        plt.tight_layout()
        plt.legend(loc='lower right', framealpha=1,prop={'size': 10})
        plt.savefig(os.path.join(self.path, 'bulk.pdf'), dpi=150, transparent=False, figsize=(16,16),  format='pdf')
        plt.show()

    def float_string(self, string):

        list_float = []
        tmp = string.split(' ')
        for i in tmp:
            if i != '':
                list_float += [float(i)]
        return list_float

    def get_elastic_matrix(self):

        #open the THERMOPW output  
        tmp = open(os.path.join(self.path, self.prefix + '.elastic.out')).readlines()

        #create a empty tensor
        tensor = [[],[],[],[],[],[]]

        for i in range(len(tmp)):
            if 'Elastic constants C_ij (kbar) ' in tmp[i]:
                tensor[0] = np.array(self.float_string(tmp[i+2][7:].rstrip('\n')))
                tensor[1] = np.array(self.float_string(tmp[i+3][7:].rstrip('\n')))
                tensor[2] = np.array(self.float_string(tmp[i+4][7:].rstrip('\n')))
                tensor[3] = np.array(self.float_string(tmp[i+5][7:].rstrip('\n')))
                tensor[4] = np.array(self.float_string(tmp[i+6][7:].rstrip('\n')))
                tensor[5] = np.array(self.float_string(tmp[i+7][7:].rstrip('\n')))

        #Convert eigenvalues from kbar to GPa        
        elastic_tensor = np.array(tensor)/10

        #Calculate the determinant of the elastic tensor
        det = np.linalg.det(elastic_tensor)
        
        return elastic_tensor, det


class Bands_DOS(QMainWindow):
    
    def __init__(self, prefix):
        super().__init__()
        self.prefix = prefix
        self.path = QFileDialog.getExistingDirectory()
        self.eFermi = self.get_Fermi()
        self.orbitals = [['H',  ['1(s)']], ['He', ['1(s)']], ['Li', ['1(s)']],['Be', ['1(s)']],
                                 ['B',  ['1(s)', '2(p)']],['C',  ['1(s)', '2(p)']],['N',  ['1(s)', '2(p)']] , ['Si', ['1(s)', '2(p)']]]
        
    def Symmetries(self):
        
        '''Imput = (prefix) | Recebe de entrada o arquivo prefix.band.out e retorna um np.array com as coordenadas dos pontos de alta simetria em unidades de 2pi/alat'''
        
        if os.path.exists(os.path.join(self.path, self.prefix+'.band.out')) == True:
            qe_file = os.path.join(self.path, self.prefix+'.band.out')
            
        else:
            print("Cannot find band file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE Band file (*.band.out);; All Files (*.*)")[0]
        temp = open(qe_file,'r')
        x = np.zeros(0)
        for i in temp:
            if "high-symmetry" in i:
                x = np.append(x,float(i.split()[-1]))
        temp.close()
        return np.unique(x)
      
    def get_Fermi(self):
        
        if os.path.exists(os.path.join(self.path, self.prefix+'.nscf.out')) == True:
            
            qe_file = os.path.join(self.path, self.prefix+'.nscf.out')
        else:
            
            if os.path.exists(os.path.join(self.path, self.prefix+'.scf.out')) == True:
                
                qe_file = os.path.join(self.path, self.prefix+'.scf.out')
                
            else: 
                print("Cannot find any scf or nscf out file! Please select the one: ")
                qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "SCF OUT file (*.scf.out);; NSCF OUT file (*.nscf.out);; All Files (*.*)")[0]
            
        temp = open(qe_file, 'r')
        file = temp.readlines()
        temp.close()
        e_fermi = 0
        lumo = 0
    
        for i in file:
            if 'Fermi' in i or 'highest occupied level' in i: 
                raw = i.rstrip('ev\n')
                raw = [k for k in raw.split(' ') if k!='']
                e_fermi = float(raw[-1])
                print('Fermi Energy = ', e_fermi)
                print('Band gap = 0.00 eV')
            if 'highest occupied, lowest unoccupied level' in i:
                raw = i.rstrip('\n')
                raw = [k for k in raw.split(' ') if k!='']
                e_fermi = float(raw[-2])
                lumo = float(raw[-1])
                print('highest occupied level = ', e_fermi, ' eV')
                print('lowest unoccupied level = ', lumo, ' eV')
                print('Band gap = %f eV' %(lumo-e_fermi))
                
        return [e_fermi, lumo]
    
    def subtract_fermi_bands(self):
        
        k,e = self.get_bands()
        fermi = self.eFermi[0]
        e = np.array(e)-fermi
        return k, e    
      
    def subtract_fermi_dos(self):
        
        dos,e,intdos = self.getDOS()      
        fermi = self.eFermi[0]
        e = np.array(e)-fermi
        return dos, e, intdos  
     
    def get_bands(self):
        
        '''Imput = (prefix | Recebe o arquivo prefix.band.gnu e manipula o arquivo para returnar 
        uma lista com os k-points e suas respectivas bandas, subtraindo a energia Fermi de todas 
        as bandas, tornando a energia de Fermi zero.'''
        
        if os.path.exists(os.path.join(self.path, self.prefix+'.band.gnu')) == True:
            qe_file = os.path.join(self.path, self.prefix+'.band.gnu')
            
        else:
            print("Cannot find band file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE Band file (*.band.gnu);; All Files (*.*)")[0]
        
        temp=open(qe_file, 'r')
        bandsfile = temp.readlines()
        temp.close()
        
        bands=[]
        for i in range(len(bandsfile)):
            bandsfile[i] = bandsfile[i].rstrip('\n')
            bandsfile[i] = bandsfile[i].lstrip('    ')
            bandsfile[i] = bandsfile[i].split(' ')     
            if len(bandsfile[i])>2:
                bands.append([float(bandsfile[i][0]),float(bandsfile[i][-1])]) 
                
        #Separar em bandas de condução e valência: usar o np.all() <fermi = valência  
        
        k = []    
        for i in range(len(bands)):
            if (bands[i][0] in k) == False:
                k.append(bands[i][0])
                
        k = np.array(k)
        energy = [[]]*len(k)
    
        for i in range(len(bands)):
            for j in range(len(k)):
                if bands[i][0]==k[j]:
                    energy[j] = energy[j] + [bands[i][1]]

        energy = [*zip(*energy)]
       
        return k, energy
    
    def getDOS(self):
        
        if os.path.exists(os.path.join(self.path, self.prefix+'.dos')) == True:
            qe_file = os.path.join(self.path, self.prefix+'.dos')
            
        else:
            print("Cannot find band file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE DOS file (*.dos);; All Files (*.*)")[0]
            
        temp = open(qe_file, 'r')
        dos = temp.readlines()
        temp.close()
        dosplt = dos[1:]
        dos,e,intdos = [],[],[]
        
        for i in range(len(dosplt)):
            dosplt[i] = dosplt[i].rstrip('\n')
            dosplt[i] = [k for k in dosplt[i].split(' ') if k!='']
            
            e.append([float(dosplt[i][0])])
            dos.append([float(dosplt[i][1])])
            intdos.append([float(dosplt[i][2])])
            
        dos = np.array(dos)
        e = np.array(e)
        intdos=np.array(intdos)
        
        return dos,e,intdos
                
    def get_pdos_file(self, filename):
        
        temp = open(os.path.join(self.path,filename), 'r')
        pdos_file = temp.readlines()[1:]
        temp.close()

        for i in range(len(pdos_file)):
            pdos_file[i]=pdos_file[i].rstrip('\n')
            pdos_file[i]=pdos_file[i].split(' ')
            pdos_file[i] = [float(j) for j in pdos_file[i] if j!='']
            
        pdos_file = np.transpose(pdos_file)
        return pdos_file
        
    def get_all_pdos(self):
            
        list_atoms = Build_Files.read_qe_in(self, self.prefix)['list_atom']
        
        energy=[]
        list_dos_s=[[]]*len(list_atoms)
            
        list_dos_p=[]

        #Adiciona, ao lado de cada elemento, seus respectivos orbitais
        for i in range(len(list_atoms)):
            for j in range(len(self.orbitals)):
                if list_atoms[i][0]==self.orbitals[j][0]: list_atoms[i].append(self.orbitals[j][-1])
        #Itera sobre todos os átomos para pegar seu índice e tipo

        for i in range(len(list_atoms)):
            file = self.prefix+".pdos.pdos_atm#"+str(i+1)+"("+list_atoms[i][0]+")_wfc#"
        #itera sobre os tipos de orbitais (1s ou 2p) e os adiciona em list_dos_s ou list_dos_p
            for j in range(len(list_atoms[i][-1])):
                
                filej = file + list_atoms[i][-1][j]
                
                energy.append(self.get_pdos_file(filej)[0]) #Cria a lista das energias
                if list_atoms[i][-1][j] == '1(s)': #Cria a lista dos orbitais 1s
                    list_dos_s[i] = self.get_pdos_file(file+list_atoms[i][-1][j])[1:]     
                    
                if list_atoms[i][-1][j] == '2(p)': #Cria a lista dos orbiais 2p
                    list_dos_p += [self.get_pdos_file(file+list_atoms[i][-1][j])[1:]]
        
        energy = np.array(energy[0]).astype(np.float)
        
        energy = energy - self.eFermi[0]
        return energy, list_dos_s, list_dos_p
    
    def soma_pdos_s(self):
    
        energy, pdosfile = self.get_all_pdos()[0], self.get_all_pdos()[1]
        dos = np.zeros(len(pdosfile[0][0]))
        for i in range(len(pdosfile)):
            for j in range(len(pdosfile[i][0])):
                dos[j] = dos[j] + float(pdosfile[i][0][j])
              
        return energy, dos
    
    def soma_pdos_p(self):
     
        energy, pdosfile = self.get_all_pdos()[0], self.get_all_pdos()[2]
         
        ldos = np.zeros(len(pdosfile[0][0]))
        xdos = np.zeros(len(pdosfile[0][1]))
        ydos = np.zeros(len(pdosfile[0][2]))
        zdos = np.zeros(len(pdosfile[0][3]))
        
        for i in range(len(pdosfile)):
            for j in range(len(pdosfile[i][0])):
                ldos[j] = ldos[j] + float(pdosfile[i][0][j])
            for j in range(len(pdosfile[i][1])):
                xdos[j] = xdos[j] + float(pdosfile[i][1][j])
            for j in range(len(pdosfile[i][2])):
                ydos[j] = ydos[j] + float(pdosfile[i][2][j])
            for j in range(len(pdosfile[i][3])):
                zdos[j] = zdos[j] + float(pdosfile[i][3][j])
                
        return energy, ldos, xdos, ydos, zdos
    
    def soma_pdos(self):
        
        energy = self.get_all_pdos()[0]
        pdos=(self.soma_pdos_s()[1],self.soma_pdos_p()[1])
        soma=np.zeros(len(pdos[0]))
        for i in range(len(pdos)): 
            for j in range(len(pdos[i])):    
                soma[j] = soma[j] +pdos[i][j]
        #pdos=get_pdos_file(prefix+'.pdos_tot')
        #return pdos[0], pdos[-1]
        return energy, soma

class Bands_Plot(QMainWindow):
    
    def __init__(self, prefix):
        super().__init__()
        self.prefix = prefix
        self.Bands_D = Bands_DOS(self.prefix)
        self.kpath = k_path_label

    def plotband(self, eFermi=True, y = None, figname=False):
     
        if eFermi == True:    
            k, energy = self.Bands_D.subtract_fermi_bands()
        if eFermi == False:
            k, energy = self.Bands_D.get_bands()
        kpoints = self.Bands_D.Symmetries()
            
        ymin, ymax = 0,0
        if y != None: ymin, ymax = y
        if y == None:
            for i in range(len(energy)):
                if max(energy[i]) > ymax: ymax = max(energy[i]) 
                if min(energy[i]) < ymin: ymin = min(energy[i])
                    
        plt.figure(dpi=170)
        for i in range(len(energy)):
            plt.plot(k,energy[i], lw=0.55)
        plt.plot([kpoints[1:-1], kpoints[1:-1]],[ymin, ymax],'--',lw=0.55,color='red',alpha=0.75)    
        plt.axis([min(k), max(k), ymin, ymax])
        if eFermi==True:
            plt.plot([min(k), max(k)], [0, 0], '-',lw=0.55,color='red',alpha=0.75)
        plt.ylabel('Energy (eV)')
        plt.xticks(kpoints, self.kpath)
        if figname!= False: plt.savefig(figname, dpi=150, transparent=True)
        plt.show()
        
       
    def plotDOS(self):
            
        x,y,z = self.Bands_D.getDOS()
        plt.plot(y,x)
        plt.fill(x, alpha=0.7) 
        plt.axis([min(y),max(y),0,max(x)*0.3])
        plt.show()
    
    def pltbandasDOS(self, eFermi = True, y = [-5,5], figname=True):
      
        if eFermi == True:    
            k, energy = self.Bands_D.subtract_fermi_bands()
        if eFermi == False:
            k, energy = self.Bands_D.get_bands()
        kpoints = self.Bands_D.Symmetries()
            
        ymin, ymax = -5,5
        if y != None: ymin, ymax = y
        if y == None:
            for i in range(len(energy)):
                if max(energy[i]) > ymax: ymax = max(energy[i]) 
                if min(energy[i]) < ymin: ymin = min(energy[i])
                 
        #plt.figure(figsize = (16,12),dpi = 150)
        plt.figure()
        ax1 = plt.subplot2grid((1,8), (0,0), colspan = 5)
        ax2 = plt.subplot2grid((1,8), (0,5), colspan = 3)           
        for i in range(len(energy)):
            #ax1.plot(k,energy[i],'.',  ms=1 , lw=0.6, color='black')
	        ax1.plot(k,energy[i], color='black')
        ax1.plot([kpoints[1:-1], kpoints[1:-1]],[ymin, ymax],'--',lw=0.55,color='red')    
        ax1.axis([min(k), max(k), ymin, ymax])
        ax1.plot([min(k), max(k)], [0, 0], '-',lw=0.55,color='red',alpha=0.75)
        ax1.set_ylabel('Energy (E-$\epsilon_f$) [eV]')
        ax1.set_xlabel('K-path')
        ax1.set_xticks(kpoints)
        ax1.set_xticklabels(self.kpath)
        ax1.set_title("Band diagram")
    
        
        dose, energy, intdos = self.Bands_D.subtract_fermi_dos()
            
        ax2.plot(dose,energy, '-',color='darkgreen',alpha=1, lw=1.)
        ax2.set_ylim(ax1.get_ylim())
        ax2.set_xlim(min(dose)[0], max(dose)[0]*1.1)
        ax2.plot([min(dose)[0],max(dose)[0]*1.1], [0, 0], '-',lw=0.55,color='red',alpha=0.75)
        ax2.fill(dose,energy, color='green',alpha=0.7)        
        ax2.set_yticks(())
        ax2.set_title("DOS distribution")
        #ax2.legend(loc='upper right', framealpha=0,prop={'size': 7.5})
        ax2.set_xlabel('DOS($\epsilon$) [arb. units]')
        plt.tight_layout()
    
        if figname!= False: 
            plt.savefig(os.path.join(self.Bands_D.path, self.prefix+'.bands.pdf'), dpi=150, transparent=False, figsize=(16,12))
        plt.show() 
        
    def pltbandaspDOS(self, eFermi = True, y = [-5,5], figname=True):
      
        if eFermi == True:    
            k, energy = self.Bands_D.subtract_fermi_bands()
        if eFermi == False:
            k, energy = self.Bands_D.get_bands()
        kpoints = self.Bands_D.Symmetries()
            
        ymin, ymax = -5,5
        if y != None: ymin, ymax = y
        if y == None:
            for i in range(len(energy)):
                if max(energy[i]) > ymax: ymax = max(energy[i]) 
                if min(energy[i]) < ymin: ymin = min(energy[i])
                 
        #plt.figure(figsize = (16,12),dpi = 150)
        plt.figure()
        ax1 = plt.subplot2grid((1,8), (0,0), colspan = 5)
        ax2 = plt.subplot2grid((1,8), (0,5), colspan = 3)             
        for i in range(len(energy)):
            ax1.plot(k,energy[i],'.',  ms=1 , lw=0.6, color='black')
            #ax1.plot(k,energy[i],lw=0.6, color='black')
        ax1.plot([kpoints[1:-1], kpoints[1:-1]],[ymin, ymax],'--',lw=0.55,color='red')    
        ax1.axis([min(k), max(k), ymin, ymax])
        ax1.plot([min(k), max(k)], [0, 0], '-',lw=0.55,color='red',alpha=0.75)
        ax1.set_ylabel('Energy (E-$\epsilon_f$) [eV]')
        ax1.set_xlabel('K-path')
        ax1.set_xticks(kpoints)
        ax1.set_xticklabels(self.kpath)
        ax1.set_title("Band diagram")
    
        
        dose, energy, intdos = self.Bands_D.subtract_fermi_dos()
        energy, pdos_s = self.Bands_D.soma_pdos_s()
        pdos_p = self.Bands_D.soma_pdos_p()
        pdos_total= self.Bands_D.soma_pdos()[1]

        ax2.plot(pdos_total, energy, label="Total pDOS", color='gray')
        ax2.plot(pdos_p[1], energy, label="pDOS Orbital P", color='r')   
        ax2.plot(pdos_s, energy, label="pDOS Orbital S", color='b')

        plt.fill_between(pdos_total, 0, energy, alpha=.4, color='gray', linewidth=0.0)
        plt.fill_between(pdos_p[1], 0, energy, alpha=.4, color='r', linewidth=0.0)
        plt.fill_between(pdos_s, 0, energy, alpha=.4, color='b', linewidth=0.0)  

        ax2.set_ylim(ax1.get_ylim())
        ax2.set_xlim(min(dose)[0], max(dose)[0]*0.45)
        ax2.plot([min(dose)[0],max(dose)[0]*.8], [0, 0], '-',lw=0.55,color='red',alpha=0.75)
        #ax2.fill(dose,energy, color='green',alpha=0.7)        
        ax2.set_yticks(())
        ax2.set_title("pDOS dist.")
        ax2.legend(loc='upper right', framealpha=0,prop={'size': 7.5})
        ax2.set_xlabel('pDOS($\epsilon$) [a. u.]')
        plt.tight_layout()
    
        if figname!= False: 
            plt.savefig(os.path.join(self.Bands_D.path, self.prefix+'.bands.pdf'), dpi=150, transparent=False, figsize=(16,12))
        plt.show() 
        
    def pltbandas_all_pDOS(self, eFermi = True, y = [-5,5], figname=True):
      
        if eFermi == True:    
            k, energy = self.Bands_D.subtract_fermi_bands()
        if eFermi == False:
            k, energy = self.Bands_D.get_bands()
        kpoints = self.Bands_D.Symmetries()
            
        ymin, ymax = -5,5
        if y != None: ymin, ymax = y
        if y == None:
            for i in range(len(energy)):
                if max(energy[i]) > ymax: ymax = max(energy[i]) 
                if min(energy[i]) < ymin: ymin = min(energy[i])
                 
        #plt.figure(figsize = (16,12),dpi = 150)
        
        plt.figure()
        ax1 = plt.subplot2grid((1,8), (0,0), colspan = 5)
        ax2 = plt.subplot2grid((1,8), (0,5), colspan = 3)               
        for i in range(len(energy)):
            #ax1.plot(k,energy[i],'.',  ms=1 , lw=0.6, color='black')
            ax1.plot(k,energy[i],lw=1.6, color='black')
        ax1.plot([kpoints[1:-1], kpoints[1:-1]],[ymin, ymax],'--',lw=0.55,color='red')    
        ax1.axis([min(k), max(k), ymin, ymax])
        ax1.plot([min(k), max(k)], [0, 0], '-',lw=0.55,color='red',alpha=0.75)
        ax1.set_ylabel('Energy (E-$\epsilon_f$) [eV]')
        ax1.set_xlabel('K-path')
        ax1.set_xticks(kpoints)
        ax1.set_xticklabels(self.kpath)
        ax1.set_title("Band diagram")
          
        dose, energy, intdos = self.Bands_D.subtract_fermi_dos()
        energy, pdos_s = self.Bands_D.soma_pdos_s()
        pdos_p = self.Bands_D.soma_pdos_p()
        pdos_total= self.Bands_D.soma_pdos()[1]
        
        ap = 0.4
        
        ax2.plot(pdos_s, energy, label="orbital s", color='b')
        ax2.plot(pdos_p[2], energy, label=r"orbital p$_x$", color='red')
        #ax2.plot(pdos_p[3], energy, label=r"orbital p$_y$ and p$_z$", color='green')
        ax2.plot(pdos_p[3], energy, label=r"orbital p$_y$", color='green')
        ax2.plot(pdos_p[4], energy, label=r"orbital p$_z$", color='yellow')
        #ax2.plot(pdos_total, 0, energy, label='Total DOS', color='gray')    

        ax2.fill_between(pdos_s, 0, energy, alpha=ap, color='blue', linewidth=0.0)
        ax2.fill_between(pdos_p[2], 0, energy, alpha=ap, color='red', linewidth=0.0)
        ax2.fill_between(pdos_p[3], 0, energy, alpha=ap, color='green', linewidth=0.0)
        ax2.fill_between(pdos_p[4], 0, energy, alpha=ap, color='yellow', linewidth=0.0)
        ax2.fill_between(pdos_total, 0, energy, alpha=ap, color='gray', linewidth=0.0)


        ax2.set_ylim(ax1.get_ylim())
        ax2.set_xlim(0, 3)
        ax2.plot([min(dose)[0],max(dose)[0]*.8], [0, 0], '-',lw=0.55,color='red',alpha=0.75)
        
        #ax2.fill(dose,energy, color='green',alpha=0.7)        
        ax2.set_yticks(())
        ax2.set_title("pDOS dist.")
        ax2.legend(loc='upper right', framealpha=1,prop={'size': 10})
        ax2.set_xlabel('pDOS($\epsilon$) [states/eV]')
        plt.tight_layout()
    

        plt.savefig(os.path.join(self.Bands_D.path, self.prefix+'.bands.pdf'), dpi=150, transparent=False, figsize=(16,12))
        plt.savefig(os.path.join(self.Bands_D.path, self.prefix+'.bands.svg'), dpi=150, transparent=True, figsize=(16,12))
        plt.savefig(os.path.join(self.Bands_D.path, self.prefix+'.bands.png'), dpi=150, transparent=True, figsize=(16,12))
        plt.show() 

    def plot_pdos(self, x=None, y=None, sump=True, figname=True):
        
        energy, pdos_s = self.Bands_D.soma_pdos_s()
        pdos_p = self.Bands_D.soma_pdos_p()
        pdos_total= self.Bands_D.soma_pdos()[1]
        
        plt.plot(energy, pdos_total, label="Total pDOS", color='gray')
        if sump==True: plt.plot(energy, pdos_p[1], label="pDOS Orbital P", color='r')   
        if sump==True: plt.plot(energy, pdos_s, label="pDOS Orbital S", color='b')

        ap = 0.4
        
        if sump!=True:
            
            plt.plot(energy, pdos_p[2], label="pDOS Orbital Px", color='r')
            plt.plot(energy, pdos_p[3], label="pDOS Orbital Py", color='green')
            plt.plot(energy, pdos_p[4], label="pDOS Orbital Pz", color='yellow')
            plt.plot(energy, pdos_s, label="pDOS Orbital S", color='b')
            
            plt.fill_between(energy, 0, pdos_p[2], alpha=ap, color='r', linewidth=0.0)
            plt.fill_between(energy, 0, pdos_p[3], alpha=ap, color='green', linewidth=0.0)
            plt.fill_between(energy, 0, pdos_p[4], alpha=ap, color='yellow', linewidth=0.0)
            plt.fill_between(energy, 0, pdos_total, alpha=ap, color='gray', linewidth=0.0)
            plt.fill_between(energy, 0, pdos_s, alpha=ap, color='gray', linewidth=0.0)
        
        if sump==True:
            plt.fill_between(energy, 0, pdos_total, alpha=ap, color='gray', linewidth=0.0)
            plt.fill_between(energy, 0, pdos_p[1], alpha=ap, color='r', linewidth=0.0)
            plt.fill_between(energy, 0, pdos_s, alpha=ap, color='b', linewidth=0.0)        
        
        plt.legend(loc='upper right', framealpha=0,prop={'size': 12.5})
        plt.tight_layout()
        
        if x!=None: plt.axis([x[0], x[1], 0, float(max(pdos_s))])
        if y!=None: plt.axis([x[0], x[1], 0, y])
        if x==None and y==None: plt.axis([float(min(energy)), float(max(energy)), 0, float(max(pdos_total))])
        
        plt.xlabel('Energy (E-$\epsilon_f$) [eV]')
        plt.ylabel('DOS($\epsilon$) [arb. units]')  

        if figname!= False: 
            plt.savefig(os.path.join(self.Bands_D.path, self.prefix + '.pdos.pdf'), dpi=150, transparent=False, figsize=(16,12))
            plt.savefig(os.path.join(self.Bands_D.path, self.prefix + '.pdos.svg'), dpi=150, transparent=False, figsize=(16,12))
        
        plt.close()
        
class Plot_Phonon(QMainWindow):
    
    def __init__(self, prefix):
        super().__init__()

        self.prefix = prefix
        self.kpath = k_path_label
        self.path = QFileDialog.getExistingDirectory()
        
    def get_sym(self):
        
        if os.path.exists(os.path.join(self.path, self.prefix+'.plotband.out')) == True:
            print('File ', self.prefix+'.plotband.out', ' found!')
            qe_file = os.path.join(self.path, self.prefix+'.plotband.out')
            
        else:
            print("Cannot find ", self.prefix+'.plotband.out', " file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE Band file (*.plotband.out);; All Files (*.*)")[0]
            
        temp = open(qe_file,'r').readlines()
        sym = np.zeros(0)
        for i in range(len(temp)):
            if "coordinate" in temp[i]:
                temp[i] = temp[i].rstrip('\n')
                temp[i] = temp[i].split('   ')
                sym = np.append(sym, float(temp[i][-1]))
        return np.unique(sym)
            
    def get_freq(self):
        
        if os.path.exists(os.path.join(self.path, 'freq.plot')) == True:
            print('File ', self.prefix+'.freq.plot', ' found!')
            qe_file = os.path.join(self.path, 'freq.plot')
            
        else:
            print("Cannot find", self.prefix+'.freq.plot' ," file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE Plot file (*.plot);; All Files (*.*)")[0]
            
        temp=open(qe_file, 'r')
        bandsfile = temp.readlines()
        temp.close()
        
        bands=[]
        for i in range(len(bandsfile)):
            bandsfile[i] = bandsfile[i].rstrip('\n')
            bandsfile[i] = bandsfile[i].lstrip('  ')
            bandsfile[i] = bandsfile[i].split(' ')     
            
            if bandsfile[i] != ['']:
                bands.append([float(bandsfile[i][0]),float(bandsfile[i][-1])]) 
    
    
        k = []    
        for i in range(len(bands)):
            if (bands[i][0] in k) == False:
                k.append(bands[i][0])
        k = np.array(k)
        
        energy = [[]]*len(k)
    
        for i in range(len(bands)):
            for j in range(len(k)):
                if bands[i][0]==k[j]:
                    energy[j] = energy[j] + [bands[i][1]]
    
        #e = np.transpose(energy)
        e = [*zip(*energy)]
        
        return k, e
    
    def plot_freq(self):
        
        k, e = self.get_freq()
        HS = self.get_sym()
        ymin, ymax = -100, 2000  
            
        for i in range(len(e)):
            plt.plot(k,e[i], lw=1, color='black')
            
        for i in HS[1:-1]:
            plt.plot([i,i], [ymin,ymax], lw=.6, color='red')
            
        plt.plot([min(k), max(k)], [0,0],'--', lw=.8, color='red')
        plt.ylim([ymin,ymax])
        plt.xlim(min(k), max(k))
        plt.xticks(HS, self.kpath, size=12)
        plt.yticks(size=12)
        
        plt.ylabel(r'Frequency (cm$^{-1}$)', size=14)
        plt.xlabel('K-path', size=14)
        plt.show()
    
    def get_DOS(self):
        
        if os.path.exists(os.path.join(self.path, self.prefix+'.phdos')) == True:
            print('File ', self.prefix+'.phdos', ' found!')
            qe_file = os.path.join(self.path, self.prefix+'.phdos')
            
        else:
            print("Cannot find ",self.prefix+'.phdos'," file! Please select the file: ")
            qe_file = QFileDialog.getOpenFileName(self, 'Open Quantum ESPRESSO file', self.path, "QE PHDos file (*.phdos);; All Files (*.*)")[0]
            
        dos = []
        temp = open(qe_file, 'r').readlines()
        temp = temp[1:]
        for i in range(len(temp)):
            temp[i] = temp[i].rstrip('\n')
            temp[i] = temp[i].split('  ')
            dos += [[float(i) for i in temp[i] if i!='']]
            
        return np.transpose(dos)[0:2]
        
    
    def pltbandasDOS(self, figname=True):
      
        k, energy = self.get_freq()
        kpoints = self.get_sym()
        
        ymin, ymax = -100, 1500       
        plt.figure(figsize = (9,5),dpi = 50)
        #plt.figure()
        ax1 = plt.subplot2grid((1,8), (0,0), colspan = 5)
        ax2 = plt.subplot2grid((1,8), (0,5), colspan = 3)   
          
        for i in range(len(energy)):
            #ax1.plot(k,energy[i],'.', ms=1 , lw=0.7,color='black')
            ax1.plot(k,energy[i],lw=1.6, color='black')
        ax1.plot([kpoints[1:-1], kpoints[1:-1]],[ymin, ymax],lw=0.6,color='red')    
        ax1.axis([min(k), max(k), ymin, ymax])
        ax1.plot([min(k), max(k)], [0,0],'--', lw=.6, color='red')
        
        if os.path.exists(os.path.join(self.path, 'Exp.csv')) == True:
            print('File ', 'Exp.csv', ' found!')
            exp =np.transpose( np.genfromtxt(str(self.path) + '/Exp.csv', delimiter=';'))
            ax1.plot(exp[0], exp[2], '.',  ms=8, lw=2, mfc='red', color='black')
        
        
        ax1.set_ylabel(r'Frequency [cm\textsuperscript{-1}]')
        
        ax1.set_xlabel('K-path')
        ax1.set_xticks(kpoints)
        ax1.set_xticklabels(self.kpath)
        ax1.set_title("Phonon Dispersion")
    
        freq, dos = self.get_DOS()
            
        ax2.plot(dos,freq, '-',color='b',alpha=1, lw=1.,label="Phonon DOS")
        ax2.set_ylim(ax1.get_ylim())
        ax2.set_xlim(min(dos), max(dos)*1.1)
        ax2.plot([min(dos), max(dos)*1.1], [0,0],'--', lw=.6, color='red')
        ax2.fill(dos , freq, alpha=0.4, color='darkblue')        
        ax2.set_yticks(())
        ax2.set_title("Vibrational DOS")
        #ax2.legend(loc='upper right', framealpha=0,prop={'size': 7.5})
        ax2.set_xlabel(r'VDOS [states/cm\textsuperscript{-1}]')
        plt.tight_layout()
        if figname!= False: 
            plt.savefig(os.path.join(self.path, self.prefix+'_phdisp.pdf'), dpi=150, transparent=False, figsize=(16,12))
            plt.savefig(os.path.join(self.path, self.prefix+'_phdisp.svg'), dpi=150, transparent=True, figsize=(16,12))
            plt.savefig(os.path.join(self.path, self.prefix+'_phdisp.png'), dpi=150, transparent=True, figsize=(16,12))
        
        plt.close()
        


print('\n')
print('----------------------------------------------------------------------------------------------------------------------- \n')
print('                             pyQE - A front-end suite for Quantum ESPRESSO data analysis\n')
print('                                         Code created by Felipe Lopes \n')
print('----------------------------------------------------------------------------------------------------------------------- \n')

a = QApplication([])
#w.show()
#w.hide()

print('What job do you want to do?')
print('1 - Plot band structure and DOS')
print('11 - Plot band structure and summed pDOS')
print('12 - Plot band structure and all pDOS')
print('2 - Plot summed pDOS')
print('21 - Plot all pDOS')
print('3 - Plot phonon dispersion')
print('4 - Plot Bulk modulus information')
print('5 - Elastic')
    
temp = input()
    
if temp == '1':
    print('Prefix:')
    prefix = input()
    band = Bands_Plot(prefix)
    band.pltbandasDOS()

if temp == '11':
    print('Prefix:')
    prefix = input()
    band = Bands_Plot(prefix)
    band.pltbandaspDOS()

if temp == '12':
    print('Prefix:')
    prefix = input()
    band = Bands_Plot(prefix)
    band.pltbandas_all_pDOS()

if temp == '2':
    print('Prefix:')
    prefix = input()
    band = Bands_Plot(prefix)
    band.plot_pdos(x=[-5,5],y=6)

if temp == '21':
    print('Prefix:')
    prefix = input()
    band = Bands_Plot(prefix)
    band.plot_pdos(sump = False, x=[-5,5],y=6)

elif temp == '3':

    print('Prefix:')
    prefix = input()
    ph = Plot_Phonon(prefix)
    print('Plotting phonon dispersion...')
    ph.pltbandasDOS()

elif temp == '4':
    print('Prefix:')
    prefix = input()
    bulk = Bulk_Modulus(prefix)
    bulk.get_bulk_data()

elif temp == '5':
    print('Prefix:')
    prefix = input()
    bulk = Bulk_Modulus(prefix)
    bulk.get_elastic_matrix()

        
