# In v8 we have added H2_dissocJ and dust_G_parameter 
# to the cross sections table. 
import numpy as np
import re
from scipy import integrate
from scipy import interpolate
import h5py
import os 
import sys 
from mpi4py import MPI

try: 
    sys.path.insert(1, os.path.join(sys.path[0], '..')) 
    import read_chimes.read_chimes as rc
except ImportError: 
    sys.path.insert(1, os.path.join(sys.path[0], '..', 'read_chimes')) 
    import read_chimes as rc

min_float_value = 1.0e-100

def read_in_V95_cross_sections(filename):
    file_id = open(filename, "r")
    C = np.zeros((6, 6, 3))
    N = np.zeros((7, 6, 3))
    O = np.zeros((8, 6, 3))
    Ne = np.zeros((10, 6, 3))
    Mg = np.zeros((12, 6, 4))
    Si = np.zeros((14, 6, 5))
    S = np.zeros((16, 6, 5))
    Ca = np.zeros((20,6, 6))
    Fe = np.zeros((26, 6, 7))

    shell_quantum_nos = np.zeros((9, 26, 7, 2))     # dimensions: element, ion, shell, quantum no.
    C_shell = np.zeros(6, dtype = np.int)
    N_shell = np.zeros(7, dtype = np.int)
    O_shell = np.zeros(8, dtype = np.int)
    Ne_shell = np.zeros(10, dtype = np.int)
    Mg_shell = np.zeros(12, dtype = np.int)
    Si_shell = np.zeros(14, dtype = np.int)
    S_shell = np.zeros(16, dtype = np.int)
    Ca_shell = np.zeros(20, dtype = np.int)
    Fe_shell = np.zeros(26, dtype = np.int)
    
    for line in file_id:
        if re.search("#", line) != None:
            continue
        else:
            values = line.split()
            values[-1].strip('\n')
        
            if int(values[0]) == 6:
                for i in range(6):
                    C[6 - int(values[1])][i][C_shell[6 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[0][6 - int(values[1])][C_shell[6 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[0][6 - int(values[1])][C_shell[6 - int(values[1])]][1] = int(values[3])
                C_shell[6 - int(values[1])] += 1
        
            elif int(values[0]) == 7:
                for i in range(6):
                    N[7 - int(values[1])][i][N_shell[7 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[1][7 - int(values[1])][N_shell[7 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[1][7 - int(values[1])][N_shell[7 - int(values[1])]][1] = int(values[3])
                N_shell[7 - int(values[1])] += 1
        
            elif int(values[0]) == 8:
                for i in range(6):
                    O[8 - int(values[1])][i][O_shell[8 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[2][8 - int(values[1])][O_shell[8 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[2][8 - int(values[1])][O_shell[8 - int(values[1])]][1] = int(values[3])
                O_shell[8 - int(values[1])] += 1
        
            elif int(values[0]) == 10:
                for i in range(6):
                    Ne[10 - int(values[1])][i][Ne_shell[10 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[3][10 - int(values[1])][Ne_shell[10 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[3][10 - int(values[1])][Ne_shell[10 - int(values[1])]][1] = int(values[3])
                Ne_shell[10 - int(values[1])] += 1
        
            elif int(values[0]) == 12:
                for i in range(6):
                    Mg[12 - int(values[1])][i][Mg_shell[12 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[4][12 - int(values[1])][Mg_shell[12 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[4][12 - int(values[1])][Mg_shell[12 - int(values[1])]][1] = int(values[3])
                Mg_shell[12 - int(values[1])] += 1
        
            elif int(values[0]) == 14:
                for i in range(6):
                    Si[14 - int(values[1])][i][Si_shell[14 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[5][14 - int(values[1])][Si_shell[14 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[5][14 - int(values[1])][Si_shell[14 - int(values[1])]][1] = int(values[3])
                Si_shell[14 - int(values[1])] += 1
        
            elif int(values[0]) == 16:
                for i in range(6):
                    S[16 - int(values[1])][i][S_shell[16 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[6][16 - int(values[1])][S_shell[16 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[6][16 - int(values[1])][S_shell[16 - int(values[1])]][1] = int(values[3])
                S_shell[16 - int(values[1])] += 1
        
            elif int(values[0]) == 20:
                for i in range(6):
                    Ca[20 - int(values[1])][i][Ca_shell[20 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[7][20 - int(values[1])][Ca_shell[20 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[7][20 - int(values[1])][Ca_shell[20 - int(values[1])]][1] = int(values[3])
                Ca_shell[20 - int(values[1])] += 1
        
            elif int(values[0]) == 26:
                for i in range(6):
                    Fe[26 - int(values[1])][i][Fe_shell[26 - int(values[1])]] = float(values[i+4])
                shell_quantum_nos[8][26 - int(values[1])][Fe_shell[26 - int(values[1])]][0] = int(values[2])
                shell_quantum_nos[8][26 - int(values[1])][Fe_shell[26 - int(values[1])]][1] = int(values[3])
                Fe_shell[26 - int(values[1])] += 1
                
    file_id.close()
    
    return np.array([C, N, O, Ne, Mg, Si, S, Ca, Fe, shell_quantum_nos]), np.array([C_shell, N_shell, O_shell, Ne_shell, Mg_shell, Si_shell, S_shell, Ca_shell, Fe_shell])

def read_in_KM93_tables(filename):
    N_Aug_C = 2
    N_Aug_N = 2
    N_Aug_O = 2
    N_Aug_Ne = 2
    N_Aug_Mg = 4
    N_Aug_Si = 5
    N_Aug_S = 6
    N_Aug_Ca = 8
    N_Aug_Fe = 10

    file_id = open(filename, "r")
    C = np.zeros((6, 3, N_Aug_C))
    N = np.zeros((7, 3, N_Aug_N))
    O = np.zeros((8, 3, N_Aug_O))
    Ne = np.zeros((10, 3, N_Aug_Ne))
    Mg = np.zeros((12, 4, N_Aug_Mg))
    Si = np.zeros((14, 5, N_Aug_Si))
    S = np.zeros((16, 5, N_Aug_S))
    Ca = np.zeros((20, 7, N_Aug_Ca))
    Fe = np.zeros((26, 7, N_Aug_Fe))

    for line in file_id:
        if re.search('#', line) != None:
            continue
        else:
            values = line.split()
            values[-1].strip('\n')
            if int(values[0]) == 6:
                for i in range(N_Aug_C):
                    C[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 7:
                for i in range(N_Aug_N):
                    N[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 8:
                for i in range(N_Aug_O):
                    O[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 10:
                for i in range(N_Aug_Ne):
                    Ne[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 12:
                for i in range(N_Aug_Mg):
                    Mg[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 14:
                for i in range(N_Aug_Si):
                    Si[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 16:
                for i in range(N_Aug_S):
                    S[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 20:
                for i in range(N_Aug_Ca):
                    Ca[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
            if int(values[0]) == 26:
                for i in range(N_Aug_Fe):
                    Fe[int(values[1]) - 1][int(values[2]) - 1][i] = float(values[i+4])
    
    return np.array([C, N, O, Ne, Mg, Si, S, Ca, Fe])


def read_in_V96_cross_sections(filename):
        
    file_id = open(filename, "r")
    H = np.zeros((1, 9))
    He = np.zeros((2, 9))
    C = np.zeros((6, 9))
    N = np.zeros((7, 9))
    O = np.zeros((8, 9))
    Ne = np.zeros((10, 9))
    Mg = np.zeros((12,9))
    Si = np.zeros((14, 9))
    S = np.zeros((16, 9))
    Ca = np.zeros((20,9))
    Fe = np.zeros((26, 9))
    
    for line in file_id:
        if re.search("#", line) != None:
            continue
        else:
            values = line.split()
            values[-1].strip('\n')
        
            if int(values[0]) == 1:
                for i in range(9):
                    H[1 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 2:
                for i in range(9):
                    He[2 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 6:
                for i in range(9):
                    C[6 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 7:
                for i in range(9):
                    N[7 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 8:
                for i in range(9):
                    O[8 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 10:
                for i in range(9):
                    Ne[10 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 12:
                for i in range(9):
                    Mg[12 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 14:
                for i in range(9):
                    Si[14 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 16:
                for i in range(9):
                    S[16 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 20:
                for i in range(9):
                    Ca[20 - int(values[1])][i] = float(values[i+2])
        
            if int(values[0]) == 26:
                for i in range(9):
                    Fe[26 - int(values[1])][i] = float(values[i+2])
                
    file_id.close()
    
    return np.array([H, He, C, N, O, Ne, Mg, Si, S, Ca, Fe])
    
def read_in_spectrum_file(filename):
    file_id = open(filename, "r")
    spectrum_E_log = []     #In Ryd
    spectrum_J_log = []     #In erg/sr/cm^2/s/Hz
  
    for line in file_id:
        if re.search("#", line) != None:
            continue
        else: 
            input_data = line.split()
            input_data[-1].strip('\n')
            spectrum_E_log.append(float(input_data[0]))
            spectrum_J_log.append(float(input_data[1]))
    file_id.close()
  
    return np.array(spectrum_E_log), np.array(spectrum_J_log)

def sigma(E, pars, l):
    # Note: we ignore l - it is taken to be zero.
    E_0 = pars[0]
    sigma_0 = pars[1]
    y_a = pars[2]
    P = pars[3]
    y_w = pars[4]
    y_0 = pars[5]
    y_1 = pars[6]
    
    x = (E / E_0) - y_0
    y = np.sqrt((x**2.0) + (y_1**2.0))
    # Returns the cross section in cm^2
    return sigma_0 * 1.0e-18 * (((x - 1.0) ** 2.0) + (y_w ** 2.0)) * (y ** ((0.5 * P) - 5.5)) * ((1.0 + np.sqrt(y / y_a)) ** -P)
      
# The following two variants ensure that zero
# is returned below the threshold energy. 
def sigma_Hydrogen(E, pars, l):
    if E < 13.6:
        return 0.0
    else:
        return sigma(E, pars, l)
      
def sigma_Helium(E, pars, l):
    if E < 24.59:
        return 0.0
    else:
        return sigma(E, pars, l)

def sigma_V95(E, pars, l):
    E_0 = pars[0]
    sigma_0 = pars[1]
    y_a = pars[2]
    P = pars[3]
    y_w = pars[4]

    x = E / E_0
    # Returns the cross section in cm^2
    return sigma_0 * 1.0e-18 * (((x - 1.0) ** 2.0) + (y_w ** 2.0)) * (x ** ((0.5 * P) - 5.5 - l)) * ((1.0 + np.sqrt(x / y_a)) ** -P)

def sigma_H2(E):
    # This is the cross section (in cm^2) for the
    # photoionisation of H2. Note that 
    # the cross section in GJ07 used an
    # old reference, and did not go above 
    # 30 eV. We therefore replace it with
    # the refernces given below.
    x = E / 15.4
    if E < 15.4:
        return 0.0
    elif E < 18.0:
        return 1.0e-16 * (-37.895 + (99.723 * x) - (87.227 * (x ** 2)) + (25.4 * (x ** 3)))  # Yan et al. (1998); eqn 17
    elif E < 30.0:
        return 2.0e-17 * ((0.071 * (x ** -0.252)) - (0.673 * (x ** -1.252)) + (1.977 * (x ** -2.252)) - (0.692 * (x ** -3.252))) # Yan et al. (1998); eqn 18
    elif E < 85.0:
        return 1.0e-18 * (0.664 - (11.768 * (x ** -1.0)) + (78.118 * (x ** -2.0)) - (231.339 * (x ** -3.0)) + (368.053 * (x ** -4.0)) - (189.953 * (x ** -5.0))) # Wilms et al. (2000); eqn 6
    else:
        return 4.557e-23 * (1.0 - (2.003 * (x ** -0.5)) - (4.806 * (x ** -1.0)) + (50.577 * (x ** -1.5)) - (171.044 * (x ** -2.0)) + (231.608 * (x ** -2.5)) - (81.885 * (x ** -3.0))) * ((E / 1.0e3) ** -3.5)# Yan et al. (1998); eqn 19
    
def sigma_Hminus(E):
    if E > 0.755:
        return 2.11e-16 * ((E - 0.755) ** 1.5) * (E ** -3.0)
    else:
        return 0.0
      
def find_photoion_reaction(species_idx, reactants):
    ind = (reactants == species_idx)
    if sum(ind) > 1:
        raise Exception("species %d found multiple times in the reactant list. Aborting")
    else:
        if sum(ind) == 0:
            return -1
        else: 
            return np.arange(len(reactants))[ind][0]
        
def find_auger_reaction(species_idx, N_elec, reactants, products):
    ind = ((reactants == species_idx) & (products[:, 0] == species_idx + N_elec)) 
    if sum(ind) > 1:
        raise Exception("species %d found multiple times in the reactant list. Aborting")
    else:
        if sum(ind) == 0:
            return -1
        else: 
            return np.arange(len(reactants))[ind][0]

def read_parameter_file(infile):
    parameters = {
        "chimes_main_data_path" : "../../chimes-data/chimes_main_data.hdf5", 
        "cross_sections_data_path" : "data",
        "spectrum_file" : "spectrum.dat", 
        "output_file" : "cross_sections.hdf5"
    }

    fd = open(infile)

    for line in fd:
        if len(line.strip()) < 1:
            continue
        elif line.strip()[0] == "#":
            continue
        else:
            values = line.split()
            if values[0] in parameters:
                parameters[values[0]] = values[1]
            else:
                raise KeyError("Parameter %s not recognised. Aborting" % (values[0], ))

    return parameters 
        
def main():
    parameter_file = sys.argv[1]
    pars = read_parameter_file(parameter_file) 
    
    # Set up MPI variables
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    N_task = comm.Get_size() 
    
    N_Aug_C = 2
    N_Aug_N = 2
    N_Aug_O = 2
    N_Aug_Ne = 2
    N_Aug_Mg = 4
    N_Aug_Si = 5
    N_Aug_S = 6
    N_Aug_Ca = 8
    N_Aug_Fe = 10
    
    log_N = np.arange(15.0, 24.1, 0.1)
    N_x = 10.0 ** log_N         # Use same bins for NHI and NHeI.

    V96_data_file = "%s/cross_section_fits_verner96.dat" % (pars["cross_sections_data_path"], )
    x_section_pars = read_in_V96_cross_sections(V96_data_file)
    spectrum_E_log, spectrum_J_log = read_in_spectrum_file(pars["spectrum_file"])  
    spectrum_function = interpolate.interp1d(spectrum_E_log, spectrum_J_log)    #NB: In log-log space. Also, as a function of energy in Ryd.

    V95_data_file = "%s/cross_section_fits_verner95.dat" % (pars["cross_sections_data_path"], )
    x_section_pars_new, no_of_shells = read_in_V95_cross_sections(V95_data_file)
    KM93_data_file = "%s/auger_ionisation_probabilities_kaastra93.dat" % (pars["cross_sections_data_path"], ) 
    auger_probs = read_in_KM93_tables(KM93_data_file)
    shell_quantum_nos_V95 = x_section_pars_new[-1]

    # Parameters controlling the accuracy of the integration
    myEpsAbs = 1.0e-30
    myLimit = 300

    grey_section_denominator = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = myEpsAbs, limit=myLimit)[0]
    
    isotropic_photon_density = (4.0 * np.pi / 6.6260755e-27) * grey_section_denominator / 3.0e10  # cm^-3
    
    G0 = integrate.quad(lambda y: np.log(10) * (10.0 ** y) * 13.6 * 1.60217657e-12 * (10**spectrum_function(y)), np.log10(6.0 / 13.6), 0.0, epsabs = myEpsAbs, limit=myLimit)[0] * (4.0 * np.pi / 6.6260755e-27) / (3.0e10 * 5.29e-14)   # energy density 6->13.6 eV, in units of the Habing field 
    
    G0_parameter = G0 / (isotropic_photon_density * 3.0e10)
    
    H2_dissocJ = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(12.24 / 13.6), np.log10(13.51 / 13.6), epsabs = myEpsAbs, limit=myLimit)[0] * (4.0 * np.pi / (6.6260755e-27 * 3.0e10)) / (isotropic_photon_density * 3.0e10)

    He_thresh_energy = x_section_pars[1][:,0]
    C_thresh_energy = x_section_pars[2][:,0]
    N_thresh_energy = x_section_pars[3][:,0]
    O_thresh_energy = x_section_pars[4][:,0]
    Ne_thresh_energy = x_section_pars[5][:,0]
    Mg_thresh_energy = x_section_pars[6][:,0]
    Si_thresh_energy = x_section_pars[7][:,0]
    S_thresh_energy = x_section_pars[8][:,0]
    Ca_thresh_energy = x_section_pars[9][:,0]
    Fe_thresh_energy = x_section_pars[10][:,0]
    H_thresh_energy = np.zeros(3)
    H_thresh_energy_temp = x_section_pars[0][:,0]
    H_thresh_energy[0] = H_thresh_energy_temp[0]
    H_thresh_energy[1] = 0.755  # H-
    H_thresh_energy[2] = 15.4   # H2

    # The hydrogen arrays will be in the order:
    # HI, H-, H2
    
    if rank == 0: 
        print("Calculating cross sections")
        sys.stdout.flush()

        # Hydrogen and Helium 
        # use only rank 0. 
        H_cross_sections_root = np.zeros(3)
    
        # HI
        H_cross_sections_root[0] = (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(x_section_pars[0][0][0] / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) 
    
        # H- - for upper limit on integrals, just use HI's upper limit
        H_cross_sections_root[1] = (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_Hminus(13.6 * (10 ** y)), np.log10(0.755 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(0.755 / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) 

        # H2 + gamma -> H2+ + e-
        H_cross_sections_root[2] = (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_H2(13.6 * (10 ** y)), np.log10(15.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(15.4 / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) 

        print("H finished")
        sys.stdout.flush()

        He_cross_sections_root = np.zeros(2)
        for i in range(2):
            He_cross_sections_root[i] = (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(x_section_pars[1][i][0] / 13.6), np.log10(x_section_pars[1][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) 

        print("He finished")
        sys.stdout.flush()

    # For the remaining elements,
    # divide the work between
    # the MPI tasks. 
    C_cross_sections = np.zeros((6, N_Aug_C))
    for i in range(6):
        for k in range(N_Aug_C):
            work_idx = (i * N_Aug_C) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[0][i]):
                    if k == 0 and j == no_of_shells[0][i] - 1:      # For ionisation from the outermost shell, use verner96
                        C_cross_sections[i][k] += auger_probs[0][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        C_cross_sections[i][k] += (auger_probs[0][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[0][i, 1:6, (no_of_shells[0][i] - j - 1)], shell_quantum_nos_V95[0][i][no_of_shells[0][i] - j - 1][1]), np.log10(x_section_pars_new[0][i][0][(no_of_shells[0][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    C_cross_sections_root = comm.reduce(C_cross_sections, op = MPI.SUM, root = 0)
    
    if rank == 0: 
        print("C finished")
        sys.stdout.flush() 
    
    N_cross_sections = np.zeros((7, N_Aug_N))
    for i in range(7):
        for k in range(N_Aug_N):
            work_idx = (i * N_Aug_N) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[1][i]):
                    if k == 0 and j == no_of_shells[1][i] - 1:
                        N_cross_sections[i][k] += auger_probs[1][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        N_cross_sections[i][k] += (auger_probs[1][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[1][i, 1:6, (no_of_shells[1][i] - j - 1)], shell_quantum_nos_V95[1][i][no_of_shells[1][i] - j - 1][1]), np.log10(x_section_pars_new[1][i][0][(no_of_shells[1][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    N_cross_sections_root = comm.reduce(N_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("N finished")
        sys.stdout.flush() 
    
    O_cross_sections = np.zeros((8, N_Aug_O))
    for i in range(8):
        for k in range(N_Aug_O):
            work_idx = (i * N_Aug_O) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[2][i]):
                    if k == 0 and j == no_of_shells[2][i] - 1:
                        O_cross_sections[i][k] += auger_probs[2][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        O_cross_sections[i][k] += (auger_probs[2][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[2][i, 1:6, (no_of_shells[2][i] - j - 1)], shell_quantum_nos_V95[2][i][no_of_shells[2][i] - j - 1][1]), np.log10(x_section_pars_new[2][i][0][(no_of_shells[2][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    O_cross_sections_root = comm.reduce(O_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("O finished")
        sys.stdout.flush() 
    
    Ne_cross_sections = np.zeros((10, N_Aug_Ne))
    for i in range(10):
        for k in range(N_Aug_Ne):
            work_idx = (i * N_Aug_Ne) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[3][i]):
                    if k == 0 and j == no_of_shells[3][i] - 1:
                        Ne_cross_sections[i][k] += auger_probs[3][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        Ne_cross_sections[i][k] += (auger_probs[3][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[3][i, 1:6, (no_of_shells[3][i] - j - 1)], shell_quantum_nos_V95[3][i][no_of_shells[3][i] - j - 1][1]), np.log10(x_section_pars_new[3][i][0][(no_of_shells[3][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    Ne_cross_sections_root = comm.reduce(Ne_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Ne finished")
        sys.stdout.flush() 
    
    Mg_cross_sections = np.zeros((12, N_Aug_Mg))
    for i in range(12):
        for k in range(N_Aug_Mg):
            work_idx = (i * N_Aug_Mg) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[4][i]):
                    if k == 0 and j == no_of_shells[4][i] - 1:
                        Mg_cross_sections[i][k] += auger_probs[4][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        Mg_cross_sections[i][k] += (auger_probs[4][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[4][i, 1:6, (no_of_shells[4][i] - j - 1)], shell_quantum_nos_V95[4][i][no_of_shells[4][i] - j - 1][1]), np.log10(x_section_pars_new[4][i][0][(no_of_shells[4][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    Mg_cross_sections_root = comm.reduce(Mg_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Mg finished")
        sys.stdout.flush() 
    
    Si_cross_sections = np.zeros((14, N_Aug_Si))
    for i in range(14):
        for k in range(N_Aug_Si):
            work_idx = (i * N_Aug_Si) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[5][i]):
                    if k == 0 and j == no_of_shells[5][i] - 1:
                        Si_cross_sections[i][k] += auger_probs[5][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        Si_cross_sections[i][k] += (auger_probs[5][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[5][i, 1:6, (no_of_shells[5][i] - j - 1)], shell_quantum_nos_V95[5][i][no_of_shells[5][i] - j - 1][1]), np.log10(x_section_pars_new[5][i][0][(no_of_shells[5][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    Si_cross_sections_root = comm.reduce(Si_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Si finished")
        sys.stdout.flush() 
    
    S_cross_sections = np.zeros((16, N_Aug_S))            
    for i in range(16):
        for k in range(N_Aug_S):
            work_idx = (i * N_Aug_S) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[6][i]):
                    if k == 0 and j == no_of_shells[6][i] - 1:
                        S_cross_sections[i][k] += auger_probs[6][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        S_cross_sections[i][k] += (auger_probs[6][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[6][i, 1:6, (no_of_shells[6][i] - j - 1)], shell_quantum_nos_V95[6][i][no_of_shells[6][i] - j - 1][1]), np.log10(x_section_pars_new[6][i][0][(no_of_shells[6][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    S_cross_sections_root = comm.reduce(S_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("S finished")
        sys.stdout.flush() 
   
    Ca_cross_sections = np.zeros((20, N_Aug_Ca))
    for i in range(20):
        for k in range(N_Aug_Ca):
            work_idx = (i * N_Aug_Ca) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[7][i]):
                    if k == 0 and j == no_of_shells[7][i] - 1:
                        Ca_cross_sections[i][k] += auger_probs[7][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        Ca_cross_sections[i][k] += (auger_probs[7][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[7][i, 1:6, (no_of_shells[7][i] - j - 1)], shell_quantum_nos_V95[7][i][no_of_shells[7][i] - j - 1][1]), np.log10(x_section_pars_new[7][i][0][(no_of_shells[7][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    Ca_cross_sections_root = comm.reduce(Ca_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Ca finished")
        sys.stdout.flush() 
   
    Fe_cross_sections = np.zeros((26, N_Aug_Fe))
    for i in range(26):
        for k in range(N_Aug_Fe):
            work_idx = (i * N_Aug_Fe) + k 
            if work_idx % N_task == rank:
                for j in range(no_of_shells[8][i]):
                    if k == 0 and j == no_of_shells[8][i] - 1:
                        Fe_cross_sections[i][k] += auger_probs[8][i][j][k] * (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]) / (integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)), np.log10(x_section_pars[0][0][0] / 13.6), 10.0, epsabs = 1e-30, limit = myLimit)[0])
                    else:
                        Fe_cross_sections[i][k] += (auger_probs[8][i][j][k] * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_V95(13.6 * (10**y), x_section_pars_new[8][i, 1:6, (no_of_shells[8][i] - j - 1)], shell_quantum_nos_V95[8][i][no_of_shells[8][i] - j - 1][1]), np.log10(x_section_pars_new[8][i][0][(no_of_shells[8][i] - j - 1)] / 13.6), 10.0, epsabs = myEpsAbs, limit = myLimit)[0]) / grey_section_denominator

    comm.Barrier()
    Fe_cross_sections_root = comm.reduce(Fe_cross_sections, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Fe finished")
        sys.stdout.flush() 

    if rank == 0: 
        print("Calculating photoheating rates")
        sys.stdout.flush() 

        epsilon_H_root = np.zeros(3)
        
        # HI 
        epsilon_H_root[0] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(x_section_pars[0][0][0] / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(x_section_pars[0][0][0] / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[0][0][0]) 

        # H- 
        epsilon_H_root[1] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma_Hminus(13.6 * (10**y)), np.log10(0.755 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_Hminus(13.6 * (10**y)), np.log10(0.755 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - 0.755) 

        # H2 + gamma -> H2+ + e-
        epsilon_H_root[2] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma_H2(13.6 * (10**y)), np.log10(15.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma_H2(13.6 * (10**y)), np.log10(15.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - 15.4) 

        print("H finished")
        sys.stdout.flush() 

        epsilon_He_root = np.zeros(2) 
        for i in range(2):
            epsilon_He_root[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(x_section_pars[1][i][0] / 13.6), np.log10(x_section_pars[1][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(x_section_pars[1][i][0] / 13.6), np.log10(x_section_pars[1][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[1][i][0]) 

        print("He finished")
        sys.stdout.flush() 
            
    epsilon_C = np.zeros(6)
    for i in range(6):
        if i % N_task == rank: 
            epsilon_C[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[2][i][0]) 

    comm.Barrier()
    epsilon_C_root = comm.reduce(epsilon_C, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("C finished")
        sys.stdout.flush() 

    epsilon_N = np.zeros(7)
    for i in range(7):
        if i % N_task == rank: 
            epsilon_N[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[3][i][0]) 
        
    comm.Barrier()
    epsilon_N_root = comm.reduce(epsilon_N, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("N finished")
        sys.stdout.flush() 

    epsilon_O = np.zeros(8) 
    for i in range(8):
        if i % N_task == rank: 
            epsilon_O[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[4][i][0]) 

    comm.Barrier()
    epsilon_O_root = comm.reduce(epsilon_O, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("O finished")
        sys.stdout.flush() 

    epsilon_Ne = np.zeros(10)
    for i in range(10):
        if i % N_task == rank: 
            epsilon_Ne[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[5][i][0]) 

    comm.Barrier()
    epsilon_Ne_root = comm.reduce(epsilon_Ne, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Ne finished")
        sys.stdout.flush() 

    epsilon_Mg = np.zeros(12) 
    for i in range(12):
        if i % N_task == rank: 
            epsilon_Mg[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[6][i][0]) 

    comm.Barrier()
    epsilon_Mg_root = comm.reduce(epsilon_Mg, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Mg finished")
        sys.stdout.flush() 

    epsilon_Si = np.zeros(14)
    for i in range(14):
        if i % N_task == rank: 
            epsilon_Si[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[7][i][0])

    comm.Barrier()
    epsilon_Si_root = comm.reduce(epsilon_Si, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Si finished")
        sys.stdout.flush() 

    epsilon_S = np.zeros(16)
    for i in range(16):
        if i % N_task == rank: 
            epsilon_S[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[8][i][0]) 

    comm.Barrier()
    epsilon_S_root = comm.reduce(epsilon_S, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("S finished")
        sys.stdout.flush() 
    
    epsilon_Ca = np.zeros(20)
    for i in range(20):
        if i % N_task == rank: 
            epsilon_Ca[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[9][i][0]) 

    comm.Barrier()
    epsilon_Ca_root = comm.reduce(epsilon_Ca, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Ca finished")
        sys.stdout.flush() 
    
    epsilon_Fe = np.zeros(26) 
    for i in range(26):
        if i % N_task == rank: 
            epsilon_Fe[i] = 1.6022e-12 * (13.6*(((integrate.quad(lambda y: np.log(10) * (10**y) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0])/(integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]))) - x_section_pars[10][i][0]) 

    comm.Barrier()
    epsilon_Fe_root = comm.reduce(epsilon_Fe, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Fe finished")
        sys.stdout.flush() 
    
    # Each species has a total of nine factors that go into calculating
    # shielding. Three of these are used to calculate the suppressed 
    # photoionisation rate and six are used to calculate the effect
    # of spectral hardening. Three are tabulated in one dimension
    # per species and six are tabulated in two dimensions per species,
    # so we group them accordingly.
    if rank == 0: 
        print("Calculating shielding factors")
        sys.stdout.flush() 
    
    shield_H_1D = np.zeros((3, 3, len(N_x)))        # Dimensions: factor, ion, NHI. 

    # HI 
    Gamma_H0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(x_section_pars[0][0][0] / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
    for i in range(len(N_x)):
        if i % N_task == rank: 
            # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
            shield_H_1D[0][0][i] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(x_section_pars[0][0][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_H0
            # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
            shield_H_1D[1][0][i] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[0][0][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(x_section_pars[0][0][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
        
            # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
            shield_H_1D[2][0][i] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(x_section_pars[0][0][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_H_1D_root = comm.reduce(shield_H_1D, op = MPI.SUM, root = 0)
    
    shield_H_2D = np.zeros((6, 3, len(N_x), len(N_x)))       # Dimensions: factor, ion, N_H, N_He

    # HI 
    for i in range(len(N_x)):
        for j in range(len(N_x)):
            work_idx = (i * len(N_x)) + j
            if work_idx % N_task == rank: 
                # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                shield_H_2D[0][0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(15.4 / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_H0

                # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                shield_H_2D[1][0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(54.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_H0

                # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                shield_H_2D[2][0][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[0][0][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(15.4 / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                # Spectral hardening, numerator, E = 54.4 eV -> inf
                shield_H_2D[3][0][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[0][0][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(54.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                
                # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                shield_H_2D[4][0][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(15.4 / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                # Spectral hardening, denominator, E = 54.4 eV -> inf
                shield_H_2D[5][0][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[0][0][2:9], 0), np.log10(54.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    # H2
    Gamma_H20 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma_H2(13.6 * (10**y)), np.log10(15.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
    for i in range(len(N_x)):
        for j in range(len(N_x)):
            work_idx = (i * len(N_x)) + j
            if work_idx % N_task == rank: 
                # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                shield_H_2D[0][2][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma_H2(13.6 * (10**y)), np.log10(15.4 / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_H20

                # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                shield_H_2D[1][2][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma_H2(13.6 * (10**y)), np.log10(54.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_H20

                # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                shield_H_2D[2][2][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (15.4 / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma_H2(13.6 * (10**y)), np.log10(15.4 / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                # Spectral hardening, numerator, E = 54.4 eV -> inf
                shield_H_2D[3][2][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (15.4 / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma_H2(13.6 * (10**y)), np.log10(54.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                shield_H_2D[4][2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma_H2(13.6 * (10**y)), np.log10(15.4 / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                # Spectral hardening, denominator, E = 54.4 eV -> inf
                shield_H_2D[5][2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[i] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[j] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma_H2(13.6 * (10**y)), np.log10(54.4 / 13.6), np.log10(x_section_pars[0][0][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_H_2D_root = comm.reduce(shield_H_2D, op = MPI.SUM, root = 0)
                
    if rank == 0: 
        print("H finished")
        sys.stdout.flush() 

    # All He species have ionisation energies > 15.4 eV
    # so their 1D factors are all zero. 
    shield_He_1D_root = np.zeros((3, 2, len(N_x)))
    
    shield_He_2D = np.zeros((6, 2, len(N_x), len(N_x))) 
    for i in range(2):
        Gamma_He0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(x_section_pars[1][i][0] / 13.6), np.log10(x_section_pars[1][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
        for j in range(len(N_x)):
            for k in range(len(N_x)):
                work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                if work_idx % N_task == rank:
                    if x_section_pars[1][i][0] < He_thresh_energy[1]:       
                        # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                        shield_He_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(x_section_pars[1][i][0] / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_He0
                        
                    # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                    shield_He_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(54.4 / 13.6), np.log10(x_section_pars[1][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_He0

                    if x_section_pars[1][i][0] < He_thresh_energy[1]:       
                        # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                        shield_He_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[1][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(x_section_pars[1][i][0] / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                    # Spectral hardening, numerator, E = 54.4 eV -> inf
                    shield_He_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[1][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(54.4 / 13.6), np.log10(x_section_pars[1][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                    if x_section_pars[1][i][0] < He_thresh_energy[1]:       
                        # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                        shield_He_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(x_section_pars[1][i][0] / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                    # Spectral hardening, denominator, E = 54.4 eV -> inf
                    shield_He_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[1][i][2:9], 0), np.log10(54.4 / 13.6), np.log10(x_section_pars[1][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_He_2D_root = comm.reduce(shield_He_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("He finished")
        sys.stdout.flush() 
    
    shield_C_1D = np.zeros((3, 6, len(N_x)))
    for i in range(6):
        if x_section_pars[2][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_C0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank: 
                    if x_section_pars[2][i][0] < 15.4:      
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_C_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_C0
                    
                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_C_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[2][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_C_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_C_1D_root = comm.reduce(shield_C_1D, op = MPI.SUM, root = 0)
                    
    shield_C_2D = np.zeros((6, 6, len(N_x), len(N_x))) 
    for i in range(6):
        if x_section_pars[2][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_C0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(x_section_pars[2][i][0] / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[2][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_C_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(max(15.4, x_section_pars[2][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_C0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_C_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(max(54.4, x_section_pars[2][i][0]) / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_C0
                        
                        if x_section_pars[2][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_C_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[2][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(max(15.4, x_section_pars[2][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_C_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[2][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(max(54.4, x_section_pars[2][i][0]) / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[2][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_C_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(max(15.4, x_section_pars[2][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_C_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[2][i][2:9], shell_quantum_nos_V95[0][i][0][1]), np.log10(max(54.4, x_section_pars[2][i][0]) / 13.6), np.log10(x_section_pars[2][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_C_2D_root = comm.reduce(shield_C_2D, op = MPI.SUM, root = 0)

    if rank == 0:       
        print("C finished")
        sys.stdout.flush()
        
    shield_N_1D = np.zeros((3, 7, len(N_x)))
    for i in range(7):      
        if x_section_pars[3][i][0] > 13.6:        # below this threshold, shielding is by dust    
            Gamma_N0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank: 
                    if x_section_pars[3][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_N_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_N0

                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_N_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[3][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_N_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_N_1D_root = comm.reduce(shield_N_1D, op = MPI.SUM, root = 0)
    
    shield_N_2D = np.zeros((6, 7, len(N_x), len(N_x)))                     
    for i in range(7):      
        if x_section_pars[3][i][0] > 13.6:        # below this threshold, shielding is by dust    
            Gamma_N0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(x_section_pars[3][i][0] / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[3][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_N_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(max(15.4, x_section_pars[3][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_N0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_N_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(max(54.4, x_section_pars[3][i][0]) / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_N0

                        if x_section_pars[3][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_N_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[3][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(max(15.4, x_section_pars[3][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_N_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[3][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(max(54.4, x_section_pars[3][i][0]) / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[3][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_N_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(max(15.4, x_section_pars[3][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_N_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[3][i][2:9], shell_quantum_nos_V95[1][i][0][1]), np.log10(max(54.4, x_section_pars[3][i][0]) / 13.6), np.log10(x_section_pars[3][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_N_2D_root = comm.reduce(shield_N_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("N finished")
        sys.stdout.flush() 
            
    shield_O_1D = np.zeros((3, 8, len(N_x)))
    for i in range(8):      
        if x_section_pars[4][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_O0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank:
                    if x_section_pars[4][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_O_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_O0

                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_O_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[4][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_O_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_O_1D_root = comm.reduce(shield_O_1D, op = MPI.SUM, root = 0)

    shield_O_2D = np.zeros((6, 8, len(N_x), len(N_x)))
    for i in range(8):      
        if x_section_pars[4][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_O0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(x_section_pars[4][i][0] / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[4][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_O_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(max(15.4, x_section_pars[4][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_O0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_O_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(max(54.4, x_section_pars[4][i][0]) / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_O0

                        if x_section_pars[4][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_O_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[4][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(max(15.4, x_section_pars[4][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_O_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[4][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(max(54.4, x_section_pars[4][i][0]) / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
                        if x_section_pars[4][i][0] < He_thresh_energy[1]: 
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_O_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(max(15.4, x_section_pars[4][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_O_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[4][i][2:9], shell_quantum_nos_V95[2][i][0][1]), np.log10(max(54.4, x_section_pars[4][i][0]) / 13.6), np.log10(x_section_pars[4][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_O_2D_root = comm.reduce(shield_O_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("O finished")
        sys.stdout.flush() 
        
    shield_Ne_1D = np.zeros((3, 10, len(N_x)))
    for i in range(10):      
        if x_section_pars[5][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_Ne0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank:
                    if x_section_pars[5][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_Ne_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Ne0
                        
                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_Ne_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[5][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_Ne_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_Ne_1D_root = comm.reduce(shield_Ne_1D, op = MPI.SUM, root = 0)
    
    shield_Ne_2D = np.zeros((6, 10, len(N_x), len(N_x)))
    for i in range(10): 
        if x_section_pars[5][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_Ne0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(x_section_pars[5][i][0] / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[5][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_Ne_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(max(15.4, x_section_pars[5][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Ne0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_Ne_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(max(54.4, x_section_pars[5][i][0]) / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Ne0
                        
                        if x_section_pars[5][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_Ne_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[5][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(max(15.4, x_section_pars[5][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_Ne_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[5][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(max(54.4, x_section_pars[5][i][0]) / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[5][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_Ne_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(max(15.4, x_section_pars[5][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_Ne_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[5][i][2:9], shell_quantum_nos_V95[3][i][0][1]), np.log10(max(54.4, x_section_pars[5][i][0]) / 13.6), np.log10(x_section_pars[5][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_Ne_2D_root = comm.reduce(shield_Ne_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Ne finished")
        sys.stdout.flush() 
        
    shield_Mg_1D = np.zeros((3, 12, len(N_x)))
    for i in range(12):      
        if x_section_pars[6][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_Mg0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank:
                    if x_section_pars[6][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_Mg_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Mg0 

                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_Mg_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[6][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_Mg_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_Mg_1D_root = comm.reduce(shield_Mg_1D, op = MPI.SUM, root = 0)

    shield_Mg_2D = np.zeros((6, 12, len(N_x), len(N_x)))                    
    for i in range(12):      
        if x_section_pars[6][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_Mg0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(x_section_pars[6][i][0] / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[6][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_Mg_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(max(15.4, x_section_pars[6][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Mg0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_Mg_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(max(54.4, x_section_pars[6][i][0]) / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Mg0

                        if x_section_pars[6][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_Mg_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[6][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(max(15.4, x_section_pars[6][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_Mg_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[6][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(max(54.4, x_section_pars[6][i][0]) / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[6][i][0] < He_thresh_energy[1]: 
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_Mg_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(max(15.4, x_section_pars[6][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_Mg_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[6][i][2:9], shell_quantum_nos_V95[4][i][0][1]), np.log10(max(54.4, x_section_pars[6][i][0]) / 13.6), np.log10(x_section_pars[6][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_Mg_2D_root = comm.reduce(shield_Mg_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Mg finished")
        sys.stdout.flush() 
        
    shield_Si_1D = np.zeros((3, 14, len(N_x)))
    for i in range(14):      
        if x_section_pars[7][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_Si0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank:
                    if x_section_pars[7][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_Si_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Si0 

                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_Si_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[7][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_Si_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_Si_1D_root = comm.reduce(shield_Si_1D, op = MPI.SUM, root = 0)

    shield_Si_2D = np.zeros((6, 14, len(N_x), len(N_x)))                    
    for i in range(14):      
        if x_section_pars[7][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_Si0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(x_section_pars[7][i][0] / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[7][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_Si_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(max(15.4, x_section_pars[7][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Si0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_Si_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(max(54.4, x_section_pars[7][i][0]) / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Si0

                        if x_section_pars[7][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_Si_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[7][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(max(15.4, x_section_pars[7][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_Si_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[7][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(max(54.4, x_section_pars[7][i][0]) / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[7][i][0] < He_thresh_energy[1]: 
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_Si_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(max(15.4, x_section_pars[7][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_Si_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[7][i][2:9], shell_quantum_nos_V95[5][i][0][1]), np.log10(max(54.4, x_section_pars[7][i][0]) / 13.6), np.log10(x_section_pars[7][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_Si_2D_root = comm.reduce(shield_Si_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Si finished")
        sys.stdout.flush() 
            
    shield_S_1D = np.zeros((3, 16, len(N_x)))
    for i in range(16):      
        if x_section_pars[8][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_S0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank:
                    if x_section_pars[8][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_S_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_S0 

                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_S_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[8][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_S_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_S_1D_root = comm.reduce(shield_S_1D, op = MPI.SUM, root = 0)
                    
    shield_S_2D = np.zeros((6, 16, len(N_x), len(N_x)))
    for i in range(16): 
        if x_section_pars[8][i][0] > 13.6:        # below this threshold, shielding is by dust
            Gamma_S0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(x_section_pars[8][i][0] / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[8][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_S_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(max(15.4, x_section_pars[8][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_S0
                        
                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_S_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(max(54.4, x_section_pars[8][i][0]) / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_S0

                        if x_section_pars[8][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_S_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[8][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(max(15.4, x_section_pars[8][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_S_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[8][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(max(54.4, x_section_pars[8][i][0]) / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[8][i][0] < He_thresh_energy[1]: 
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_S_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(max(15.4, x_section_pars[8][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_S_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[8][i][2:9], shell_quantum_nos_V95[6][i][0][1]), np.log10(max(54.4, x_section_pars[8][i][0]) / 13.6), np.log10(x_section_pars[8][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

    comm.Barrier()
    shield_S_2D_root = comm.reduce(shield_S_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("S finished")
        sys.stdout.flush() 
        
    shield_Ca_1D = np.zeros((3, 20, len(N_x)))
    for i in range(20):      
        if x_section_pars[9][i][0] > 13.6:        # below this threshold, shielding is by dust    
            Gamma_Ca0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank:
                    if x_section_pars[9][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_Ca_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Ca0 
                        
                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_Ca_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[9][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_Ca_1D[2][i][j] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_Ca_1D_root = comm.reduce(shield_Ca_1D, op = MPI.SUM, root = 0)
                   
    shield_Ca_2D = np.zeros((6, 20, len(N_x), len(N_x)))
    for i in range(20):  
        if x_section_pars[9][i][0] > 13.6:        # below this threshold, shielding is by dust  
            Gamma_Ca0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(x_section_pars[9][i][0] / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[9][i][0] < He_thresh_energy[1]:     
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_Ca_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(max(15.4, x_section_pars[9][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Ca0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_Ca_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(max(54.4, x_section_pars[9][i][0]) / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Ca0

                        if x_section_pars[9][i][0] < He_thresh_energy[1]:     
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_Ca_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[9][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(max(15.4, x_section_pars[9][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_Ca_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[9][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(max(54.4, x_section_pars[9][i][0]) / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[9][i][0] < He_thresh_energy[1]: 
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_Ca_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(max(15.4, x_section_pars[9][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_Ca_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[9][i][2:9], shell_quantum_nos_V95[7][i][0][1]), np.log10(max(54.4, x_section_pars[9][i][0]) / 13.6), np.log10(x_section_pars[9][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_Ca_2D_root = comm.reduce(shield_Ca_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Ca finished")
        sys.stdout.flush() 


        
    shield_Fe_1D = np.zeros((3, 26, len(N_x)))
    for i in range(26):      
        if x_section_pars[10][i][0] > 13.6:       # below this threshold, shielding is by dust
            Gamma_Fe0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                work_idx = (i * len(N_x)) + j
                if work_idx % N_task == rank:
                    if x_section_pars[10][i][0] < 15.4: 
                        # Ionisation suppression from HI for E = 13.6 eV -> 15.4 eV
                        shield_Fe_1D[0][i][j] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Fe0 

                        # Spectral hardening, numerator, E = 13.6 eV -> 15.4 eV
                        shield_Fe_1D[1][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[10][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 13.6 eV -> 15.4 eV
                        shield_Fe_1D[2][i][j] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(15.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_Fe_1D_root = comm.reduce(shield_Fe_1D, op = MPI.SUM, root = 0)
    
    shield_Fe_2D = np.zeros((6, 26, len(N_x), len(N_x)))
    for i in range(26):      
        if x_section_pars[10][i][0] > 13.6:       # below this threshold, shielding is by dust
            Gamma_Fe0 = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(x_section_pars[10][i][0] / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
            for j in range(len(N_x)):
                for k in range(len(N_x)):
                    work_idx = (i * len(N_x) * len(N_x)) + (j * len(N_x)) + k
                    if work_idx % N_task == rank: 
                        if x_section_pars[10][i][0] < He_thresh_energy[1]:    
                            # Ionisation suppression from (HI + 3 * H2) and HeI for E = 15.4 eV -> 54.4 eV
                            shield_Fe_2D[0][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(max(15.4, x_section_pars[10][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Fe0

                        # Ionisation suppression from (HI + 3 * H2) and (HeI + 0.8 * HeII) for E = 54.4 eV -> inf
                        shield_Fe_2D[1][i][j][k] = integrate.quad(lambda y: np.log(10) * 4.0 * (np.pi / 6.62606957e-27) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(max(54.4, x_section_pars[10][i][0]) / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0] / Gamma_Fe0

                        if x_section_pars[10][i][0] < He_thresh_energy[1]:    
                            # Spectral hardening, numerator, E = 15.4 eV -> 54.4 eV
                            shield_Fe_2D[2][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[10][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(max(15.4, x_section_pars[10][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, numerator, E = 54.4 eV -> inf
                        shield_Fe_2D[3][i][j][k] = 1.6022e-12 * 13.6 * integrate.quad(lambda y: np.log(10) * ((10**y) - (x_section_pars[10][i][0] / 13.6)) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(max(54.4, x_section_pars[10][i][0]) / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        if x_section_pars[10][i][0] < He_thresh_energy[1]: 
                            # Spectral hardening, denominator, E = 15.4 eV -> 54.4 eV
                            shield_Fe_2D[4][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(max(15.4, x_section_pars[10][i][0]) / 13.6), np.log10(54.4 / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]

                        # Spectral hardening, denominator, E = 54.4 eV -> inf
                        shield_Fe_2D[5][i][j][k] = integrate.quad(lambda y: np.log(10) * (10**spectrum_function(y)) * np.exp(- N_x[j] * sigma_Hydrogen(13.6 * (10**y), x_section_pars[0][0][2:9], 0)) * np.exp(- N_x[k] * sigma_Helium(13.6 * (10**y), x_section_pars[1][0][2:9], 0)) * sigma(13.6 * (10**y), x_section_pars[10][i][2:9], shell_quantum_nos_V95[8][i][0][1]), np.log10(max(54.4, x_section_pars[10][i][0]) / 13.6), np.log10(x_section_pars[10][i][1] / 13.6), epsabs = myEpsAbs, limit = myLimit)[0]
                        
    comm.Barrier()
    shield_Fe_2D_root = comm.reduce(shield_Fe_2D, op = MPI.SUM, root = 0)

    if rank == 0: 
        print("Fe finished")
        sys.stdout.flush() 

    # For species with E_thresh < 13.6 eV, replace
    # gamma_HI with gamma_d from van Dishoeck et al. (2006). 
    # NOTE: CaI & CaII are not included, so we
    # use the species with the closest E_thresh
    # (MgI & CI respectively).
    if rank == 0: 
        shield_H_1D_root[0,1,0] = 0.5   # H-
        shield_C_1D_root[0,0,0] = 3.33   # CI
        shield_Mg_1D_root[0,0,0] = 2.08  # MgI
        shield_Si_1D_root[0,0,0] = 2.27  # SiI
        shield_S_1D_root[0,0,0] = 3.08   # SI
        shield_Ca_1D_root[0,0,0] = 2.08  # CaI
        shield_Ca_1D_root[0,1,0] = 3.33  # CaII
        shield_Fe_1D_root[0,0,0] = 2.20  # FeI

    # Protect against log10(0)
    if rank == 0: 
        shield_H_1D_root[(shield_H_1D_root < min_float_value)] = min_float_value 
        shield_He_1D_root[(shield_He_1D_root < min_float_value)] = min_float_value 
        shield_C_1D_root[(shield_C_1D_root < min_float_value)] = min_float_value 
        shield_N_1D_root[(shield_N_1D_root < min_float_value)] = min_float_value 
        shield_O_1D_root[(shield_O_1D_root < min_float_value)] = min_float_value 
        shield_Ne_1D_root[(shield_Ne_1D_root < min_float_value)] = min_float_value 
        shield_Mg_1D_root[(shield_Mg_1D_root < min_float_value)] = min_float_value 
        shield_Si_1D_root[(shield_Si_1D_root < min_float_value)] = min_float_value 
        shield_S_1D_root[(shield_S_1D_root < min_float_value)] = min_float_value 
        shield_Ca_1D_root[(shield_Ca_1D_root < min_float_value)] = min_float_value 
        shield_Fe_1D_root[(shield_Fe_1D_root < min_float_value)] = min_float_value
        
        shield_H_2D_root[(shield_H_2D_root < min_float_value)] = min_float_value 
        shield_He_2D_root[(shield_He_2D_root < min_float_value)] = min_float_value 
        shield_C_2D_root[(shield_C_2D_root < min_float_value)] = min_float_value 
        shield_N_2D_root[(shield_N_2D_root < min_float_value)] = min_float_value 
        shield_O_2D_root[(shield_O_2D_root < min_float_value)] = min_float_value 
        shield_Ne_2D_root[(shield_Ne_2D_root < min_float_value)] = min_float_value 
        shield_Mg_2D_root[(shield_Mg_2D_root < min_float_value)] = min_float_value 
        shield_Si_2D_root[(shield_Si_2D_root < min_float_value)] = min_float_value 
        shield_S_2D_root[(shield_S_2D_root < min_float_value)] = min_float_value 
        shield_Ca_2D_root[(shield_Ca_2D_root < min_float_value)] = min_float_value 
        shield_Fe_2D_root[(shield_Fe_2D_root < min_float_value)] = min_float_value

    # Arrange the cross sections, photoheating
    # rates and shielding factors into the
    # reaction groups from chimes_main_data.hdf5.
    if rank == 0:
        main_data_file = "../../../../chimes-data/chimes_main_data.hdf5"

        with h5py.File(pars["chimes_main_data_path"], "r") as h5_main_data:
            photoion_fuv_N = np.array(h5_main_data["photoion_fuv/N_reactions"])[1]
            photoion_fuv_reactants = np.array(h5_main_data["photoion_fuv/reactants"])
            
            photoion_euv_N = np.array(h5_main_data["photoion_euv/N_reactions"])[1]
            photoion_euv_reactants = np.array(h5_main_data["photoion_euv/reactants"])
            
            photoion_auger_fuv_N = np.array(h5_main_data["photoion_auger_fuv/N_reactions"])[1]
            photoion_auger_fuv_reactants = np.array(h5_main_data["photoion_auger_fuv/reactants"])
            photoion_auger_fuv_products = np.array(h5_main_data["photoion_auger_fuv/products"])
            
            photoion_auger_euv_N = np.array(h5_main_data["photoion_auger_euv/N_reactions"])[1]
            photoion_auger_euv_reactants = np.array(h5_main_data["photoion_auger_euv/reactants"])
            photoion_auger_euv_products = np.array(h5_main_data["photoion_auger_euv/products"])

        photoion_fuv_sigmaPhot = np.zeros(photoion_fuv_N, dtype = np.float32)
        photoion_fuv_epsilonPhot = np.zeros(photoion_fuv_N, dtype = np.float32)

        photoion_euv_sigmaPhot = np.zeros(photoion_euv_N, dtype = np.float32)
        photoion_euv_shieldFactor_1D = np.zeros((photoion_euv_N, 3, len(N_x)), dtype = np.float32)
        photoion_euv_shieldFactor_2D = np.zeros((photoion_euv_N, 6, len(N_x), len(N_x)), dtype = np.float32)

        photoion_auger_fuv_sigmaPhot = np.zeros(photoion_auger_fuv_N, dtype = np.float32)
        
        photoion_auger_euv_sigmaPhot = np.zeros(photoion_auger_euv_N, dtype = np.float32) 

        N_fuv_found = 0
        N_euv_found = 0 

        # HI
        reaction_idx = find_photoion_reaction(rc.chimes_dict["HI"], photoion_fuv_reactants)
        if reaction_idx < 0:
            reaction_idx = find_photoion_reaction(rc.chimes_dict["HI"], photoion_euv_reactants)
            photoion_euv_sigmaPhot[reaction_idx] = H_cross_sections_root[0]
            photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_H_1D_root[:, 0, :])
            photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_H_2D_root[:, 0, :, :])
            N_euv_found += 1
        else:
            photoion_fuv_sigmaPhot[reaction_idx] = H_cross_sections_root[0]
            photoion_fuv_epsilonPhot[reaction_idx] = epsilon_H_root[0]
            N_fuv_found += 1
            
        # H-
        reaction_idx = find_photoion_reaction(rc.chimes_dict["Hm"], photoion_fuv_reactants)
        if reaction_idx < 0:
            reaction_idx = find_photoion_reaction(rc.chimes_dict["Hm"], photoion_euv_reactants)
            photoion_euv_sigmaPhot[reaction_idx] = H_cross_sections_root[1]
            photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_H_1D_root[:, 1, :])
            photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_H_2D_root[:, 1, :, :])
            N_euv_found += 1
        else:
            photoion_fuv_sigmaPhot[reaction_idx] = H_cross_sections_root[1]
            photoion_fuv_epsilonPhot[reaction_idx] = epsilon_H_root[1]
            N_fuv_found += 1
            
        # H2
        reaction_idx = find_photoion_reaction(rc.chimes_dict["H2"], photoion_fuv_reactants)
        if reaction_idx < 0:
            reaction_idx = find_photoion_reaction(rc.chimes_dict["H2"], photoion_euv_reactants)
            photoion_euv_sigmaPhot[reaction_idx] = H_cross_sections_root[2]
            photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_H_1D_root[:, 2, :])
            photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_H_2D_root[:, 2, :, :])
            N_euv_found += 1
        else:
            photoion_fuv_sigmaPhot[reaction_idx] = H_cross_sections_root[2]
            photoion_fuv_epsilonPhot[reaction_idx] = epsilon_H_root[2]
            N_fuv_found += 1

        # Helium
        for ion_idx in range(2): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["HeI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["HeI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = He_cross_sections_root[ion_idx]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_He_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_He_2D_root[:, ion_idx, :, :])
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = He_cross_sections_root[ion_idx]
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_He_root[ion_idx]
                N_fuv_found += 1
                
        # Carbon
        for ion_idx in range(6): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["CI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["CI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = C_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_C_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_C_2D_root[:, ion_idx, :, :])
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = C_cross_sections_root[ion_idx, 0]
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_C_root[ion_idx]
                N_fuv_found += 1
                
        # Nitrogen
        for ion_idx in range(7): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["NI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["NI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = N_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_N_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_N_2D_root[:, ion_idx, :, :])
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = N_cross_sections_root[ion_idx, 0]
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_N_root[ion_idx]
                N_fuv_found += 1
                
        # Oxygen
        for ion_idx in range(8): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["OI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["OI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = O_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_O_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_O_2D_root[:, ion_idx, :, :])
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = O_cross_sections_root[ion_idx, 0]
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_O_root[ion_idx]
                N_fuv_found += 1
                
        # Neon
        for ion_idx in range(10): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["NeI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["NeI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = Ne_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_Ne_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_Ne_2D_root[:, ion_idx, :, :])
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = Ne_cross_sections_root[ion_idx, 0]
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_Ne_root[ion_idx]
                N_fuv_found += 1
                
        # Magnesium
        for ion_idx in range(12): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["MgI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["MgI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = Mg_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_Mg_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_Mg_2D_root[:, ion_idx, :, :])
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = Mg_cross_sections_root[ion_idx, 0]
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_Mg_root[ion_idx]
                N_fuv_found += 1
                
        # Silicon
        for ion_idx in range(14): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["SiI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["SiI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = Si_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_Si_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_Si_2D_root[:, ion_idx, :, :]) 
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = Si_cross_sections_root[ion_idx, 0] 
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_Si_root[ion_idx]
                N_fuv_found += 1
                
        # Sulphur
        for ion_idx in range(16): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["SI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["SI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = S_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_S_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_S_2D_root[:, ion_idx, :, :]) 
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = S_cross_sections_root[ion_idx, 0] 
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_S_root[ion_idx] 
                N_fuv_found += 1
                
        # Calcium
        for ion_idx in range(20): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["CaI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["CaI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = Ca_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_Ca_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_Ca_2D_root[:, ion_idx, :, :]) 
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = Ca_cross_sections_root[ion_idx, 0] 
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_Ca_root[ion_idx]
                N_fuv_found += 1
                
        # Iron
        for ion_idx in range(26): 
            reaction_idx = find_photoion_reaction(rc.chimes_dict["FeI"] + ion_idx, photoion_fuv_reactants)
            if reaction_idx < 0:
                reaction_idx = find_photoion_reaction(rc.chimes_dict["FeI"] + ion_idx, photoion_euv_reactants)
                photoion_euv_sigmaPhot[reaction_idx] = Fe_cross_sections_root[ion_idx, 0]
                photoion_euv_shieldFactor_1D[reaction_idx, :, :] = np.log10(shield_Fe_1D_root[:, ion_idx, :])
                photoion_euv_shieldFactor_2D[reaction_idx, :, :, :] = np.log10(shield_Fe_2D_root[:, ion_idx, :, :]) 
                N_euv_found += 1
            else:
                photoion_fuv_sigmaPhot[reaction_idx] = Fe_cross_sections_root[ion_idx, 0] 
                photoion_fuv_epsilonPhot[reaction_idx] = epsilon_Fe_root[ion_idx] 
                N_fuv_found += 1

        # Check that we have found the correct
        # number of photoion_fuv and
        # photoion_euv reactions.
        if N_fuv_found != photoion_fuv_N:
            raise Exception("ERROR: %d photoion_fuv reactions in chimes_main_data.hdf5, but %d photoion_fuv reactions were found in the cross sections. Aborting." % (photoion_fuv_N, N_fuv_found))
        
        if N_euv_found != photoion_euv_N:
            raise Exception("ERROR: %d photoion_euv reactions in chimes_main_data.hdf5, but %d photoion_euv reactions were found in the cross sections. Aborting." % (photoion_fuv_N, N_fuv_found))

        N_auger_fuv_found = 0
        N_auger_euv_found = 0 

        # Carbon Auger
        for ion_idx in range(6):
            for auger_idx in range(1, N_Aug_C): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["CI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["CI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = C_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = C_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Nitrogen Auger
        for ion_idx in range(7):
            for auger_idx in range(1, N_Aug_N): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["NI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["NI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = N_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = N_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Oxygen Auger
        for ion_idx in range(8):
            for auger_idx in range(1, N_Aug_O): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["OI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["OI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = O_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = O_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Neon Auger
        for ion_idx in range(10):
            for auger_idx in range(1, N_Aug_Ne): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["NeI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["NeI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = Ne_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = Ne_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Magnesium Auger
        for ion_idx in range(12):
            for auger_idx in range(1, N_Aug_Mg): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["MgI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["MgI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = Mg_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = Mg_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Silicon Auger
        for ion_idx in range(14):
            for auger_idx in range(1, N_Aug_Si): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["SiI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["SiI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = Si_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = Si_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Sulphur Auger
        for ion_idx in range(16):
            for auger_idx in range(1, N_Aug_S): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["SI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["SI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = S_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = S_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Calcium Auger
        for ion_idx in range(20):
            for auger_idx in range(1, N_Aug_Ca): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["CaI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["CaI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = Ca_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = Ca_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        # Iron Auger
        for ion_idx in range(26):
            for auger_idx in range(1, N_Aug_Fe): 
                reaction_idx = find_auger_reaction(rc.chimes_dict["FeI"] + ion_idx, auger_idx + 1, photoion_auger_fuv_reactants, photoion_auger_fuv_products)
                if reaction_idx < 0:
                    reaction_idx = find_auger_reaction(rc.chimes_dict["FeI"] + ion_idx, auger_idx + 1, photoion_auger_euv_reactants, photoion_auger_euv_products)
                    if reaction_idx >= 0: 
                        photoion_auger_euv_sigmaPhot[reaction_idx] = Fe_cross_sections_root[ion_idx, auger_idx]
                        N_auger_euv_found += 1
                else:
                    photoion_auger_fuv_sigmaPhot[reaction_idx] = Fe_cross_sections_root[ion_idx, auger_idx]
                    N_auger_fuv_found += 1
                    
        if N_auger_fuv_found != photoion_auger_fuv_N:
            raise Exception("ERROR: %d photoion_auger_fuv reactions in chimes_main_data.hdf5, but %d photoion_auger_fuv reactions were found in the cross sections. Aborting." % (photoion_auger_fuv_N, N_auger_fuv_found))
        
        if N_auger_euv_found != photoion_auger_euv_N:
            raise Exception("ERROR: %d photoion_auger_euv reactions in chimes_main_data.hdf5, but %d photoion_auger_euv reactions were found in the cross sections. Aborting." % (photoion_auger_euv_N, N_auger_euv_found))
    
    if rank == 0: 
        # Write arrays to output table
        with h5py.File(pars["output_file"], "w") as h5file:
            h5file["H2_dissocJ"] = np.array([H2_dissocJ], dtype = np.float32)
            h5file["G0_parameter"] = np.array([G0_parameter], dtype = np.float32)
            h5file["isotropic_photon_density"] = np.array([isotropic_photon_density], dtype = np.float32)

            h5file["photoion_fuv/sigmaPhot"] = photoion_fuv_sigmaPhot
            h5file["photoion_fuv/epsilonPhot"] = photoion_fuv_epsilonPhot
            h5file["photoion_fuv/N_reactions"] = photoion_fuv_N

            h5file["photoion_euv/sigmaPhot"] = photoion_euv_sigmaPhot
            h5file["photoion_euv/shieldFactor_1D"] = photoion_euv_shieldFactor_1D
            h5file["photoion_euv/shieldFactor_2D"] = photoion_euv_shieldFactor_2D
            h5file["photoion_euv/N_reactions"] = photoion_euv_N

            h5file["photoion_auger_fuv/sigmaPhot"] = photoion_auger_fuv_sigmaPhot
            h5file["photoion_auger_fuv/N_reactions"] = photoion_auger_fuv_N

            h5file["photoion_auger_euv/sigmaPhot"] = photoion_auger_euv_sigmaPhot
            h5file["photoion_auger_euv/N_reactions"] = photoion_auger_euv_N

            h5file["TableBins/Column_densities"] = np.log10(N_x, dtype = np.float32) 
            h5file["TableBins/N_Column_densities"] = len(N_x)
            
    comm.Barrier() 
        
    return

if __name__ == "__main__": 
    main()
