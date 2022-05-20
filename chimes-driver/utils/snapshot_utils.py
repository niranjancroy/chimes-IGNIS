import h5py
import numpy as np 
import os
import sys 
import time

sys.path.append('./ext-lib/pfh_python')
sys.path.append('./utils')
sys.path.append('/mnt/home/nroy/source_code/chimes-IGNIS/meshoid/lib_meshoid/lib/python3.8/site-packages/meshoid-1.41-py3.8.egg') #add path to meshoid installation directory here for your own machine

import gadget as g
from gizmopy.load_fire_snap import load_fire_snap 
from gizmopy.load_from_snapshot import load_from_snapshot 
from meshoid import *
import matplotlib.pyplot as plt

from chimes_utils import set_initial_chemistry_abundances 
from phys_const import proton_mass_cgs, boltzmann_cgs, seconds_in_a_Myr 
from shielding_utils import compute_jeans_shield_length, compute_colibre_shield_length, compute_sobolev_shield_length 


class SnapshotData: 
    def __init__(self, driver_pars, global_pars, gas_pars): 
        self.driver_pars = driver_pars 
        self.global_pars = global_pars 
        self.gas_pars = gas_pars 

        # Declare particle arrays that 
        # will need to be read in. 
        self.nH_arr = None 
        self.temperature_arr = None 
        self.metallicity_arr = None 
        self.shieldLength_arr = None 
        self.init_chem_arr = None 
        self.ChimesFluxIon_arr = None 
        self.ChimesFluxG0_arr = None 
        self.gas_coords_arr = None 
        self.star_coords_arr = None 
        self.star_mass_arr = None 
        self.star_age_Myr_arr = None 
        self.HIIregion_delay_time = None
        #niranjan 2021: adding the following
        self.gas_mass = None
        self.gas_hsml = None #smoothing length
        self.gas_density = None 
        return 

    def set_shielding_array(self): 
        if self.driver_pars['shield_mode'] == None: 
            self.shieldLength_arr = np.ones(len(self.nH_arr), dtype = np.float64) 
        elif self.driver_pars['shield_mode'] == "read-in": 
            # Reads in the shielding column density, from which 
            # it calculates the corresponding shielding lengths. 
            with h5py.File(self.driver_pars['input_file'], 'r') as h5file:
                self.shieldLength_arr = np.array(h5file[self.driver_pars["snapshot_column_density_array"]])/ self.nH_arr 
        elif self.driver_pars['shield_mode'] == "Jeans": 
            # Compute Jeans length assuming 
            # hydrogen mass fraction XH = 0.7 
            # mean molecular weight mu = 1 
            self.shieldLength_arr = np.zeros(len(self.nH_arr), dtype = np.float64) 
     
            for i in range(len(self.nH_arr)): 
                self.shieldLength_arr[i] = compute_jeans_shield_length(self.temperature_arr[i], self.nH_arr[i], self.driver_pars["shield_length_factor"], self.driver_pars["max_shield_length"]) 
                
        elif self.driver_pars['shield_mode'] == "Colibre": 
            self.shieldLength_arr = np.zeros(len(self.nH_arr), dtype = np.float64) 
            
            for i in range(len(self.nH_arr)): 
                XH = 1.0 - (self.metallicity_arr[i, 0] + self.metallicity_arr[i, 1]) 
                self.shieldLength_arr[i] = compute_colibre_shield_length(self.temperature_arr[i], self.nH_arr[i], XH, self.driver_pars["shield_length_factor"], self.driver_pars["max_shield_length"], self.driver_pars["colibre_log_T_min"], self.driver_pars["colibre_log_T_max"])


        #niranjan Dec2021: adding Sobolev self-shielding scheme
        elif self.driver_pars['shield_mode'] == "Sobolev":
             self.shieldLength_arr = np.zeros(len(self.nH_arr), dtype = np.float64)     
             print("Meshoid calculating density gradient...")        
             #print("Shape of Coord, Mass, HSML array = ", np.shape(self.gas_coords_arr), np.shape(self.gas_mass), np.shape(self.gas_hsml))
             M = Meshoid(self.gas_coords_arr, self.gas_mass, self.gas_hsml)
             density_gradient = M.D(self.gas_density) #calculating density gradient using meshoid
             density_gradient_magnitude = np.sqrt((density_gradient * density_gradient).sum(axis=1))       
             print("Density gradient calc done. Calculating Sobolev Shielding length...")      
             for i in range(len(self.nH_arr)):
                   self.shieldLength_arr[i] = compute_sobolev_shield_length(self.gas_density[i], density_gradient_magnitude[i], self.gas_hsml[i])

        else:
            raise Exception("ERROR: shield_mode %d not recognised. Aborting." % (self.driver_pars['shield_mode'], ))

        # HII regions 
        if self.driver_pars["disable_shielding_in_HII_regions"] == 1: 
            ind_HII = (self.HIIregion_delay_time > 0.0)
            #print('non-zero indices =', ind_HII.nonzero())
            self.shieldLength_arr[ind_HII] = 1.0 
        #print ('MIN shieldLength_arr =', np.min(self.shieldLength_arr), 'MAX shieldLength_arr =', np.max(self.shieldLength_arr)) 

        return




    def load_GIZMO(self):
        
        with h5py.File(self.driver_pars['input_file'], 'r') as h5file:
            # Define unit conversions. 
            hubble = h5file['Header'].attrs['HubbleParam']

            if self.driver_pars["snapshot_cosmo_flag"] == 0: 
                expansion_factor = 1.0 
            elif self.driver_pars["snapshot_cosmo_flag"] == 1: 
                expansion_factor = h5file['Header'].attrs['Time']

            # Gizmo units include factors of h^-1
            unit_mass_in_cgs = self.driver_pars["snapshot_unitMass_cgs"] / hubble  
            unit_length_in_cgs = self.driver_pars["snapshot_unitLength_cgs"] * expansion_factor / hubble  # co-moving to physical 
            unit_velocity_in_cgs = self.driver_pars["snapshot_unitVelocity_cgs"] 
            unit_internal_energy_in_cgs = unit_velocity_in_cgs ** 2.0 
            unit_time_in_cgs = unit_length_in_cgs / unit_velocity_in_cgs 

            print("Reading in particle data\n" )
            sys.stdout.flush() 

            # Read in metallicity array. 
            # Given as mass fractions relative 
            # to total in the order: 
            # All_metals, He, C, N, O, Ne, Mg, Si, S, Ca, Fe 
            self.metallicity_arr = np.array(h5file['PartType0/Metallicity'])

            # Calculate nH from the density array 
            density_arr = np.array(h5file['PartType0/Density'])
            XH = 1.0 - (self.metallicity_arr[:, 0] + self.metallicity_arr[:, 1]) 
            self.nH_arr = (unit_mass_in_cgs / (unit_length_in_cgs ** 3)) * density_arr * XH / proton_mass_cgs

            # Either read in the initial CHIMES 
            # abundance array, or set it by hand. 
            if self.driver_pars["snapshot_chemistry_array"] == None: 
                self.init_chem_arr = set_initial_chemistry_abundances(self.metallicity_arr, self.global_pars, self.gas_pars["InitIonState"]) 
            else: 
                try: 
                    self.init_chem_arr = np.array(h5file[self.driver_pars["snapshot_chemistry_array"]]) 
                except KeyError: 
                    raise Exception("ERROR: Chemistry array not found. The %s array is not present in the snapshot." % (self.driver_pars["snapshot_chemistry_array"], )) 

            # Calculate temperature from 
            # the internal energy array                     
            internal_energy_arr = np.array(h5file['PartType0/InternalEnergy']) 
            internal_energy_arr *= unit_internal_energy_in_cgs   # cgs 

            # If the simulation was run with CHIMES, we can use the mean molecular 
            # weights from the non-eqm chemical abundances to calculate the 
            # temperature. Otherwise, use mu assuming neutral gas. 
            try: 
                mmw_mu_arr = np.array(h5file['PartType0/ChimesMu']) 
            except KeyError: 
                helium_mass_fraction = self.metallicity_arr[:,1]
                y_helium = helium_mass_fraction / (4*(1-helium_mass_fraction))
                mmw_mu_arr = (1.0 + 4*y_helium) / (1+y_helium+ElectronAbundance) 
              
            self.temperature_arr = (2.0 / 3.0) * mmw_mu_arr * proton_mass_cgs * internal_energy_arr / boltzmann_cgs 

            # Read in stellar fluxes, if needed. 
            if self.driver_pars['UV_field'] == "StellarFluxes": 
                if self.driver_pars["compute_stellar_fluxes"] == 0: 
                    try: 
                        self.ChimesFluxIon_arr = np.array(h5file[self.driver_pars["snapshot_flux_ion_array"]]) 
                    except KeyError: 
                        raise Exception("ERROR: could not find array %s in the snapshot. You will need to compute the stellar fluxes, using the compute_stellar_fluxes parameter. Aborting." % (self.driver_pars["snapshot_flux_ion_array"], ))

                    try: 
                        self.ChimesFluxG0_arr = np.array(h5file[self.driver_pars["snapshot_flux_G0_array"]]) 
                    except KeyError: 
                        raise Exception("ERROR: could not find array %s in the snapshot. You will need to compute the stellar fluxes, using the compute_stellar_fluxes parameter. Aborting." % (self.driver_pars["snapshot_flux_G0_array"], ))
                elif self.driver_pars["compute_stellar_fluxes"] == 1: 
                    # Read in the star and gas particle data 
                    # needed to compute stellar fluxes. 
                    self.gas_coords_arr = np.array(h5file['PartType0/Coordinates']) * unit_length_in_cgs 

                    if self.driver_pars["snapshot_cosmo_flag"] == 0: 
                        time_Myr = h5file['Header'].attrs['Time'] * unit_time_in_cgs / seconds_in_a_Myr
                    
                        try:
                            coords_type4 = np.array(h5file['PartType4/Coordinates']) * unit_length_in_cgs
                            mass_type4 = np.array(h5file['PartType4/Masses']) * unit_mass_in_cgs
                            age_type4 = time_Myr - (np.array(h5file['PartType4/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr) 
                            self.star_coords_arr = coords_type4
                            self.star_mass_arr = mass_type4
                            self.star_age_Myr_arr = age_type4
                        except KeyError: 
                            print("Type 4 star particles are not present. Continuing.")
                            sys.stdout.flush()
                            self.star_coords_arr = np.empty((0, 3), dtype = np.float32) 
                            self.star_mass_arr = np.empty(0, dtype = np.float32) 
                            self.star_age_Myr_arr = np.empty(0, dtype = np.float32) 

                        try:
                            coords_type2 = np.array(h5file['PartType2/Coordinates']) * unit_length_in_cgs
                            mass_type2 = np.array(h5file['PartType2/Masses']) * unit_mass_in_cgs
                            age_type2 = time_Myr - (np.array(h5file['PartType2/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr)
                            self.star_coords_arr = np.concatenate((self.star_coords_arr, coords_type2)) 
                            self.star_mass_arr = np.concatenate((self.star_mass_arr, mass_type2)) 
                            self.star_age_Myr_arr = np.concatenate((self.star_age_Myr_arr, age_type2)) 
                        except KeyError: 
                            print("Type 2 star particles are not present. Continuing.")
                            sys.stdout.flush()

                        try:
                            coords_type3 = np.array(h5file['PartType3/Coordinates']) * unit_length_in_cgs
                            mass_type3 = np.array(h5file['PartType3/Masses']) * unit_mass_in_cgs
                            age_type3 = time_Myr - (np.array(h5file['PartType3/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr)
                            self.star_coords_arr = np.concatenate((self.star_coords_arr, coords_type3)) 
                            self.star_mass_arr = np.concatenate((self.star_mass_arr, mass_type3)) 
                            self.star_age_Myr_arr = np.concatenate((self.star_age_Myr_arr, age_type3)) 
                        except KeyError:                               
                            print("Type 3 star particles are not present. Continuing.")
                            sys.stdout.flush()

                    else: 
                        omega0 = h5file['Header'].attrs['Omega0'] 
                        H0_cgs = hubble * 3.2407789e-18            # Convert HubbleParam (in 100 km/s/Mpc) to s^-1 
                        try: 
                            self.star_coords_arr = np.concatenate((self.star_coords_arr, np.array(h5file['PartType4/Coordinates']) * unit_length_in_cgs)) 
                            self.star_mass_arr = np.concatenate((self.star_mass_arr, np.array(h5file['PartType4/Masses']) * unit_mass_in_cgs)) 
                            a_form = np.array(h5file['PartType4/StellarFormationTime']) 
                            x_form = (omega0 / (1.0 - omega0)) / (a_form ** 3.0) 
                            x_now = (omega0 / (1.0 - omega0)) / (expansion_factor ** 3.0) 
                            
                            self.star_age_Myr_arr = (2. / (3. * np.sqrt(1 - omega0))) * np.log(np.sqrt(x_form * x_now)/ ((np.sqrt(1.0 + x_now) - 1.0) * (np.sqrt(1.0 + x_form) + 1.0)))
                            self.star_age_Myr_arr /= H0_cgs 
                            self.star_age_Myr_arr /= seconds_in_a_Mpc 
                        except KeyError: 
                            print("Type 4 star particles are not present. Continuing.")
                            sys.stdout.flush()
                            self.star_coords_arr = np.empty((0, 3), dtype = np.float32) 
                            self.star_mass_arr = np.empty(0, dtype = np.float32) 
                            self.star_age_Myr_arr = np.empty(0, dtype = np.float32) 

                else: 
                    raise Exception("compute_stellar_fluxes == %d not recognised. Aborting." % (self.driver_pars["compute_stellar_fluxes"], )) 
            elif self.driver_pars["UV_field"] == "S04": 
                self.gas_coords_arr = np.array(h5file['PartType0/Coordinates']) * unit_length_in_cgs 

            if self.driver_pars["disable_shielding_in_HII_regions"] == 1: 
                try: 
                    self.HIIregion_delay_time = np.array(h5file[self.driver_pars["snapshot_HIIregion_array"]]) 
                except KeyError: 
                    raise Exception("ERROR: could not find array %s in the snapshot." % (self.driver_pars["snapshot_HIIregion_array"], )) 
            
            # Set the shielding length array 
            self.set_shielding_array() 

        return 

    
    
    
    
    
##################################################################################################
#niranjan-2021: adding a new function to load multiple hdf5 files in a single snapshot
##################################################################################################

    def load_GIZMO_MultiFile(self):
 
         #with g.readsnap(input_dir, snapnum, part_type, header_only=0) as h5file:
         #niranjan-2021: adding the directory and snapnum
         input_file = self.driver_pars['input_file']

         LastSlash = input_file.rindex('/')
         input_dir = input_file[:LastSlash+1:]
         
         LastUnderscore = input_file.rindex('_')
         snapnum  = int(input_file[LastUnderscore+1:]) 
 
         #with h5py.File(self.driver_pars['input_file'], 'r') as h5file:
         header = g.readsnap(input_dir, snapnum, 0, header_only=1)
        
         # Define unit conversions. 
         #hubble = h5file['Header'].attrs['HubbleParam']
         hubble = header['hubble']

         if self.driver_pars["snapshot_cosmo_flag"] == 0: 
             expansion_factor = 1.0 
         elif self.driver_pars["snapshot_cosmo_flag"] == 1: 
             #expansion_factor = h5file['Header'].attrs['Time']
             expansion_factor = header['time']

         # Gizmo units include factors of h^-1
         ##unit_mass_in_cgs = self.driver_pars["snapshot_unitMass_cgs"] / hubble  
         unit_mass_in_cgs = self.driver_pars["snapshot_unitMass_cgs"]

         ##unit_length_in_cgs = self.driver_pars["snapshot_unitLength_cgs"] * expansion_factor / hubble  # co-moving to physical 
         unit_length_in_cgs = self.driver_pars["snapshot_unitLength_cgs"]

         unit_velocity_in_cgs = self.driver_pars["snapshot_unitVelocity_cgs"] 
         unit_internal_energy_in_cgs = unit_velocity_in_cgs ** 2.0 
         unit_time_in_cgs = unit_length_in_cgs / unit_velocity_in_cgs 

         print("Reading in particle data\n" )
         sys.stdout.flush() 

         # Read in metallicity array. 
         # Given as mass fractions relative 
         # to total in the order: 
         # All_metals, He, C, N, O, Ne, Mg, Si, S, Ca, Fe 
         ##self.metallicity_arr = np.array(h5file['PartType0/Metallicity'])
         self.gas_mass = unit_mass_in_cgs * load_from_snapshot( 'Masses', 0, input_dir, snapnum)
         self.gas_hsml = unit_length_in_cgs * load_from_snapshot( 'SmoothingLength', 0, input_dir, snapnum)
         self.gas_density = (unit_mass_in_cgs / (unit_length_in_cgs ** 3)) *  load_from_snapshot( 'Density', 0, input_dir, snapnum)
         self.metallicity_arr = load_from_snapshot( 'Metallicity', 0, input_dir, snapnum)

         # Calculate nH from the density array 
         ##density_arr = np.array(h5file['PartType0/Density'])
         density_arr = load_from_snapshot( 'Density', 0, input_dir, snapnum)
 
         XH = 1.0 - (self.metallicity_arr[:, 0] + self.metallicity_arr[:, 1]) 
         self.nH_arr = (unit_mass_in_cgs / (unit_length_in_cgs ** 3)) * density_arr * XH / proton_mass_cgs

         # Either read in the initial CHIMES 
         # abundance array, or set it by hand. 
         if self.driver_pars["snapshot_chemistry_array"] == None: 
             self.init_chem_arr = set_initial_chemistry_abundances(self.metallicity_arr, self.global_pars, self.gas_pars["InitIonState"]) 
         else: 
             try: 
                 ##self.init_chem_arr = np.array(h5file[self.driver_pars["snapshot_chemistry_array"]]) 
                 self.init_chem_arr = load_from_snapshot( 'snapshot_chemistry_array', 0, input_dir, snapnum)
             except KeyError: 
                 raise Exception("ERROR: Chemistry array not found. The %s array is not present in the snapshot." % (self.driver_pars["snapshot_chemistry_array"], )) 

         # Calculate temperature from 
         # the internal energy array                     
         ##internal_energy_arr = np.array(h5file['PartType0/InternalEnergy']) 
         internal_energy_arr = load_from_snapshot( 'InternalEnergy', 0, input_dir, snapnum)
         internal_energy_arr *= unit_internal_energy_in_cgs   # cgs 

         # If the simulation was run with CHIMES, we can use the mean molecular 
         # weights from the non-eqm chemical abundances to calculate the 
         # temperature. Otherwise, use mu assuming neutral gas. 
         ##try: 
         ##    ##mmw_mu_arr = np.array(h5file['PartType0/ChimesMu']) 
         ##    mmw_mu_arr = load_from_snapshot('ChimesMu', 0, input_dir, snapnum)
         ##except KeyError:

         mmw_mu_arr = load_from_snapshot('ChimesMu', 0, input_dir, snapnum)
         if (np.size(mmw_mu_arr) == 1): #if 'ChimesMu' is not in the snapshot, load_from_snap will return 0, not an array
             print ('mmw_mu_arr not found in snapshot, calculating manually...') 

             helium_mass_fraction = self.metallicity_arr[:,1]
             #print('helium_mass_fraction min', np.min(helium_mass_fraction), 'helium_mass_fraction max', np.max(helium_mass_fraction))

             y_helium = helium_mass_fraction / (4*(1-helium_mass_fraction))
             ElectronAbundance = load_from_snapshot('ElectronAbundance', 0, input_dir, snapnum)
             #print('ElectronAbundance min', np.min(ElectronAbundance), 'ElectronAbundance max', np.max(ElectronAbundance))

             mmw_mu_arr = (1.0 + 4*y_helium) / (1+y_helium+ElectronAbundance) 


         self.temperature_arr = (2.0 / 3.0) * mmw_mu_arr * proton_mass_cgs * internal_energy_arr / boltzmann_cgs 

        # print ('internal_energy_arr min =', np.min(internal_energy_arr), 'internal_energy_arr max =', np.max(internal_energy_arr))
        # print ('mmw_mu_arr min', np.min(mmw_mu_arr), 'mmw_mu_arr max', np.max(mmw_mu_arr))
        # print ('temp min =', np.min(self.temperature_arr), 'temp max =', np.max(self.temperature_arr))
        # print ('ElectronAbundance min =', np.min(ElectronAbundance), 'ElectronAbundance max =', np.max(ElectronAbundance))
        # print ('helium_mass_fraction min =', np.min(helium_mass_fraction), 'helium_mass_fraction max =', np.max(helium_mass_fraction))

         # Read in stellar fluxes, if needed. 
         if self.driver_pars['UV_field'] == "StellarFluxes": 
             if self.driver_pars["compute_stellar_fluxes"] == 0: 
                 try: 
                     self.ChimesFluxIon_arr = np.array(h5file[self.driver_pars["snapshot_flux_ion_array"]]) 
                 except KeyError: 
                     raise Exception("ERROR: could not find array %s in the snapshot. You will need to compute the stellar fluxes, using the compute_stellar_fluxes parameter. Aborting." % (self.driver_pars["snapshot_flux_ion_array"], ))

                 try: 
                     self.ChimesFluxG0_arr = np.array(h5file[self.driver_pars["snapshot_flux_G0_array"]]) 
                 except KeyError: 
                     raise Exception("ERROR: could not find array %s in the snapshot. You will need to compute the stellar fluxes, using the compute_stellar_fluxes parameter. Aborting." % (self.driver_pars["snapshot_flux_G0_array"], ))
             elif self.driver_pars["compute_stellar_fluxes"] >= 1: 
                 # Read in the star and gas particle data 
                 # needed to compute stellar fluxes. 
                 ##self.gas_coords_arr = np.array(h5file['PartType0/Coordinates']) * unit_length_in_cgs 
                 self.gas_coords_arr = load_from_snapshot('Coordinates', 0, input_dir, snapnum) * unit_length_in_cgs 

                 if self.driver_pars["snapshot_cosmo_flag"] == 0: 
                     ##time_Myr = h5file['Header'].attrs['Time'] * unit_time_in_cgs / seconds_in_a_Myr
                     time_Myr = header['time'] * unit_time_in_cgs / seconds_in_a_Myr

                     try:
                         ##coords_type4 = np.array(h5file['PartType4/Coordinates']) * unit_length_in_cgs
                         coords_type4 = load_from_snapshot('Coordinates', 4, input_dir, snapnum) * unit_length_in_cgs

                         ##mass_type4 = np.array(h5file['PartType4/Masses']) * unit_mass_in_cgs
                         mass_type4 = load_from_snapshot('Masses', 4, input_dir, snapnum) * unit_mass_in_cgs

                         ##age_type4 = time_Myr - (np.array(h5file['PartType4/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr) 
                         age_type4 = time_Myr - (load_from_snapshot('StellarFormationTime', 4, input_dir, snapnum) * unit_time_in_cgs / seconds_in_a_Myr )

                         self.star_coords_arr = coords_type4
                         self.star_mass_arr = mass_type4
                         self.star_age_Myr_arr = age_type4
                     except KeyError: 
                         print("Type 4 star particles are not present. Continuing.")
                         sys.stdout.flush()
                         self.star_coords_arr = np.empty((0, 3), dtype = np.float32) 
                         self.star_mass_arr = np.empty(0, dtype = np.float32) 
                         self.star_age_Myr_arr = np.empty(0, dtype = np.float32) 

                     try:
                         ##coords_type2 = np.array(h5file['PartType2/Coordinates']) * unit_length_in_cgs
                         ##mass_type2 = np.array(h5file['PartType2/Masses']) * unit_mass_in_cgs
                         ##age_type2 = time_Myr - (np.array(h5file['PartType2/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr)
                         coords_type2 = load_from_snapshot('Coordinates', 2, input_dir, snapnum) * unit_length_in_cgs
                         mass_type2 = load_from_snapshot('Masses', 2, input_dir, snapnum) * unit_mass_in_cgs
                         age_type2 = time_Myr - (load_from_snapshot('StellarFormationTime', 2, input_dir, snapnum) * unit_time_in_cgs / seconds_in_a_Myr )

                         self.star_coords_arr = np.concatenate((self.star_coords_arr, coords_type2)) 
                         self.star_mass_arr = np.concatenate((self.star_mass_arr, mass_type2)) 
                         self.star_age_Myr_arr = np.concatenate((self.star_age_Myr_arr, age_type2)) 
                     except KeyError: 
                         print("Type 2 star particles are not present. Continuing.")
                         sys.stdout.flush()

                     try:
                         ##coords_type3 = np.array(h5file['PartType3/Coordinates']) * unit_length_in_cgs
                         ##mass_type3 = np.array(h5file['PartType3/Masses']) * unit_mass_in_cgs
                         ##age_type3 = time_Myr - (np.array(h5file['PartType3/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr)
                         coords_type3 = load_from_snapshot('Coordinates', 3, input_dir, snapnum) * unit_length_in_cgs
                         mass_type3 = load_from_snapshot('Masses', 3, input_dir, snapnum) * unit_mass_in_cgs
                         age_type3 = time_Myr - (load_from_snapshot('StellarFormationTime', 3, input_dir, snapnum) * unit_time_in_cgs / seconds_in_a_Myr )

                         self.star_coords_arr = np.concatenate((self.star_coords_arr, coords_type3)) 
                         self.star_mass_arr = np.concatenate((self.star_mass_arr, mass_type3)) 
                         self.star_age_Myr_arr = np.concatenate((self.star_age_Myr_arr, age_type3)) 
                     except KeyError:                               
                         print("Type 3 star particles are not present. Continuing.")
                         sys.stdout.flush()

                 else: 
                     ##omega0 = h5file['Header'].attrs['Omega0'] 
                     omega0 = 0.2821
                     #omega0 = Omega0 #from daa_constants #check for this using the load_from_snapshot
                     H0_cgs = hubble * 3.2407789e-18            # Convert HubbleParam (in 100 km/s/Mpc) to s^-1 
                     
                     try: 
                         #self.star_coords_arr = np.concatenate((self.star_coords_arr, load_from_snapshot('Coordinates', 4, input_dir, snapnum) * unit_length_in_cgs)) 
                         self.star_coords_arr = load_from_snapshot('Coordinates', 4, input_dir, snapnum) * unit_length_in_cgs
                         
                         #self.star_mass_arr = np.concatenate((self.star_mass_arr, load_from_snapshot('Masses', 4, input_dir, snapnum) * unit_mass_in_cgs)) 
                         self.star_mass_arr = load_from_snapshot('Masses', 4, input_dir, snapnum) * unit_mass_in_cgs
                         
                         #a_form = np.array(h5file['PartType4/StellarFormationTime'])
                         a_form = load_from_snapshot('StellarFormationTime', 4, input_dir, snapnum) 
                         x_form = (omega0 / (1.0 - omega0)) / (a_form ** 3.0) 
                         x_now = (omega0 / (1.0 - omega0)) / (expansion_factor ** 3.0) 

                         self.star_age_Myr_arr = (2. / (3. * np.sqrt(1 - omega0))) * np.log(np.sqrt(x_form * x_now)/ ((np.sqrt(1.0 + x_now) - 1.0) * (np.sqrt(1.0 + x_form) + 1.0)))
                         self.star_age_Myr_arr /= H0_cgs 
                         self.star_age_Myr_arr /= seconds_in_a_Myr #niranjan: previously seconds_in_a_Mpc(?)
                     except KeyError: 
                         print("Type 4 star particles are not present. Continuing.")
                         sys.stdout.flush()
                         self.star_coords_arr = np.empty((0, 3), dtype = np.float32) 
                         self.star_mass_arr = np.empty(0, dtype = np.float32) 
                         self.star_age_Myr_arr = np.empty(0, dtype = np.float32) 

             else: 
                 raise Exception("compute_stellar_fluxes == %d not recognised. Aborting." % (self.driver_pars["compute_stellar_fluxes"], )) 
         elif self.driver_pars["UV_field"] == "S04": 
             self.gas_coords_arr = load_from_snapshot('Coordinates', 0, input_dir, snapnum) * unit_length_in_cgs 

         if self.driver_pars["disable_shielding_in_HII_regions"] == 1: 
             #try: 
             #    self.HIIregion_delay_time = load_from_snapshot('snapshot_HIIregion_array', 0, input_dir, snapnum)
             #    #np.array(h5file[self.driver_pars["snapshot_HIIregion_array"]]) 
             #except: 
             #    #raise Exception("ERROR: could not find array %s in the snapshot." % (self.driver_pars["snapshot_HIIregion_array"], ))

             self.HIIregion_delay_time = load_from_snapshot('snapshot_HIIregion_array', 0, input_dir, snapnum)

             if (np.size(self.HIIregion_delay_time) == 1):
                 print("could not find array %s in the snapshot, calculating manually...\n." %(self.driver_pars["snapshot_HIIregion_array"],))
                 neutral_H = load_from_snapshot('NeutralHydrogenAbundance', 0, input_dir, snapnum)
                 self.HII_region_selection(neutral_H)                 
 
         #reading the black hole coordinates to define the center
         bh_center = load_from_snapshot('Coordinates', 5, input_dir, snapnum)
         filtering_radius = self.driver_pars["filtering_radius"]

         self.distance_filter(bh_center, filtering_radius)

         # Set the shielding length array 
         self.set_shielding_array() 


         #reading the black hole coordinates to define the center
         ##bh_center = load_from_snapshot('Coordinates', 5, input_dir, snapnum)
         ##filtering_radius = self.driver_pars["filtering_radius"] 

         #filtering the data within cutoff radius wrt black hole
         ##self.distance_filter(bh_center, filtering_radius)
         
         #print ("INDECES OF HII_REGION PARTICLES AFTER DISTANCE FILTERING : ", np.nonzero(self.HIIregion_delay_time)[0])         
         #print ("Shape of HIIregion_delay_time after distance filtering =", np.shape(self.HIIregion_delay_time), "non-zero items =", np.count_nonzero(self.HIIregion_delay_time))

         return 

    def HII_region_selection(self, neutral_H):
         '''selecting the gas particles based on temp~10,000 and HNumberDensity > 1 and IonizationFactor = 1'''
         mask = ((self.temperature_arr > 9e3) & (self.temperature_arr < 15e3) & (self.nH_arr > 0.1) & (neutral_H == 0))
         self.HIIregion_delay_time = mask.astype(int)
         return 



    def distance_filter(self, center, radius):
         ''' 'center' is the central coordinate, 'radius' is the radius of desired region in kpc. 
         IMPORTANT::When the black will be used as an ionization source, we need to convert its coordinate to be ZERO'''

         unit_length_in_cgs = self.driver_pars["snapshot_unitLength_cgs"]
         
         radius *= unit_length_in_cgs  #converting radius into CGS unit cm
         center *= unit_length_in_cgs
         
         self.gas_coords_arr -= center
         self.star_coords_arr -= center

         #creating mask for gas particles
         R_gas = np.sqrt((self.gas_coords_arr * self.gas_coords_arr).sum(axis=1))
         gas_mask = R_gas < radius

         #creating mask for star particles
         R_star = np.sqrt((self.star_coords_arr * self.star_coords_arr).sum(axis=1))
         star_mask = R_star < radius
    
         #applying the masks. If property has more than one column, applying mask just on the number of particles
         self.metallicity_arr = self.metallicity_arr[gas_mask,:] #multicolumn
         self.nH_arr = self.nH_arr[gas_mask]
         self.init_chem_arr = self.init_chem_arr[gas_mask,:] #multicolumn
         self.temperature_arr = self.temperature_arr[gas_mask]
         self.gas_coords_arr = self.gas_coords_arr[gas_mask,:] #multicolumn
         self.gas_mass = self.gas_mass[gas_mask]
         self.gas_density = self.gas_density[gas_mask]
         self.gas_hsml = self.gas_hsml[gas_mask]
         self.star_coords_arr = self.star_coords_arr[star_mask,:] #multicolumn
         self.star_mass_arr = self.star_mass_arr[star_mask]
         self.star_age_Myr_arr = self.star_age_Myr_arr[star_mask]

         if (self.driver_pars["disable_shielding_in_HII_regions"] == 1) : 
             self.HIIregion_delay_time  = self.HIIregion_delay_time[gas_mask]

         print ('Applied Distance Filter')         

         return

##################################################################################################


    def load_AREPO(self): 
        with h5py.File(self.driver_pars['input_file'], 'r') as h5file:
            # Define unit conversions. 
            hubble = h5file['Header'].attrs['HubbleParam']

            if self.driver_pars["snapshot_cosmo_flag"] == 0: 
                expansion_factor = 1.0 
            elif self.driver_pars["snapshot_cosmo_flag"] == 1: 
                redshift = h5file['Header'].attrs['Redshift'] 
                expansion_factor = 1.0 / (1.0 + redshift) 

            # Arepo units include factors of h^-1
            unit_mass_in_cgs = self.driver_pars["snapshot_unitMass_cgs"] / hubble  
            unit_length_in_cgs = self.driver_pars["snapshot_unitLength_cgs"] * expansion_factor / hubble  # co-moving to physical 
            unit_velocity_in_cgs = self.driver_pars["snapshot_unitVelocity_cgs"] 
            unit_internal_energy_in_cgs = unit_velocity_in_cgs ** 2.0  

            print("Reading in particle data\n" )
            sys.stdout.flush()

            # Read in metallicity array. 
            # Given as mass fractions relative 
            # to total in the order: 
            # H, He, C, N, O, Ne, Mg, Si, Fe, Other 
            GFM_metals = np.array(h5file["PartType0/GFM_Metals"]) 

            # Total metallicity 
            GFM_Z = np.array(h5file["PartType0/GFM_Metallicity"]) 
            
            # For chimes-driver, we the metallicity 
            # array to be in the order: 
            # All_metals, He, C, N, O, Ne, Mg, Si, S, Ca, Fe 
            N_part = len(GFM_Z) 
            self.metallicity_arr = np.zeros((N_part, 11)) 
            self.metallicity_arr[:, 0] = GFM_Z 
            self.metallicity_arr[:, 1] = GFM_metals[:, 1] 
            self.metallicity_arr[:, 2] = GFM_metals[:, 2] 
            self.metallicity_arr[:, 3] = GFM_metals[:, 3] 
            self.metallicity_arr[:, 4] = GFM_metals[:, 4] 
            self.metallicity_arr[:, 5] = GFM_metals[:, 5] 
            self.metallicity_arr[:, 6] = GFM_metals[:, 6] 
            self.metallicity_arr[:, 7] = GFM_metals[:, 7] 
            self.metallicity_arr[:, 10] = GFM_metals[:, 8] 

            # S and Ca are not explicitly tracked in 
            # Arepo. Assume Solar S/Si and Ca/Si ratios 
            self.metallicity_arr[:, 8] = GFM_metals[:, 7] * 0.464955 
            self.metallicity_arr[:, 9] = GFM_metals[:, 7] * 0.0964663 

            # Calculate nH from the density array 
            density_arr = np.array(h5file['PartType0/Density'])
            XH = GFM_metals[:, 0] 
            self.nH_arr = (unit_mass_in_cgs / (unit_length_in_cgs ** 3)) * density_arr * XH / proton_mass_cgs

            # Either read in the initial CHIMES 
            # abundance array, or set it by hand. 
            if self.driver_pars["snapshot_chemistry_array"] == None: 
                self.init_chem_arr = set_initial_chemistry_abundances(self.metallicity_arr, self.global_pars, self.gas_pars["InitIonState"]) 
            else: 
                try: 
                    self.init_chem_arr = np.array(h5file[self.driver_pars["snapshot_chemistry_array"]]) 
                except KeyError: 
                    raise Exception("ERROR: Chemistry array not found. The %s array is not present in the snapshot." % (self.driver_pars["snapshot_chemistry_array"], )) 

            # Calculate temperature from 
            # the internal energy array                     
            internal_energy_arr = np.array(h5file['PartType0/InternalEnergy']) 
            internal_energy_arr *= unit_internal_energy_in_cgs   # cgs 

            # If the simulation was run with CHIMES, we can use the mean molecular 
            # weights from the non-eqm chemical abundances to calculate the 
            # temperature. Otherwise, use mu assuming neutral gas. 
            try: 
                mmw_mu_arr = np.array(h5file['PartType0/ChimesMu']) 
            except KeyError: 
                helium_mass_fraction = self.metallicity_arr[:,1]
                y_helium = helium_mass_fraction / (4*(1-helium_mass_fraction))
                mmw_mu_arr = (1.0 + 4*y_helium) / (1+y_helium+ElectronAbundance) 
              
            self.temperature_arr = (2.0 / 3.0) * mmw_mu_arr * proton_mass_cgs * internal_energy_arr / boltzmann_cgs 

            # Read in stellar fluxes, if needed. 
            if self.driver_pars['UV_field'] == "StellarFluxes": 
                if driver_pars["compute_stellar_fluxes"] == 0: 
                    try: 
                        self.ChimesFluxIon_arr = np.array(h5file[self.driver_pars["snapshot_flux_ion_array"]]) 
                    except KeyError: 
                        raise Exception("ERROR: could not find array %s in the snapshot. You will need to compute the stellar fluxes, using the compute_stellar_fluxes parameter. Aborting." % (self.driver_pars["snapshot_flux_ion_array"], ))

                    try: 
                        self.ChimesFluxG0_arr = np.array(h5file[self.driver_pars["snapshot_flux_G0_array"]]) 
                    except KeyError: 
                        raise Exception("ERROR: could not find array %s in the snapshot. You will need to compute the stellar fluxes, using the compute_stellar_fluxes parameter. Aborting." % (self.driver_pars["snapshot_flux_G0_array"], ))
                    
            elif self.driver_pars["UV_field"] == "S04": 
                self.gas_coords_arr = np.array(h5file['PartType0/Coordinates']) * unit_length_in_cgs 

            if self.driver_pars["disable_shielding_in_HII_regions"] == 1: 
                try: 
                    self.HIIregion_delay_time = np.array(h5file[self.driver_pars["snapshot_HIIregion_array"]]) 
                except KeyError: 
                    raise Exception("ERROR: could not find array %s in the snapshot." % (self.driver_pars["snapshot_HIIregion_array"], )) 
                
            # Set the shielding length array 
            self.set_shielding_array() 

        return 

    def load_USER(self): 
        try: 
            from utils.user_def_utils import user_load_function 
        except ImportError: 
            raise Exception("ERROR: trying to load a user-defined snapshot file, but cannot find user_load_function in utils.user_def_utils. Aborting.") 

        # Read in snapshot data 
        user_load_function(self) 

        # Set the shielding length array 
        self.set_shielding_array() 

        return 
        
# Routine to check whether output 
# arrays already exist in the 
# output file, for IO_mode == snapshot. 
def snapshot_check_output_arrays(driver_pars): 
    return_value = False 

    if os.path.exists(driver_pars["output_file"]): 
        with h5py.File(driver_pars["output_file"], 'r') as h5file: 
            if driver_pars["driver_mode"] == "eqm_state": 
                array_name_list = ["%s/EqmChemistryAbundances" % (driver_pars["hdf5_output_group"], )] 
            elif driver_pars["driver_mode"] == "cooling_rates": 
                array_name_list = ["%s/log_cooling_rate" % (driver_pars["hdf5_output_group"], ), "%s/log_heating_rate" % (driver_pars["hdf5_output_group"], )] 
            elif driver_pars["driver_mode"] == "noneq_evolution": 
                array_name_list = ["%s/AbundanceEvolution" % (driver_pars["hdf5_output_group"], ), "%s/TemperatureEvolution" % (driver_pars["hdf5_output_group"], ), "%s/TimeArray_seconds" % (driver_pars["hdf5_output_group"], )]
            else: 
                raise Exception("ERROR: IO_mode %s cannot be run with driver_mode %s. Aborting." % (driver_pars["IO_mode"], driver_pars["driver_mode"])) 

            if driver_pars["UV_field"] == "StellarFluxes" and driver_pars["compute_stellar_fluxes"] == 1: 
                if driver_pars["snapshot_flux_ion_array"] != None: 
                    array_name_list.append(driver_pars["snapshot_flux_ion_array"]) 

                if driver_pars["snapshot_flux_G0_array"] != None: 
                    array_name_list.append(driver_pars["snapshot_flux_G0_array"]) 

            for array_name in array_name_list: 
                node_test = array_name in h5file 
                
                if node_test == True: 
                    return_value = True
                    print("%s already present in file %s" % (array_name, driver_pars["output_file"]))
                    sys.stdout.flush()

    return return_value 
