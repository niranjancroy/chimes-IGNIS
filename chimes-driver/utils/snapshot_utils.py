import h5py
import numpy as np 
import os
import sys 
sys.path.append('/home/jovyan/code/pfh_python')
sys.path.append('/home/jovyan/code/daa_python')
#sys.path.append('/home/jovyan/code/utilities')
#sys.path.append('/home/jovyan/code/gizmo_analysis')
import gadget as g
from gizmopy.load_fire_snap import load_fire_snap #go othrough the ocde
from gizmopy.load_from_snapshot import load_from_snapshot #go othrough the ocde
from daa_python import *
from daa_constants import *

from chimes_utils import set_initial_chemistry_abundances 
from phys_const import proton_mass_cgs, boltzmann_cgs, seconds_in_a_Myr 
from shielding_utils import compute_jeans_shield_length, compute_colibre_shield_length 


#niranjan-2021: adding the directory and snapnum
nofeedback_dir = '/home/jovyan/data/h113_HR_sn152/nof_s8e1_n128'
snapnum = 252
ptype = 0

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
        else: 
            raise Exception("ERROR: shield_mode %d not recognised. Aborting." % (self.driver_pars['shield_mode'], ))

        # HII regions 
        if self.driver_pars["disable_shielding_in_HII_regions"] == 1: 
            ind_HII = (self.HIIregion_delay_time > 0.0) 
            self.shieldLength_arr[ind_HII] = 1.0 
        
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
        
        #with g.readsnap(nofeedback_dir, snapnum, part_type, header_only=0) as h5file:
        
        #with h5py.File(self.driver_pars['input_file'], 'r') as h5file:
        header = g.readsnap(nofeedback_dir, snapnum, 0, header_only=1)
        
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
        self.metallicity_arr = load_from_snapshot( 'Metallicity', 0, nofeedback_dir, snapnum)

        # Calculate nH from the density array 
        ##density_arr = np.array(h5file['PartType0/Density'])
        density_arr = load_from_snapshot( 'Density', 0, nofeedback_dir, snapnum)

        XH = 1.0 - (self.metallicity_arr[:, 0] + self.metallicity_arr[:, 1]) 
        self.nH_arr = (unit_mass_in_cgs / (unit_length_in_cgs ** 3)) * density_arr * XH / proton_mass_cgs

        # Either read in the initial CHIMES 
        # abundance array, or set it by hand. 
        if self.driver_pars["snapshot_chemistry_array"] == None: 
            self.init_chem_arr = set_initial_chemistry_abundances(self.metallicity_arr, self.global_pars, self.gas_pars["InitIonState"]) 
        else: 
            try: 
                ##self.init_chem_arr = np.array(h5file[self.driver_pars["snapshot_chemistry_array"]]) 
                self.init_chem_arr = load_from_snapshot( 'snapshot_chemistry_array', 0, nofeedback_dir, snapnum)
            except KeyError: 
                raise Exception("ERROR: Chemistry array not found. The %s array is not present in the snapshot." % (self.driver_pars["snapshot_chemistry_array"], )) 

        # Calculate temperature from 
        # the internal energy array                     
        ##internal_energy_arr = np.array(h5file['PartType0/InternalEnergy']) 
        internal_energy_arr = load_from_snapshot( 'InternalEnergy', 0, nofeedback_dir, snapnum)
        internal_energy_arr *= unit_internal_energy_in_cgs   # cgs 

        # If the simulation was run with CHIMES, we can use the mean molecular 
        # weights from the non-eqm chemical abundances to calculate the 
        # temperature. Otherwise, use mu assuming neutral gas. 
        try: 
            ##mmw_mu_arr = np.array(h5file['PartType0/ChimesMu']) 
            mmw_mu_arr = load_from_snapshot('ChimesMu', 0, nofeedback_dir, snapnum)
        except KeyError: 
            helium_mass_fraction = self.metallicity_arr[:,1]
            y_helium = helium_mass_fraction / (4*(1-helium_mass_fraction))
            ElectronAbundance = load_from_snapshot('ElectronAbundance', 0, nofeedback_dir, snapnum)
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
                ##self.gas_coords_arr = np.array(h5file['PartType0/Coordinates']) * unit_length_in_cgs 
                self.gas_coords_arr = load_from_snapshot('Coordinates', 0, nofeedback_dir, snapnum) * unit_length_in_cgs 

                if self.driver_pars["snapshot_cosmo_flag"] == 0: 
                    ##time_Myr = h5file['Header'].attrs['Time'] * unit_time_in_cgs / seconds_in_a_Myr
                    time_Myr = header['time'] * unit_time_in_cgs / seconds_in_a_Myr

                    try:
                        ##coords_type4 = np.array(h5file['PartType4/Coordinates']) * unit_length_in_cgs
                        coords_type4 = load_from_snapshot('Coordinates', 4, nofeedback_dir, snapnum) * unit_length_in_cgs

                        ##mass_type4 = np.array(h5file['PartType4/Masses']) * unit_mass_in_cgs
                        mass_type4 = load_from_snapshot('Masses', 4, nofeedback_dir, snapnum) * unit_mass_in_cgs

                        ##age_type4 = time_Myr - (np.array(h5file['PartType4/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr) 
                        age_type4 = time_Myr - (load_from_snapshot('StellarFormationTime', 4, nofeedback_dir, snapnum) * unit_time_in_cgs / seconds_in_a_Myr )

                        self.star_coords_arr = coords_type4
                        self.star_mass_arr = mass_type4
                        self.star_age_Myr_arr = age_type4
                    except KeyError: 
                        print("Type 4 star particles are not present. Continuing.")
                        sys.stdout.flush()
                        self.star_coords_arr = np.empty((0, 3), dtype = np.float32) 
                        self.star_mass_arr = np.empty(0, dtype = np.float32) 
                        self.star_age_Myr_arr = np.empty(0, dtype = np.float32) 

##                    try:
##                        ##coords_type2 = np.array(h5file['PartType2/Coordinates']) * unit_length_in_cgs
##                        ##mass_type2 = np.array(h5file['PartType2/Masses']) * unit_mass_in_cgs
##                        ##age_type2 = time_Myr - (np.array(h5file['PartType2/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr)
##                        coords_type2 = load_from_snapshot('Coordinates', 2, nofeedback_dir, snapnum) * unit_length_in_cgs
##                        mass_type2 = load_from_snapshot('Masses', 2, nofeedback_dir, snapnum) * unit_mass_in_cgs
##                        age_type2 = time_Myr - (load_from_snapshot('StellarFormationTime', 2, nofeedback_dir, snapnum) * unit_time_in_cgs / seconds_in_a_Myr )
##
##                        self.star_coords_arr = np.concatenate((self.star_coords_arr, coords_type2)) 
##                        self.star_mass_arr = np.concatenate((self.star_mass_arr, mass_type2)) 
##                        self.star_age_Myr_arr = np.concatenate((self.star_age_Myr_arr, age_type2)) 
##                    except KeyError: 
##                        print("Type 2 star particles are not present. Continuing.")
##                        sys.stdout.flush()
##
##                    try:
##                        ##coords_type3 = np.array(h5file['PartType3/Coordinates']) * unit_length_in_cgs
##                        ##mass_type3 = np.array(h5file['PartType3/Masses']) * unit_mass_in_cgs
##                        ##age_type3 = time_Myr - (np.array(h5file['PartType3/StellarFormationTime']) * unit_time_in_cgs / seconds_in_a_Myr)
##                        coords_type3 = load_from_snapshot('Coordinates', 3, nofeedback_dir, snapnum) * unit_length_in_cgs
##                        mass_type3 = load_from_snapshot('Masses', 3, nofeedback_dir, snapnum) * unit_mass_in_cgs
##                        age_type3 = time_Myr - (load_from_snapshot('StellarFormationTime', 3, nofeedback_dir, snapnum) * unit_time_in_cgs / seconds_in_a_Myr )
##
##                        self.star_coords_arr = np.concatenate((self.star_coords_arr, coords_type3)) 
##                        self.star_mass_arr = np.concatenate((self.star_mass_arr, mass_type3)) 
##                        self.star_age_Myr_arr = np.concatenate((self.star_age_Myr_arr, age_type3)) 
##                    except KeyError:                               
##                        print("Type 3 star particles are not present. Continuing.")
##                        sys.stdout.flush()

                else: 
                    ##omega0 = h5file['Header'].attrs['Omega0'] 
                    omega0 = Omega0 #from daa_constants #check for this using the load_from_snapshot
                    H0_cgs = hubble * 3.2407789e-18            # Convert HubbleParam (in 100 km/s/Mpc) to s^-1 
                    
                    try: 
                        #self.star_coords_arr = np.concatenate((self.star_coords_arr, load_from_snapshot('Coordinates', 4, nofeedback_dir, snapnum) * unit_length_in_cgs)) 
                        self.star_coords_arr = load_from_snapshot('Coordinates', 4, nofeedback_dir, snapnum) * unit_length_in_cgs
                        
                        #self.star_mass_arr = np.concatenate((self.star_mass_arr, load_from_snapshot('Masses', 4, nofeedback_dir, snapnum) * unit_mass_in_cgs)) 
                        self.star_mass_arr = load_from_snapshot('Masses', 4, nofeedback_dir, snapnum) * unit_mass_in_cgs
                        
                        #a_form = np.array(h5file['PartType4/StellarFormationTime'])
                        a_form = load_from_snapshot('StellarFormationTime', 4, nofeedback_dir, snapnum) 
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
            self.gas_coords_arr = load_from_snapshot('Coordinates', ptype, nofeedback_dir, snapnum) * unit_length_in_cgs 

        if self.driver_pars["disable_shielding_in_HII_regions"] == 1: 
            try: 
                self.HIIregion_delay_time = np.array(h5file[self.driver_pars["snapshot_HIIregion_array"]]) 
            except KeyError: 
                raise Exception("ERROR: could not find array %s in the snapshot." % (self.driver_pars["snapshot_HIIregion_array"], ))

        # Set the shielding length array 
        self.set_shielding_array() 

        return 

###################################################################################################



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
