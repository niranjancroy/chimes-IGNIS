#ifdef CHIMES_ENABLE_GNU_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif

#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <math.h> 

#include "chimes_proto.h" 
#include "chimes_vars.h" 

int main(int argc, const char *argv[]) 
{ 
  int i; 
  struct gasVariables ChimesGasVars; 
  struct globalVariables ChimesGlobalVars; 

  /* Check the number of arguments */
  if (argc < 2) 
    {
      printf("Usage: chimes-debugger /path/to/chimes-data \n"); 
      exit(1); 
    }

  /* Define global variables. */ 
  sprintf(ChimesGlobalVars.MainDataTablePath, "%s/chimes_main_data.hdf5", argv[1]); 
  sprintf(ChimesGlobalVars.PhotoIonTablePath[0], "%s/cross_sections_B87.hdf5", argv[1]); 
  sprintf(ChimesGlobalVars.EqAbundanceTablePath, "None"); 
  
  ChimesGlobalVars.cellSelfShieldingOn = 1; 
  ChimesGlobalVars.N_spectra = 1;
  ChimesGlobalVars.redshift_dependent_UVB_index = -1; 
  ChimesGlobalVars.use_redshift_dependent_eqm_tables = 0; 
  ChimesGlobalVars.redshift = 0.0; 
  ChimesGlobalVars.reionisation_redshift = 0.0;   
  ChimesGlobalVars.StaticMolCooling = 1; 
  ChimesGlobalVars.T_mol = 1.0e5; 
  ChimesGlobalVars.grain_temperature = 10.0; 
  ChimesGlobalVars.cmb_temperature = 2.725; 
  ChimesGlobalVars.relativeTolerance = 1.0e-4; 
  ChimesGlobalVars.absoluteTolerance = 1.0e-10; 
  ChimesGlobalVars.explicitTolerance = 0.05; 
  ChimesGlobalVars.scale_metal_tolerances = 1; 
  ChimesGlobalVars.chimes_debug = 0; 
  ChimesGlobalVars.hybrid_cooling_mode = 0; 
    
  for (i = 0; i < 9; i++) 
    ChimesGlobalVars.element_included[i] = 1; 


  printf("Initialising CHIMES\n"); 
  init_chimes(&ChimesGlobalVars); 

  printf("Allocating memory to abundances array. \n"); 
  allocate_gas_abundances_memory(&ChimesGasVars, &ChimesGlobalVars); 

  /* Define gas variables. */ 
  double log_Z_over_Zsol = 0.0; 
  int n_iterations = 10; 

  ChimesGasVars.nH_tot = 1.0; 
  ChimesGasVars.temperature = 1.0e4; 
  ChimesGasVars.TempFloor = 10.1; 
  ChimesGasVars.divVel = 1.0; 
  ChimesGasVars.doppler_broad = 7.1;
  ChimesGasVars.cr_rate = 1.8e-16;   
  ChimesGasVars.hydro_timestep = 3.16e16; 
  ChimesGasVars.ForceEqOn = 0; 
  ChimesGasVars.ThermEvolOn = 0; 
  ChimesGasVars.InitIonState = 1; 
  ChimesGasVars.constant_heating_rate = 0.0; 
  
  // Compute element abundances  
  double atomic_masses[10] = {4., 12., 14., 16., 20., 24., 28., 32., 40., 56.}; 
  double solar_abundances[11] = {0.0129, 0.2806, 2.07e-3, 8.36e-4, 5.49e-3, 1.41e-3, 5.91e-4, 6.83e-4, 4.09e-4, 6.44e-5, 1.1e-3}; 
  double metallicity_array[11]; 
  
  for (i = 0; i < 11; i++) 
    metallicity_array[i] = solar_abundances[i] * pow(10.0, log_Z_over_Zsol); 
  
  // Account for a primordial He abundance 
  metallicity_array[1] = 0.25 + ((solar_abundances[1] - 0.25) * pow(10.0, log_Z_over_Zsol)); 

  double XH = 1.0 - (metallicity_array[0] + metallicity_array[1]);  

  for (i = 0; i < 10; i++) 
    ChimesGasVars.element_abundances[i] = metallicity_array[i + 1] / (atomic_masses[i] * XH); 
  
  ChimesGasVars.metallicity = metallicity_array[0] / 0.0129; 
  ChimesGasVars.dust_ratio = ChimesGasVars.metallicity; 

  /* Take the parameters for the strength 
   * of the radiation field from the 
   * cross-sections table. */ 
  ChimesGasVars.G0_parameter[0] = chimes_table_spectra.G0_parameter[0]; 
  ChimesGasVars.H2_dissocJ[0] = chimes_table_spectra.H2_dissocJ[0]; 
  ChimesGasVars.isotropic_photon_density[0] = chimes_table_spectra.isotropic_photon_density[0]; 

  /* In this example, we set the shielding 
   * length equal to the Jeans length, 
   * limited to a maximum of 100 kpc. */ 
  ChimesGasVars.cell_size = chimes_min(sqrt((5.0 / 3.0) * PI * BOLTZMANNCGS * ChimesGasVars.temperature * XH / (6.67408e-8 * (ChimesGasVars.nH_tot * PROTON_MASS / XH) * PROTON_MASS)), 3.086e23); 

  printf("Initialising gas abundances. \n"); 
  initialise_gas_abundances(&ChimesGasVars, &ChimesGlobalVars); 
  
  // Integrate chemistry 
  printf("Integrating chemistry...\n"); 
  for (i = 0; i < n_iterations; i++) 
    {
      printf("Iteration %d of %d\n", i + 1, n_iterations); 
      chimes_network(&ChimesGasVars, &ChimesGlobalVars); 
    }

  printf("Integration complete.\n"); 

  return 1; 
} 
