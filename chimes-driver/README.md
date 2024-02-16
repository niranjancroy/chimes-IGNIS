# CHIMES Driver

Welcome to the CHIMES Driver repository. This python script can be used to run CHIMES as a stand-alone module.

For more information on CHIMES, please see the [main CHIMES home page](https://richings.bitbucket.io/chimes/home.html).

For more information on how to use CHIMES Driver, please see the [CHIMES User Guide](https://richings.bitbucket.io/chimes/user_guide/index.html).


New features added in chimes-driver by Niranjan and Daniel:
- Option to read-in multifile snapshots from FIRE-like simulations
	- By selecting 'GIZMO_MultiFile' as the 'snapshot_type' in parameter file
- HII region selection when not available from simulation data
	- If disable_shielding_in_HII_regions is set to 1 and HII region array is not available in the simulation data
- Sobolev approximation of shielding length
	- If shield_mode is set to "Sobolev"
- Radius filtering option
	- By setting filtering_radius in kpc in the parameter file.
- Significantly faster incident flux calculation from star particles using tree method (Using algorithms from Mike Grudic's pytreegrav https://github.com/mikegrudic/pytreegrav.git)
	- By setting "compute_stellar_fluxes"  to  '2'  ( '1'-Chimes default; '2'-Using Tree; '3'-Using faster bruteforce)


