   This is a github repo of Chimes (originally developed by Alex Richings https://richings.bitbucket.io/chimes/home.html) which has all the dependencies in one place with new features added by Niranjan and Daniel.
   
   INSTALLATION NOTES (also very useful notes by Alex Richings here https://richings.bitbucket.io/chimes/user_guide/GettingStarted/index.html based on which these instructions are compiled): 
   - Clone the repo from https://github.com/NiranjanCRoy/chimes-IGNIS.git
   - Clone the data repo from https://bitbucket.org/richings/chimes-data.git inside chimes-IGNIS. This is large in size so couldn't be added to chimes-IGNIS and any new data can be updated as well by directly pulling from this repo.
	- In case the files are not being downloaded properly, perform a "git lfs install". Sometimes just the file pointers are downloaded when a repo with large files is cloned. 
   - untar the 'sundials-5.8.0.tar.gz' and 'cd' to the sundials-5.8.0 directory to run the following commands
   	- mkdir build
   	- cd build
   	- cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir/ -DBUILD_ARKODE=OFF -DBUILD_CVODE=ON -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=ON -DCMAKE_C_FLAGS="-O2" -DEXAMPLES_ENABLE_C=OFF -DSUNDIALS_PRECISION=double ../
   	- If you had an error about cmake version, do a "module load cmake"
   	- make
   	- make install
   	- export LD_LIBRARY_PATH=/path/to/install/dir/lib64:$LD_LIBRARY_PATH
   	- The last command can also be added to the .bashrc fmile to avoid having to run it each time while compiling chimes
   - Add the path to sundials installation directory in the chimes Makefile parameter "SUNDIALS_PATH" 
   	- hdf5 paths are left blank, if the system fails to automatically detect it, proper paths should be added
   	- run 'make'
   	- "libchimes.so" should be created which is the file we wanted to build.
   - We now need to install Meshoid (developed by Mike Grudic https://github.com/mikegrudic/meshoid.git) which is needed to calculate the density gradient needed for the Sobolev shielding length calculations if chosen.
   	-  cd to the meshoid directory available in chimes-IGNIS and execute the following commands
   	- python setup.py install --prefix=/desired/path/to/meshoidlib
   	- add this '/path/+/lib/python3.8/site-packages/meshoid-1.41-py3.8.egg' to the PYTHONPATH (and in bash to avoid having to run it each time while using chimes)
   
   
   TO RUN: 
   Inside chimes-driver directory:
   - Add the paths of 'libchimes.so' and chimes-data directory to the chimes parameter file.
        - Add the snapshot path to input_file parameter in the parameter file
   - Please manually provide the correct redshift of the snapshot in the parameter file.
   - To run in parallel:
   	- mpirun -np 40 python chimes-driver.py GIZMO_snapshot_eqm_state.param 
   - To run in serial:
	- mpirun -np 1 python chimes-driver.py GIZMO_snapshot_eqm_state.param 
		- Using mpirun because sometimes executables don't run without that on interactive nodes. One can also run "unset SLURM_JOBID" before to run the python command in a usual way without mpirun.  
