The order in which to navigate this directory:

owlsnest_scripts
	- contains the scripts that were used to calculate the SASA values from the MD simulations 
	- the output of these scripts are stored in this folder on an external drive because they are too large in size to upload: Passport_for_Mac/files_from_owlsnest_May16_2022/p16925/curating_dataset/SASA_calculations_May15_2022/SASA_data_files

1_calculate_total_SASA.ipynb
	-  The output files from the owlsnest_scripts are loaded, the average SASA per atom per trajectory is calculated, then the sum is computed for a final, total SASA value per molarity and temperature system. 
	- The total SASA values are stored in /total_SASA_values_per_system

2_graph_total_SASA.ipynb 
	- This notebook graphs the total SASA values per molarity and temperature
	- The graph is saved in /graphs

3_surface_and_nonsurface_residues.ipynb 
	- This notebook will take the average of the SASA values per atom per trajectory, find the surface exposed residues, and save them into a numpy file. Then, it will save the non-surface residues into a numpy file. The non-surface residues are any residues that don't show up in the surface exposed residue list. 
	- The lists of surface and nonsurface residues per system are in /surface_nonsurface_residues

4_graph_number_of_surface_nonsurface_residues.ipynb
	- This notebook graphs the number of surface and nonsurface residues in each system
	- The graphs are saved in /graphs
