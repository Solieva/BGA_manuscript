The order in which to navigate this directory:

owlsnest_scripts
	- contains the scripts that were used to calculate the RMSF values from the MD simulations 
	- the output of these scripts are stored in this folder on an external drive because they are too large in size to upload: Passport_for_Mac/files_from_owlsnest_May16_2022/p16925/curating_dataset/rmsf_values

1_calculate_average_RMSF_values.ipynb
	- The output files from the owlsnest_scripts are loaded, then the  average RMSF values per residue and the average RMSF value for each of the 3 domains are calculated.  
	- These calculated values are stored in /average_RMSF_per_residue

2_graph_average_RMSF_per_residue_and_domains.ipynb
	- This notebook will graph the RMSF values per residue and per domain for temperatures 10C and 50C. 
	- the graph is in /graphs

3_calculate_and_graph_RMSF_of_surface_nonsurface_residues.ipynb
	- This notebook takes the surface and nonsurface residues lists from /surface_nonsurface_residues and calculates the average RMSF for these groups. 
	- It also graphs the average RMSFs of these groups
	- The % increase values are calculated betweeen surface and nonsurface residues, 1M and 4M, and 10C and 50C. 
	- The graphs are in /graphs and the RMSF values are printed inside of the notebook. 
