Scoring function GUI, used to determine the crystal structure from a 
potential ensemble of crystal structures with the aid of solution chemical 
shifts. A full guide to the GUI can be seen in Documentation, with a 
flowchart describing the procedure and which section of the GUI to use. 

A list of labels that match the user labels (used in experiments) compared 
to the ones generated by magres must be given, see Magres labels .csv for 
an example.

Magres files for both solution and solid state are provided for furosemide 
to be used as an example.

Chemical shifts can be extracted from the magres file by using the 
chemical shifts GUI.

Calculated solution shifts would need to be averaged once all the atoms 
have converged, which can be tested using the solution ensemble GUI. The 
output is provided as Calculated_solution_chemical_shits.csv.

Experimental solution chemical shifts are provided also, with solid 
chemical shifts obtained from Widdifield et al, 2016 (DOI	
https://doi.org/10.1039/C6CC02171A).

Using the calculated and experimental solution and solid chemical shifts, 
the scoring function can be performed in the scoring function GUI.
