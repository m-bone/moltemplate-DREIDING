The xyzToMoltemplate tool converts an .xyz molecule file to moltemplate .lt files. It produces a 'molecule'.lt file and system.lt file. Output files follow the LAMMPS 'full' unit style and are built with cubic unit cells. 

The tool is typically run with run_xyz_to_moltemplate.sh and takes three user arguments. The file directory is found automatically; the tool is designed to be accessable globally within a Linux OS so that the user can call the tool from their desired directory. If not run globally, all these files will need to be in the same directory as the .xyz file.
The first argument is the .xyz file name without the extension e.g. methanol
The second argument is the cubic unit cell size in Angstroms e.g. 20.0
The third argument is the forcefield without an extension e.g. dreiding

Bonds are found by assigning atoms within 1.6 Angstroms of each other as bonded. If your system has free atoms within 1.6 Angstroms of a neighbour, this tool will cause errors. This can be changed by expanding bondCalculation in xyzFunctions.py

The atom dictionary is not exhaustive and can be expanded by adding atoms and masses to the atomDictionary in xyzFunctions.py. The code will return an error if it encounters unknown atoms.

This tool does not come with any warrenty but has been extensively tested and is still in regular use. 
