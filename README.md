# VASPBAUM 
BAnd Unfolding Machinery for VASP output (VASP-BAUM).
VASPBAUM is written for the post-processing purpose of the VASP outputs, i.e., WAVECAR the Bloch wavefunction information. VASPBAUM can compute band structure unfolding via k-projection method [See described in Phys. Rev. B 85, 085201 (2012).] In addition Circular dichroism also can be evaluated. 

# Compile
* Serial version : set MPI_USE = NO in makefile
    > make vaspbaum.serial or
	> make
* Multicore version : set MPI_USE = YES in makefile
    > make vaspbaum.mpi     or
	> make

# Features
* Band structure unfolding along primitive unit cell brillouin zone k-path
* Circular dichroism (optical selectivity response to the circulary polarized light)

# Usage
* Instruction and possible options
> ./vaspbaum.mpi -h

* generate supercell K-point that primitive k-point along high symmetry path is folded into.
> 1. Prepare KPOINTS_PC : k-point path of primitive Brillouin zone (BZ) that to be unfolded onto
> 2. Prepare POSCAR_PC : primitive cell lattice information (VASP POSCAR format)
> 3. Prepare POSCAR_SC : Super     cell lattice information (VASP POSCAR format)
> 4. run "vaspbaum" : vaspbaum -set_unfold
> 5. KPOINTS_SC file will be generated --> copy KPOINTS_SC into KPOINTS to calculate WAVECAR

* run unfolding
> 1. Once you generated "WAVECAR" with KPOINTS_SC
> 2. run "vaspbaum" to unfold:  vaspbaum -unfold -nosoc -sigma 0.10 -nediv 4000 -norm T -ef -1.2
>    For example, -nosoc : without SOC (if SOC, set -soc)
>				  -nediv : division of energy window for spectral function
>				  -norm  : normalize wavefunction
>				  -ef    : set Fermi level       
>				  -sigma : full width at half maximum for the Lorenztian line shape function (used for smearing of spectral function)

* run unfolding with circular dichroism
> same with "unfolding" functionality but include cicularly polarization 
> 1. run "vaspbaum" : vaspbaum -unfold -nosoc -sigma 0.10 -nediv 4000 -norm T -cd 2  -ien 0.0 -fen 4.0
> 				  -ien   : incident photon energy (minimum)
> 				  -fen   : incident photon energy (maximum)
> 2. You can also controle the angle of incident photon with respect to the surface by -theta and -phi tag
> 		    : vaspbaum -unfold -nosoc -sigma 0.10 -nediv 4000 -norm T -cd 2  -ien 0.0 -fen 4.0 -theta 65.0
> 		                  -theta : angle along z-axis describing the direction of the injecting light
> 		                  -phi   : angle along x-axis describing the direction of the injecting light
> 		                   * Example:
> 		                   Light from z-axis(surface normal): (theta,phi) = (0.0,0.0)
>               		   Light from x-axis                : (theta,phi) = (90.0,0.0)
>              		      	   Light from y-axis                : (theta,phi) = (90.0,90.0)

# Example
* WSe2 twisted bilayer system 

# Contributors
* Hyun-Jung Kim: Main developer (h.kim@fz-juelich.de, PGI-1/IAS-1, Forschungszentrum Jülich)
