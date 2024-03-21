# mpiramS
fork of mpiram from [Brian Dushaw](https://staff.washington.edu/dushaw/AcousticsCode/RamFortranCode.html)
multi-processing version of RAM

## Steps to compile, run mpiram, and plot the output in a nutshell:

### LAPACK:
If you don't have ready access to lapack libraries:
(a) cd to the lapack directory 
   1. Edit the make file for your compiler.
   2. Compile the lapack libraries there, i.e., just "make".
      This should make a library smlapack.a (sm for "small"...)

### COMPILE:
(b) Back in the main directory, edit the Makefile for your compiler.
   1. Makefile is for single processor. 
   2. Makefile.mpi is for MPI parallel version.

(c) Compile:  Just "make" for single processor version (s_mpiram), 
    or "make -f Makefile.mpi" for MPI version (mpiram).

### RUN:
(d) The program opens and reads the input file "in.pe".  You may edit 
    it to suit your taste, but the read is formatted so be careful
    to respect the format.  Example sound speed and bathymetry
    files are included (test.ssp, test.bth) for convenience.

(e) Execute single processor version by just "s_mpiram".

(f) Execute the MPI version:
   1.  Start up the mpd daemons on all nodes (mpich in mind here).
       (See the util/startemup script)
   2.  Edit the "runit" script to set the number of desired processors.
   3.  Then just run it:  "runit".

(g) Don't forget that if you modify the mpi executable, then that 
    executable has to be copied/updated to all the remote machines.  
    This isn't an issue if the system is mounted NFS as aogreef is.

### PLOT RESULTS:
(h) The program (single or parallel) writes out the direct access
    binary file "psif.dat" which contains all the information needed
    for plotting the result.

(i) Start up matlab, and execute "plotram" to see that result.
    * If that fails, then your compiler likely gives a different result
    for what the record length means than my compiler.  See the top of 
    the plotram.m script for notes.  Multiply the recl number by 4?


B.D.  May 2008
