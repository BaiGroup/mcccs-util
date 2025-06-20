SUBDIRS = Fortran C

#noinst_SCRIPTS = contrib tests
#noinst_MANS = doc

# most fortran compilers can directly preprocess c-like directives: use
# 	$(MPIF90) $(F90FLAGS) -c $<
# if explicit preprocessing by the C preprocessor is needed, use:
# 	$(CPP) $(CPPFLAGS) $< -o $*.F90 
#	$(MPIF90) $(F90FLAGS) -c $*.F90 -o $*.o
# remember the tabulator in the first column !!!

.f.o:
	@f90rule@

mostlyclean-generic:
	-rm -f *.mod

# DFLAGS  = precompilation options (possible arguments to -D and -U)
#           used by the C compiler and preprocessor
# FDFLAGS = as DFLAGS, for the f90 compiler
# With the exception of IBM xlf, FDFLAGS = $(DFLAGS)
# For IBM xlf, FDFLAGS is the same as DFLAGS with separating commas 

DFLAGS         = @dflags@
FDFLAGS        = @fdflags@

# IFLAGS = how to locate directories where files to be included are
# In most cases, IFLAGS = -I../include

IFLAGS         = @iflags@

# MOD_FLAGS = flag used by f90 compiler to locate modules
# Each Makefile defines the list of needed modules in MODFLAGS

MOD_FLAG      = @imod@

# Compilers: fortran-90, fortran-77, C
# If a parallel compilation is desired, MPIF90 should be a fortran-90 
# compiler that produces executables for parallel execution using MPI
# (such as for instance mpif90, mpf90, mpxlf90,...);
# otherwise, an ordinary fortran-90 compiler (f90, g95, xlf90, ifort,...)
# If you have a parallel machine but no suitable candidate for MPIF90,
# try to specify the directory containing "mpif.h" in IFLAGS
# and to specify the location of MPI libraries in MPI_LIBS

MPIF90         = @mpif90@
F90            = $(MPIF90)
CC             = @cc@
F77            = $(MPIF90)

# C preprocessor and preprocessing flags - for explicit preprocessing, 
# if needed (see the compilation rules above)
# preprocessing flags must include DFLAGS and IFLAGS

CPP            = @cpp@
CPPFLAGS       = @cppflags@ $(DFLAGS) $(IFLAGS)

# compiler flags: C, F90, F77
# F90 flags must include MODFLAGS, IFLAGS, and FDFLAGS with appropriate syntax

CFLAGS         = @cflags@ $(DFLAGS) $(IFLAGS)
F90FLAGS       = @f90flags@ @pre_fdflags@$(FDFLAGS) $(IFLAGS) $(MODFLAGS)
FFLAGS         = @fflags@

# compiler flags without optimization for fortran-77
# the latter is NEEDED to properly compile dlamch.f, used by lapack

FFLAGS_NOOPT   = @fflags_noopt@

# Linker, linker-specific flags (if any)
# Typically LD coincides with F90 or MPIF90, LD_LIBS is empty

LD             = @ld@
LDFLAGS        = @ldflags@
LD_LIBS        = @ld_libs@

# For parallel execution, the correct path to MPI libraries must
# be specified in MPI_LIBS (except for IBM if you use mpxlf)

MPI_LIBS       = @mpi_libs@

# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv

AR             = @ar@
ARFLAGS        = @arflags@

# ranlib command. If ranlib is not needed (it isn't in most cases) use
# RANLIB = echo

RANLIB         = @ranlib@

LIBS           = $(MPI_LIBS) $(LD_LIBS)
