F90=pgf90
F90FLAGS= -gopt -fastsse -fast -Mipa=fast,inline -Msmartalloc -Mfprelaxed -Mstack_arrays -O4 -Mvect=prefetch -Mprefetch -mcmodel=medium -Minfo -mp
MODINC="-I ./"
LDFLAGS= -gopt -fastsse -fast -Mipa=fast,inline -Msmartalloc -Mfprelaxed -Mstack_arrays -O4 -Mvect=prefetch -Mprefetch -mcmodel=medium -Minfo -mp