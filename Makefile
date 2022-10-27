include make.vars

MAIN = optimwidths.o

MODS = m_read.o m_misc.o m_func.o m_optimizers.o

%.o : %.f90
	$(FC) $(FFLAG) -c $<

optimwidths.x: $(MAIN) $(MODS) 
	$(FC) $(FFLAG) -o $@ $(MAIN) $(MODS) $(LFLAG)

m_optimizers.o: m_func.o

optimwidths.o: $(MODS) 

install: optimwidths.x
