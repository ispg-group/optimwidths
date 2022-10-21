#!/bin/bash
ifort -g -c m_read.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o m_read.o
ifort -g -c m_misc.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o m_misc.o
ifort -g -c m_func.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o m_func.o
ifort -g -c m_optimizers.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o m_optimizers.o
ifort -g m_read.o m_misc.o m_func.o m_optimizers.o optimwidths.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o optimwidths.x
#./optimwidths.x
