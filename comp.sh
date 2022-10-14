#!/bin/bash
ifort -c m_func.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o m_func.o
ifort -c m_optimizers.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o m_optimizers.o
ifort  m_func.o m_optimizers.o optimwidths.f90 -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -o optimwidths.x
./optimwidths.x
