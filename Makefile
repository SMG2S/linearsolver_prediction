.PHONY: compile run clean

CC := /usr/bin/mpicxx
EX := /usr/bin/mpiexec

N := 6

OPTIONFLAG := -lpetsc ${SMG2SFLAG} ${PETSCFLAG} -lm
petsc.pc := ${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig/petsc.pc
PETSCFLAG  := $(shell pkg-config --libs-only-l ${petsc.pc})
SGM2SFLAG  := -lsmg2s

PETSC_ARCH_LIB := ${PETSC_DIR}/${PETSC_ARCH}/include

PETSC_HEADER   := ${PETSC_DIR}/include
SMG2S_HEADER   := ${SMG2S_DIR}/include/smg2s

SMG2S_PETSC_INTERFACE := ${SMG2S_DIR}/interface/PETSc
PETSC_ARCH_SO    := ${PETSC_DIR}/${PETSC_ARCH}/lib

EXEC := main.exe

all : compile run clean

compile : 
	@echo 'Compiling ...'
	@${CC} -c main.cpp SpectrumGenerator/generator.cpp -I${PETSC_ARCH_LIB} -I${PETSC_HEADER} -I${SMG2S_PETSC_INTERFACE} -L${PETSC_ARCH_SO} ${OPTIONFLAG}
	@${CC} main.o generator.o -o ${EXEC} -I${PETSC_ARCH_LIB} -I${PETSC_HEADER} -I${SMG2S_PETSC_INTERFACE} -L${PETSC_ARCH_SO} ${OPTIONFLAG}
run :
	@echo 'Launching ...'
	@${EX} -n ${N} ./main.exe

valgrind-run : compile
	@echo 'Launching ...'
	@valgrind ${EX} -n ${N} ./main.exe

check : compile valgrind-run clean

clean : 
	@echo 'Cleaning ...'
	@rm main.exe
	@rm *.o
