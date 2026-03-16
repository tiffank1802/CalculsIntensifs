OPTI=-O
#OPTI=-g -O0
OPT=O_c++

CXX=$(MPI_BIN)/mpicxx $(OPTI)

MPI_INC= -I$(MPI_INCLUDE)
MPI_LIBRARY= -L$(MPI_LIB)
MPI_LIBS= -lmpi

EXE= somme

C++SOURCES= principal somme

SOURCESFILES= $(C++SOURCES:=.cc)
OBJECTSFILES= $(C++SOURCES:=.o)

default : all
all :
	@echo 'Compilation de ' : 
	$(CXX) $(MPI_INC) -c $(SOURCESFILES)

	@echo 'Edition de lien'
	$(CXX) -o $(EXE) $(OBJECTSFILES) $(MPI_INC) $(MPI_LIBRARY) $(MPI_LIBS) -ldl

clean :
	rm -rf *.o
	rm -rf somme

