CC = g++
CFLAGS = -std=c++17 -Wall

main: Main.o Geometry.o Initial_Boundary_Conditions.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o Norm.o
	$(CC) $(CFLAGS) -o Main Main.o Geometry.o Initial_Boundary_Conditions.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o Norm.o

Unit_Testing: Unit_Testing.o Geometry.o Initial_Boundary_Conditions.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o Norm.o
	$(CC) $(CFLAGS) -o Unit_Testing Unit_Testing.o Geometry.o Initial_Boundary_Conditions.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o Norm.o

Geometry.o:	Geometry.hpp Geometry.cpp
	$(CC) $(CFLAGS) -c Geometry.cpp
 
Initial_Boundary_Conditions.o: Initial_Boundary_Conditions.hpp Initial_Boundary_Conditions.cpp
	$(CC) $(CFLAGS) -c Initial_Boundary_Conditions.cpp 

Exact_Isentropic.o: Exact_Isentropic.hpp Exact_Isentropic.cpp
	$(CC) $(CFLAGS) -c Exact_Isentropic.cpp

SoundSpeed.o: SoundSpeed.hpp SoundSpeed.cpp
	$(CC) $(CFLAGS) -c SoundSpeed.cpp 

TimeStep.o: TimeStep.hpp TimeStep.cpp
	$(CC) $(CFLAGS) -c TimeStep.cpp 

VariableSwap.o: VariableSwap.hpp VariableSwap.cpp
	$(CC) $(CFLAGS) -c VariableSwap.cpp

Flux.o: Flux.hpp Flux.cpp
	$(CC) $(CFLAGS) -c Flux.cpp 

Artifical_Dissipation.o: Artificial_Dissipation.hpp Artificial_Dissipation.cpp
	$(CC) $(CFLAGS) -c Artificial_Dissipation.cpp
	
Source_Term.o: Source_Term.hpp Source_Term.cpp
	$(CC) $(CFLAGS) -c Source_Term.cpp

Norm.o: Norm.hpp Norm.cpp
	$(CC) $(CFLAGS) -c Norm.cpp

clean:
	rm -f core *.o Main Unit_Testing 
 

