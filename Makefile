CC = g++
CFLAGS = -std=c++17 -Wall

main: Main.o Geometry.o Initial_Conditions.o Boundary_Conditions.o Isentropic_Flow.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o L2_Norm.o WriteFile.o
	$(CC) $(CFLAGS) -o Main Main.o Geometry.o Initial_Conditions.o Boundary_Conditions.o Isentropic_Flow.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o L2_Norm.o WriteFile.o

Unit_Testing: Unit_Testing.o Geometry.o Initial_Conditions.o Boundary_Conditions.o Isentropic_Flow.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o L2_Norm.o WriteFile.o
	$(CC) $(CFLAGS) -o Unit_Testing Unit_Testing.o Geometry.o Initial_Conditions.o Boundary_Conditions.o Isentropic_Flow.o Exact_Isentropic.o SoundSpeed.o TimeStep.o VariableSwap.o Flux.o Artificial_Dissipation.o Source_Term.o L2_Norm.o WriteFile.o

Geometry.o:	Geometry.hpp Geometry.cpp
	$(CC) $(CFLAGS) -c Geometry.cpp
 
Initial_Conditions.o: Initial_Conditions.hpp Initial_Conditions.cpp
	$(CC) $(CFLAGS) -c Initial_Conditions.cpp 

Boundary_Conditions.o: Boundary_Conditions.hpp Boundary_Conditions.cpp
	$(CC) $(CFLAGS) -c Boundary_Conditions.cpp 

Isentropic_Flow.o: Isentropic_Flow.hpp Isentropic_Flow.cpp
	$(CC) $(CFLAGS) -c Isentropic_Flow.cpp 

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

L2_Norm.o: L2_Norm.hpp L2_Norm.cpp
	$(CC) $(CFLAGS) -c L2_Norm.cpp
WriteFile.o: WriteFile.hpp WriteFile.cpp
	$(CC) $(CFLAGS) -c WriteFile.cpp

clean:
	rm -f core *.o Main Unit_Testing 
 

