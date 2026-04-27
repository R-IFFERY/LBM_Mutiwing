f90=mpif90
#f90 = ifort
#opt= -O1 -msse4.2  
#opt= -ipo -O1 -no-prec-div -fp-model fast=2 -xHost
opt= -O3
acc= -acc
Minfo= -Minfo
srs= main.f90 zmq_comm.f90 driver.f90 airfoil.f90 boundary_condition.f90 D2G9model.f90 ini_mesh.f90 D2Q9model.f90 mode.f90 parallel.f90 solution.f90 move_mesh.f90
OBJS=$(srs:.f90=.o)
LIB= -lzmq
#library = -L/opt/intel/impi/2018.1.163/intel64/lib

%.o:%.f90
	$(f90) $(acc) $(Minfo) $(opt) -gpu=deepcopy -c $<

default: $(OBJS)
	$(f90) $(acc) $(Minfo) $(opt) -gpu=deepcopy -cuda -cudalib=cublas -cudalib=cusolver -o exe $(OBJS) $(LIB)
clean:
	rm -f *.out *.o exe
