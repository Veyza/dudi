# This file is a part of DUDI, the Fortran-90 implementation 
# of the two-body model for dust dynamics
# Version 1.1.0
# This is free software. You can use and redistribute it 
# under the terms of the GNU General Public License (http://www.gnu.org/licenses/)
# If you do, please cite the following paper
# Anastasiia Ershova and Jürgen Schmidt, 
# Two-body model for the spatial distribution of dust ejected from
# an atmosphereless body, 2021, A&A, 650, A186 

# Author: Anastasiia Ershova
# E-mail: vveyzaa@gmail.com
objs = const.o define_types.o help.o distributions_fun.o inputdata.o dataoutmod.o gu.o twobody_fun.o integrator.o
mods = const.mod define_types.mod help.mod distributions_fun.mod inputdata.mod dataoutmod.mod gu.mod twobody_fun.mod integrator.mod
f95s = const.f90 define_types.f90 help.f90 distributions_fun.f90 inputdata.f90 dataoutmod.f90 gu.f90 twobody_fun.f90 integrator.f90

main : main_program.o $(objs)
	gfortran -fopenmp $(objs) main_program.o -o dudi 

enceladus_comp : enceladus_example.o $(objs)
	gfortran -fopenmp $(objs) enceladus_example.o -o enceladus_dudi

europa_comp : europa_example.o $(objs)
	gfortran -fopenmp $(objs) europa_example.o -o europa_dudi

io_comp : io_example.o $(objs) image_construction.o
	gfortran -fopenmp $(objs) image_construction.o io_example.o -o io_dudi

enceladus : enceladus_comp
	./enceladus_dudi
	python3 e2plot.py 

europa : europa_comp
	./europa_dudi
	python3 deposition.py

io : io_comp
	./io_dudi
	python3 volcano_image.py

insitusampling.o : insitusampling.f90 $(objs) trajectory.o
	gfortran -c $< -fopenmp

enceladus_example.o : enceladus_example.f90 $(objs)
	gfortran -c $< -fopenmp

europa_example.o : europa_example.f90 $(objs)
	gfortran -c $< -fopenmp

io_example.o : io_example.f90 $(objs) image_construction.o
	gfortran -c $< -fopenmp

main_program.o : main_program.f90 $(objs)
	gfortran -c $< -fopenmp

image_construction.o : image_construction.f90
	gfortran -c -O3 $< -fopenmp

integrator.o : integrator.f90 const.o twobody_fun.o help.o define_types.o
	gfortran -c -O3 $< -fopenmp

twobody_fun.o : twobody_fun.f90 const.o help.o define_types.o distributions_fun.o
	gfortran -c -O3 $< -fopenmp

distributions_fun.o : distributions_fun.f90 const.o
	gfortran -c -O3 $< -fopenmp

dataoutmod.o : dataoutmod.f90 const.o
	gfortran -c $< -fopenmp

inputdata.o : inputdata.f90 const.o gu.o
	gfortran -c $< -fopenmp

gu.o : gu.f90 const.o distributions_fun.o
	gfortran -c $< -fopenmp

help.o : help.f90 const.o
	gfortran -c -O3 $< -fopenmp

define_types.o : define_types.f90 const.o
	gfortran -c $<

const.o : const.f90 
	gfortran -c $< -fopenmp

clean :
	rm *.mod *.o *dudi
