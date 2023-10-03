# This file is a part of DUDI, the Fortran-95 implementation 
# of the two-body model for dust dynamics
# Version 1.0.1
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
f95s = const.f95 define_types.f95 help.f95 distributions_fun.f95 inputdata.f95 dataoutmod.f95 gu.f95 twobody_fun.f95 integrator.f95

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

insitusampling.o : insitusampling.f95 $(objs) trajectory.o
	gfortran -c $< -fopenmp

enceladus_example.o : enceladus_example.f95 $(objs)
	gfortran -c $< -fopenmp

europa_example.o : europa_example.f95 $(objs)
	gfortran -c $< -fopenmp

io_example.o : io_example.f95 $(objs) image_construction.o
	gfortran -c $< -fopenmp

main_program.o : main_program.f95 $(objs)
	gfortran -c $< -fopenmp

image_construction.o : image_construction.f95
	gfortran -c -O3 $< -fopenmp

integrator.o : integrator.f95 const.o twobody_fun.o help.o define_types.o
	gfortran -c -O3 $< -fopenmp

twobody_fun.o : twobody_fun.f95 const.o help.o define_types.o distributions_fun.o
	gfortran -c -O3 $< -fopenmp

distributions_fun.o : distributions_fun.f95 const.o
	gfortran -c -O3 $< -fopenmp

dataoutmod.o : dataoutmod.f95 const.o
	gfortran -c $< -fopenmp

inputdata.o : inputdata.f95 const.o gu.o
	gfortran -c $< -fopenmp

gu.o : gu.f95 const.o distributions_fun.o
	gfortran -c $< -fopenmp

help.o : help.f95 const.o
	gfortran -c -O3 $< -fopenmp

define_types.o : define_types.f95 const.o
	gfortran -c $<

const.o : const.f95 
	gfortran -c $< -fopenmp

clean :
	rm *.mod *.o *dudi
