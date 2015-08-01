#compiler=ifort -traceback -check
compiler=ifort -xHOST -O3 -static -parallel
compiler1=ifort -xHOST -O3 -static

#core = prmtr.o grid.o variable.o condition.o OUTPUT.o muscl.o lusgs.o euler_exp.o scheme.o viscous.o
core = prmtr.o OUTPUT.o muscl.o euler_exp.o scheme.o lusgs.o

all: a.out

prmtr.o: prmtr.f90
	$(compiler) -c $<
muscl.o: muscl.f90
	$(compiler) -c $<
#variable.o: variable.f90
#	$(compiler) -c $<
#grid.o: grid.f90 variable.o
#	$(compiler) -c $<
#condition.o: condition.f90 variable.o
#	$(compiler) -c $<
#lusgs.o: lusgs.f90 variable.o
#	$(compiler) -c $<
euler_exp.o: euler_exp.f90
	$(compiler) -c $<
lusgs.o: lusgs.f90
	$(compiler) -c $<
scheme.o: scheme.f90
	$(compiler) -c $<
OUTPUT.o: OUTPUT.f90
	$(compiler1) -c $<

a.out: main.f90 $(core)
	 $(compiler) $^
#	ifort 
#	ifort -fast -parallel 
#  ifort -traceback -idb -check all
clean:remove
	rm -rf *.mod a.out *.o data/*.d *.d *.bin *.png
remove:
tar:
	tar zcvf data.tar.gz data/*
clean_all:
	rm -rf *.mod a.out *.o data/*.d *.d *.bin *.gz
