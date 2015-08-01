compiler = ifort -fast -parallel
directory=../
core = prmtr.o variable.o grid.o
all: a.out

prmtr.o:prmtr.f90
	$(compiler) -c $<
variable.o:variable.f90
	$(compiler) -c $<
grid.o:grid.f90
	$(compiler) -c $<

a.out: postprocess.f90 $(core)
	$(compiler) $^

clean:remove
	rm -rf *.mod a.out *.o *.vtk result/*.vtk *.bin
remove:
tar:
	tar zcvf result.tar.gz result/*
clean_all:
	rm -rf *.mod a.out *.o *.vtk result/*.vtk *.bin *.gz
