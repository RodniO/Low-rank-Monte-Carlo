.SUFFIXES: .o .exe .f90
.f90.o:
	$(compiler) $(opt) -c $<

compiler = ifort
opt = -module ./obj_ifort -O3 -xHost -no-wrap-margin
#opt = -module ./obj_ifort -g -traceback -check all -xHost -qopenmp
obj = ModVec.o ModMtrx.o Main.o
exe = Main.exe

all:	exe

compile:	$(obj)

help:
	@ echo "No help for now."

lib:    compile
	@ ar cr ./obj/*.o

exe:    compile
	$(compiler) $(opt) -o $(exe) $(obj) -mkl

 #clean: @ rm ./obj/*.o
clean:
	rm -f $(obj)
	rm -f $(exe)
