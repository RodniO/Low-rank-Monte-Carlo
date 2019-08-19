source = src
obj = ModIntVec.o ModVec.o ModMtrx.o ModSparse.o ModAppr.o ModRecon.o Main.o
incfiles = $(wildcard $(source)/incfiles/*.f90)

all:	intel

run:	intel_run

compiler_intel = ifort
opt_intel = -module ./obj_intel -O3 -xHost -no-wrap-margin
obj_intel = $(addprefix obj_intel/, $(obj))
exe_intel = Main.exe

intel_run:	intel
	@./$(exe_intel)
	
intel:	$(exe_intel)

$(exe_intel):	$(obj_intel)
	$(compiler_intel) $(opt_intel) -o $(exe_intel) $(obj_intel) -mkl
	
obj_intel/%.o:	$(source)/%.f90
	$(compiler_intel) $(opt_intel) -c $< -o $@
	
$(obj_intel): | obj_intel

obj_intel:
	mkdir obj_intel
	
obj_intel/Main.o:	$(source)/Main.f90	$(incfiles)
	$(compiler_intel) $(opt_intel) -c $< -o $@
	
compiler_gnu = gfortran
opt_gnu = -J./obj_gnu -I./obj_gnu -O3 -Wall -Wno-uninitialized -Wno-unused-function -Wno-unused-dummy-argument -fimplicit-none -pedantic
obj_gnu = $(addprefix obj_gnu/, $(obj))
exe_gnu = Main_gnu.exe

gnu_run:	gnu
	@./$(exe_gnu)

gnu:    $(exe_gnu)
	
obj_gnu/%.o:	$(source)/%.f90
	$(compiler_gnu) $(opt_gnu) -c $< -o $@
	
$(exe_gnu):	$(obj_gnu)
	$(compiler_gnu) $(opt_gnu) -o $(exe_gnu) $(obj_gnu) -llapack -lblas

$(obj_gnu): | obj_gnu

obj_gnu:
	mkdir obj_gnu
	
obj_gnu/Main.o:	$(source)/Main.f90	$(incfiles)
	$(compiler_gnu) $(opt_gnu) -c $< -o $@
	
help:
	@ echo ""
	@ echo "Use 'make gnu' for gfortran (requires BLAS and LAPACK)"
	@ echo "then run '$(exe_gnu)'"
	@ echo "Or use 'make gnu_run' to compile and run"
	@ echo ""
	@ echo "Use 'make'/'make intel' for ifort (requires mkl)"
	@ echo "then run '$(exe_intel)'"
	@ echo "Or use 'make run'/'make intel_run' to compile and run"
	@ echo ""

clear:
	rm -rf obj_gnu
	rm -rf obj_intel
	rm -f *.exe

#lib:    compile
#	@ ar cr ./obj/*.o
