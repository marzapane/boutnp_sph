#----------------------------------------
# Makefile per 
#----------------------------------------
# seriale
#----------------------------------------
COMP = GNU
OPTM = DEBUG
#----------------------------------------
# FC = f77
# FC = xlf
ifeq ($(COMP),GNU)
	FC = gfortran
	OP_COMP = -fdefault-real-8
	OP_COMP += -mcmodel=large -Mlarge_arrays
else
	FC = ifx
	OP_COMP = -r8
	OP_COMP = -mcmodel=large -heap-arrays
endif
#----------------------------------------
ifeq ($(OPTM),DEBUG)
# DEBUG OPTIONS
	ifeq ($(COMP),GNU)
		OP_COMP += -ggdb -O0 -Wall -fbacktrace -recursive -ffpe-trap=invalid,zero #,overflow,underflow
	else
		OP_COMP = -g -O0 -traceback -warn all,extra -check all -check pointers -init=snan -fpe0
	endif
else
# parallelo e ottimizzato
	ifeq ($(COMP),GNU)
		FC += -fopenmp
		OP_COMP += -O3
	else
		FC += -parallel -qopenmp
		OP_COMP += -fast
		OP_COMP += -xHost
		OP_COMP += -mfma
	endif
endif
#----------------------------------------
#FC = guidef77 -WGstats -WGstatic_library 
#FC = xlf_r  -qsmp=noauto  -qnosave
#FC = assuref77 -WAproject_name=paperino
#----------------------------------------
# IBM
#OP_COMP  = -qautodbl=dbl4 -O4 -qhot -qalign=4k -qtune=pwr4 -qarch=pwr4 -qnosave -qposition=append -q64 -b64
#OP_COMP2 = -qautodbl=dbl4 -O4 -qnoipa -qhot -qalign=4k -qtune=pwr4 -qarch=pwr4 -qnosave -qposition=append -q64 -b64
#OP_COMP3 = -qautodbl=dbl4 -O4 -qnohot  -qtune=pwr4 -qarch=pwr4 -qalign=4k -qnosave -q64 -b64
#OP_COMP = -qautodbl=dbl4 -O4 -qhot -qalign=4k -qtune=pwr3 -qarch=pwr3  -bmaxdata:1500000000
#----------------------------------------
# DIGITAL
#OP_COMP = -pg -fast -arch host -tune host -real_size 64  -align dcommon -automatic -convert big_endian
#OP_COMP2 = -pg -fast -arch host -tune host -real_size 64  -align dcommon -automatic -convert big_endian
#OP_COMP3 = -pg -fast -arch host -tune host -real_size 64  -align dcommon -automatic -convert big_endian
#OP_COMP = -r8 -fast -arch ev6 -tune ev6  -align dcommon  -automatic -convert big_endian
#OP_COMP2 = -r8 -fast -arch ev6 -tune ev6  -align dcommon  -automatic -convert big_endian
#----------------------------------------
# LINUX
# OP_COMP += -fdefault-real-8
# OP_COMP += -mcmodel=large
# OP_COMP += -O3
#----------------------------------------
BUINCLUDE = src/param.f
OBJS = obj/bounn.o obj/boutn.o obj/bouhd.o obj/pressnu.o obj/tsch.o obj/spect.o obj/ffts.o obj/fishpack.o

outnp_sph: directories $(OBJS)
	$(FC) $(OP_COMP) -o boutnp_sph $(OBJS)

obj/ffts.o : src/ffts.f $(BUINCLUDE)
	$(FC) -c -w $(OP_COMP) $< -o $@

obj/fishpack.o : src/fishpack.f $(BUINCLUDE)
	$(FC) -c -w $(OP_COMP) $< -o $@

obj/%.o : src/%.f $(BUINCLUDE)
	$(FC) -c $(OP_COMP) $< -o $@

directories:
	mkdir -p obj
	mkdir -p dati

clean:
	rm -rf obj/*.o
	rm boutnp_sph