#
# Makefile for assl (Anomaly Detection Using Sparse Structure Learning)
#

prefix = /home/yamamori/opt
exec_prefix = /home/yamamori/opt
HOME   = /home/yamamori/opt
bindir = ${exec_prefix}/bin
LIBDIR = -L/home/yamamori/softwares/anom_ssl2/lib
exdir  = ./../ex/t4l_open-to-close/ene/

FC   = gfortran

FLAG = -O2

SRC = main.f90 graphical_lasso.f90 standardlize.f90

OBJS = main.o graphical_lasso.o standardlize.o

MOD_FILES = graphical_lasso.mod standardlize.mod

TARGET = assl2

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) \
	-llapack -lblas $(LIBDIR);

main.o: graphical_lasso.o standardlize.o
graphical_lasso.o:
standardlize.o:

.f90.o:
	$(FC) $(FLAG) -c $<;

.PHONY:clean install

clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET) 

distclean: clean
	rm -f Makefile config.status config.cache config.log

install:
	cp $(TARGET) $(bindir)

test:
	$(bindir)/$(TARGET) $(exdir)/ene-res_open.txt \
	$(exdir)/ene-res_close.txt > $(exdir)/ex.out
	diff $(exdir)/ex.out $(exdir)/out_e_res_res_open-close.txt
