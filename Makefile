#
# Makefile for assl
#

TARGET = assl

FC   = gfortran

FLAG = -O2

SRC = main.f90 graphical_Lasso.f90 standardlize.f90

OBJS = main.o graphical_Lasso.o standardlize.o

MOD_FILES = graphical_Lasso.mod standardlize.mod

.SUFFIXES: .f90

All:$(TARGET)

${TARGET}:${OBJS}
	$(FC) $(FLAG) -o $@ $(OBJS) \
	-llapack -lblas -L${HOME}/lib;

main.o: graphical_Lasso.o standardlize.o
graphical_Lasso.o:
standardlize.o:

.f90.o:
	$(FC) -c $<;

.PHONY:clean install
clean:
	rm $(OBJS) $(MOD_FILES) $(TARGET) 
install:
	cp $(TARGET) ../bin
