.EXPORT_ALL_VARIABLES:
.SUFFIXES:	.gz

# the operating system independent macros

# the operating system dependent macros
CC	= cc
CFLAGS	= -I$$HOME/work/software/cfitsio -I$$HOME/work/software/fftw-2.1.5/rfftw -g -I$$HOME/work/software/fftw-2.1.5/fftw -g
LFLAGS	= -L$$HOME/work/software/cfitsio -L$(FFTW2)/rfftw/.libs -L$(FFTW2)/fftw/.libs 

# more operating system independent macros
LD	= $(CC)

# other stufqf
TAR	= makefile ppmtompeg *.[ch] *.sm
CLEAN	= core gmon.out *.o *~ \\\#*
TIDY	= .smhist test* $(PROGS)


# program names and corresponding object files
PROGS	= beyond combine reversi raw2fits
MAIN	= beyond
MAINOBJS	= main.o org.o field.o minkowski.o
MAINLIBS	= -lcfitsio -lm -lrfftw -lfftw 

default:	$(PROGS)

$(MAIN):	$(MAINOBJS)	
		$(LD) $^ $(LFLAGS) $(MAINLIBS) -o $@

test   :	$(MAIN)
		./$(MAIN) -x64 -y64 -z64 -i1 -s0 -c0 -t-3.6 -T3.6 -l1 -n24 -ddatatest.fits -mmask_grid64.fits -omftest.dat >> /dev/null

.c.o:	
	$(CC) -c $(CFLAGS) $<
.gz:	
	gunzip $<

# dependence of code on headers
main.o:		main.h org.h field.h minkowski.h
org.o:		main.h
field.o:	main.h org.h
minkowski.o:	main.h org.h

# small utilities built from just a single file each
combineLIB	= -lm
raw2fitsLIB	= -lcfitsio -lm
rversiLIB	= 
%:	%.c
	$(CC) $(CFLAGS) $@.c $(LFLAGS) $($(@)LIB) -o $@

# cleanup
.PHONY:	tar clean tidy
tar:	Source.tar.gz
Source.tar.gz:	$(TAR)
	tar cvf - $(TAR) | gzip --best > Source.tar.gz
	cp Source.tar.gz Source-`date -u +%d-%B-%C%y`.tar.gz
clean:	
	rm -f $(CLEAN)
tidy:	clean
	rm -f $(TIDY)
