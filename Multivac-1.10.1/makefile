##### Compiler, tools and options.

CC	=	g++
CFLAGS	=	-O2
INCPATH	=	-Iincludes
LINK	=	g++
LFLAGS	=	
LIBS	=

####### Targets.

TARGETS	= 	track track_extension M M_fast M_lax_friedrichs \
		track_eno2_eo merge_and_island piecewise_constant image

##### Implicit rules.

.SUFFIXES: .cpp .cxx .cc .C .c

.cpp.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

##### Build rules.

track:	track.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

track_extension:	track_extension.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

M:	M.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

M_fast:	M_fast.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

M_lax_friedrichs:	M_lax_friedrichs.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

track_eno2_eo:	track_eno2_eo.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

merge_and_island:	merge_and_island.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

piecewise_constant:	piecewise_constant.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

image:	image.o
	$(LINK) $(LFLAGS) -o $@ $@.o $(LIBS)

all: $(TARGETS)

clean:
	\rm -f $(TARGETS:=.o) $(TARGETS) core
