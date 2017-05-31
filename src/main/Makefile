## Do not run 'make' directly in this directory,
##  instead run it two directories above.

ifeq ($(DEBUG),yes)
ADMB          := $(ADMB_HOME_DEBUG)
DESTDIR       := ../../build/debug/
OBJDIR        := ../../build/debug/objects/
BINDIR        := ../../build/debug/bin/
LINKERFLAGS   :=
LINKERLIBS    :=  $(ADMB)/lib/libadmb.a $(ADMB)/lib/libadmb-contrib.a
COMPILERFLAGS := -g -D__GNUDOS__ -Dlinux -DUSE_LAPLACE -std=c++11 -I. -I$(ADMB)/include -I$(ADMB)/contrib/include
LIBOBJS       := $(wildcard ../../build/debug/objects/*.o)
else
ADMB          := $(ADMB_HOME)
DESTDIR       := ../../build/dist/
OBJDIR        := ../../build/dist/objects/
BINDIR        := ../../build/dist/bin/
LINKERFLAGS   :=
LINKERLIBS    := -ladmb-contrib
COMPILERFLAGS := -c -O3 -D_FILE_OFFSET_BITS=64 -Wall -DSAFE_ALL -D__GNUDOS__ -Dlinux -DUSE_LAPLACE  -std=c++11 -I. -I$(ADMB)/include -I$(ADMB)/contrib/include
LIBOBJS       := $(wildcard ../../build/dist/objects/*.o)
endif
# LIBOBJS may contain iscam.o, if so remove it from the list so linker call does not have two iscam.o's.
LIBOBJS := $(filter-out ../../build/dist/objects/iscam.o, $(LIBOBJS))

# Must come after the if statement so debug or dist dirs are correctly prepended
EXEC := iscam
SRCS := $(EXEC).cpp
OBJS := $(EXEC).o

LINKER := g++
OBJECTS := $(addprefix $(OBJDIR), $(OBJS))

all: $(EXEC)
	@echo
	@echo iSCAM main code built successfully, DEBUG = $(DEBUG)
	@echo

$(EXEC): $(OBJECTS)
	$(LINKER) $(LINKERFLAGS) -o $@ $^ $(LIBOBJS) -L$(ADMB)/lib $(LINKERLIBS)
	@if [ ! -d "$(BINDIR)" ]; then mkdir -p $(BINDIR); fi
	-cp $(EXEC) $(BINDIR)
	-mv $(EXEC).exe $(BINDIR)
	-mv $(EXEC).cpp $(BINDIR)
	-mv $(EXEC).htp $(BINDIR)

$(EXEC).cpp:
	$(ADMB)/bin/tpl2cpp $(EXEC)

$(OBJDIR)%.o: %.cpp
	@if [ ! -d "$(DESTDIR)" ]; then mkdir -p $(DESTDIR); fi
	@if [ ! -d "$(OBJDIR)" ]; then mkdir -p $(OBJDIR); fi
	$(COMPILER) $(COMPILERFLAGS) -c $< -o $@

# Show variable values and compiler information
rules:
	@echo EXEC = $(EXEC)
	@echo DESTDIR = $(DESTDIR)
	@echo OBJDIR = $(OBJDIR)
	@echo OBJS = $(OBJS)
	@echo OBJECTS = $(OBJECTS)
	@echo Library objects = $(LIBOBJS)
	@echo SRCS = $(SRCS)
	@echo ADMB = $(ADMB)
	@echo MAKE = $(MAKE)
	@echo Debug = $(DEBUG)
	@echo Compiler version = $(COMPILER) --version
	@echo Compiler flags = $(COMPILERFLAGS)
	@echo Linker = $(LINKER)
	@echo Linker flags = $(LINKERFLAGS)
	@echo Library object files = $(LIBOBJS)
