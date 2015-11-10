## C. Grandin, Nov 2015. Tested for g++/clang on Linux and g++ on Windows (mingw).
##
## The makefiles for iscam assume you have the ADMB_HOME and ADMB_HOME_DEBUG
## environment variables set. They can be the same, but if you make debug
## using the non-debug (dist) build of ADMB, you will not be able to step into the
## ADMB source code during a debug session. You will however, still be able to
## debug the iscam.cpp file as it will be compiled with debug symbols.
## make clean will clean only the build/dist directory
## make dclean will clean only the build/debug directory

COMPILER := g++

.PHONY: clean dclean rules testcompiler

all: testcompiler
	$(MAKE) --directory=src DEBUG=no COMPILER=$(COMPILER)

debug: testcompiler
	$(MAKE) debug --directory=src COMPILER=$(COMPILER)

clean:
	$(MAKE) clean --directory=src

dclean:
	$(MAKE) dclean --directory=src

rules: testcompiler
	$(MAKE) rules --directory=src COMPILER=$(COMPILER)

testcompiler:
ifeq (, $(shell which $(COMPILER)))
	$(error "Using shell: $(SHELL), compiler $(COMPILER) not found in path.")
	exit 1
else
	@echo "Using shell: $(SHELL), compiler: $(COMPILER) found in path."
endif

