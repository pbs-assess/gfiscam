## Build iscam executables
##
## To build optimized, distribution version:
## Add a ADMB_HOME environment variable pointing to your ADMB install directory
## Command to build: make dist
## Command to clean: make clean-dist
##
## To build non-optimized, debug version (for use with GDB debugger):
## Add a ADMB_HOME_DEBUG environment variable pointing to your ADMB install
## directory. This can be the same as the ADMB_HOME directory.
## Command to build: make debug
## Command to clean: make clean-debug
##
## You can build both by running 'make' without arguments.
## You can clean both by running 'make clean'. This will also delete the
##  'build' directory and all of its contents.

COMPILER := g++

.PHONY: clean clean-debug rules testcompiler

all: debug dist

dist: testcompiler
	$(MAKE) dist --directory=src COMPILER=$(COMPILER)
	@echo
	@echo iSCAM distribution version built successfully.
	@echo

debug: testcompiler
	$(MAKE) debug --directory=src COMPILER=$(COMPILER)
	@echo
	@echo iSCAM debug version built successfully.
	@echo

clean-dist:
	$(MAKE) clean-dist --directory=src

clean-debug:
	$(MAKE) clean-debug --directory=src

clean: clean-dist clean-debug
	-rm -rf build

rules: testcompiler
	$(MAKE) rules --directory=src COMPILER=$(COMPILER)

testcompiler:
ifeq (, $(shell which $(COMPILER)))
	$(error "Using shell: $(SHELL), compiler $(COMPILER) not found in path.")
	exit 1
else
	@echo "Using shell: $(SHELL), compiler: $(COMPILER) found in path."
endif
