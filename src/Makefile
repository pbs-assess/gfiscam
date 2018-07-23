## Do not run 'make' directly in this directory,
##  instead run it in the directory above.

dist:
	$(MAKE) --directory=libs DEBUG=no COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=no COMPILER=$(COMPILER)
	-rm -rf ../build/dist/examples
	-cp -rf ../examples ../build/dist/examples
	-cp ../build/dist/bin/iscam.* ../build/dist/examples/arrowtooth-flounder
	-cp ../build/dist/bin/iscam.* ../build/dist/examples/herringSOG
	-cp ../build/dist/bin/iscam.* ../build/dist/examples/pacific-cod
	-cp ../build/dist/bin/iscam.* ../build/dist/examples/petrale-sole

debug:
	$(MAKE) --directory=libs DEBUG=yes COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=yes COMPILER=$(COMPILER)
	-rm -rf ../build/debug/examples
	-cp -rf ../examples ../build/debug/examples
	-cp ../build/debug/bin/iscam.* ../build/debug/examples/arrowtooth-flounder
	-cp ../build/debug/bin/iscam.* ../build/debug/examples/herringSOG
	-cp ../build/debug/bin/iscam.* ../build/debug/examples/pacific-cod
	-cp ../build/debug/bin/iscam.* ../build/debug/examples/petrale-sole

clean-dist:
ifneq ("$(wildcard main/iscam.cpp)","")
	-rm main/iscam.cpp
endif
ifneq ("$(wildcard main/iscam.htp)","")
	-rm main/iscam.htp
endif
ifneq ("$(wildcard main/iscam.exe)","")
	-rm main/iscam.exe
endif
ifneq ("$(wildcard main/iscam)","")
	-rm main/iscam
endif
	-rm -rf ../build/dist/

clean-debug:
ifneq ("$(wildcard main/iscam.cpp)","")
	-rm main/iscam.cpp
endif
ifneq ("$(wildcard main/iscam.htp)","")
	-rm main/iscam.htp
endif
ifneq ("$(wildcard main/iscam.exe)","")
	-rm main/iscam.exe
endif
ifneq ("$(wildcard main/iscam)","")
	-rm main/iscam
endif
	-rm -rf ../build/debug/

rules:
	$(MAKE) rules --directory=libs COMPILER=$(COMPILER)
	$(MAKE) rules --directory=main COMPILER=$(COMPILER)
