## C. Grandin, Nov 2015. Tested for g++/clang on Linux and g++ on Windows (mingw).
##
## **Make sure to set up the ADMB_HOME and ADMB_HOME_DEBUG
##  environment variables for debug and dist builds of ADMB on your system

dist:
	$(MAKE) --directory=libs DEBUG=no COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=no COMPILER=$(COMPILER)
	-rm -rf ../build/dist/examples
	-cp -rf ../examples ../build/dist/examples
	-mv ../build/dist/examples/*.* ../build/dist/bin
	-cp ../build/dist/bin/iscam_arrowtooth.dat ../build/dist/bin/iscam.dat

debug:
	$(MAKE) --directory=libs DEBUG=yes COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=yes COMPILER=$(COMPILER)
	-rm -rf ../build/debug/examples
	-cp -rf ../examples ../build/debug/examples
	-mv ../build/debug/examples/*.dat ../build/debug/bin
	-cp ../build/debug/bin/iscam_arrowtooth.dat ../build/debug/bin/iscam.dat

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
