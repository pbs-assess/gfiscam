all:
	$(MAKE) --directory=libs DEBUG=no COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=no COMPILER=$(COMPILER)

debug:
	$(MAKE) --directory=libs DEBUG=yes COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=yes COMPILER=$(COMPILER)

clean:
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

dclean:
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
