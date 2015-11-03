all:
	$(MAKE) --directory=libs DEBUG=no COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=no COMPILER=$(COMPILER)

debug:
	$(MAKE) --directory=libs DEBUG=yes COMPILER=$(COMPILER)
	$(MAKE) --directory=main DEBUG=yes COMPILER=$(COMPILER)

clean:
	-rm main/iscam.cpp
	-rm main/iscam.htp
	-rm -rf ../build/dist/

dclean:
	-rm main/iscam.cpp
	-rm main/iscam.htp
	-rm -rf ../build/debug/

rules:
	$(MAKE) rules --directory=libs COMPILER=$(COMPILER)
	$(MAKE) rules --directory=main COMPILER=$(COMPILER)
