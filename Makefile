COMPILER:=g++
all:
	$(MAKE) --directory=src DEBUG=no COMPILER=$(COMPILER)

debug:
	$(MAKE) debug --directory=src COMPILER=$(COMPILER)

clean:
	$(MAKE) clean --directory=src

dclean:
	$(MAKE) dclean --directory=src

rules:
	$(MAKE) rules --directory=src COMPILER=$(COMPILER)
