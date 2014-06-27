include Makefile.cf.linux

.PHONY: all prepare build clean prepare-install install

all:
	$(call SEQ_MAKE, $(TARGETS) )

prepare:
	$(call RECURSIVE_MAKE, $(SUBDIRS), prepare)

build:
	$(call RECURSIVE_MAKE, $(SUBDIRS), all)

clean:
	$(call RECURSIVE_MAKE, $(SUBDIRS), clean)

install: prepare-install
	$(call RECURSIVE_MAKE, $(SUBDIRS), install)

prepare-install:
	mkdir -p $(INSTALL_INCLUDE_DIR)
	mkdir -p $(INSTALL_BIN_DIR)

