MODULE_TOPDIR = $(HOME)/dev/grass

PGM = i.atcorr2

# ── Sources ─────────────────────────────────────────────────────────────────
MOD_OBJS = main.o

# ── Compiler / linker options ────────────────────────────────────────────────
# Uses the system-installed libsixsv (package libsixsv-dev / libsixsv1).
# Headers: /usr/include/sixsv/   Library: /usr/lib/x86_64-linux-gnu/libsixsv.so
EXTRA_INC     = -I/usr/include/sixsv
EXTRA_CFLAGS  = -O3 -ffast-math -fopenmp -std=c11 \
                -Wall -Wextra -Wno-unused-parameter
EXTRA_LDFLAGS = -fopenmp
LIBES         = -lsixsv $(RASTERLIB) $(GISLIB) $(MATHLIB)

include $(MODULE_TOPDIR)/include/Make/Module.make

default: cmd

# ── Installation ─────────────────────────────────────────────────────────────
# Overrides Module.make's install; warning about "ignoring old recipe" is expected.
install: install_py
	$(INSTALL) $(ARCH_DISTDIR)/bin/$(PGM)$(EXE) $(INST_DIR)/bin/
	$(INSTALL_DATA) i.atcorr2.html $(INST_DIR)/docs/html/
	$(INSTALL_DATA) $(ARCH_DISTDIR)/docs/man/man1/$(PGM).1 \
	    $(INST_DIR)/docs/man/man1/ 2>/dev/null || true

# ── Python module installation ────────────────────────────────────────────────
# Installs the Python front-end, the sensors helper, and the SRF CSV data.
# sensors.py locates sensors_csv/ relative to __file__, so both must land
# in $(INST_DIR)/scripts/.
install_py:
	$(INSTALL) i.atcorr2.py $(INST_DIR)/scripts/$(PGM)
	$(INSTALL_DATA) sensors.py $(INST_DIR)/scripts/sensors.py
	$(INSTALL) -d $(INST_DIR)/scripts/sensors_csv
	$(INSTALL_DATA) sensors_csv/*.csv $(INST_DIR)/scripts/sensors_csv/
	$(INSTALL_DATA) sensors_csv/README $(INST_DIR)/scripts/sensors_csv/
	$(INSTALL_DATA) i.atcorr2.html $(INST_DIR)/docs/html/

.PHONY: install install_py
