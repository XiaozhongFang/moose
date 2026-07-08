# kpp/build/Makefile — Build KPP shared libraries for box-mode chemistry
# Usage: make -f kpp/build/Makefile MECH=path/to/mechanism.kpp
SHELL := /bin/bash
KPP_BIN ?= $(shell command -v kpp 2>/dev/null || echo "kpp")

MECH_NAME := $(basename $(notdir $(MECH)))
MECH_DIR := $(dir $(realpath $(MECH)))
OUT_DIR ?= $(MECH_DIR)kpp_build_$(MECH_NAME)
LIB := $(OUT_DIR)/libkpp_$(MECH_NAME).so
ADAPTER_SRC := $(realpath $(dir $(lastword $(MAKEFILE_LIST)))/../adapter/kpp_adapter.c)
KPP_INTEGRATOR ?= $(shell awk 'tolower($$1)=="#integrator"{print tolower($$2); exit}' $(MECH) 2>/dev/null)
KPP_INTEGRATOR := $(if $(KPP_INTEGRATOR),$(KPP_INTEGRATOR),rosenbrock)
KPP_ADAPTER_CPPFLAGS :=
ifneq ($(KPP_INTEGRATOR),rosenbrock)
  KPP_ADAPTER_CPPFLAGS += -DKPP_USE_GENERATED_INTEGRATE
endif

.PHONY: all clean
all: $(LIB)

$(LIB): $(ADAPTER_SRC) $(realpath $(MECH))
	@echo "=== KPP Build: $(MECH_NAME) ==="
	rm -rf $(OUT_DIR) && mkdir -p $(OUT_DIR)
	set -e; \
	model=$$(awk 'tolower($$1)=="#model"{print $$2; exit}' $(realpath $(MECH))); \
	files="$(realpath $(MECH))"; \
	if [ -n "$$model" ]; then files="$$files $(MECH_DIR)$$model.spc $(MECH_DIR)$$model.eqn"; fi; \
	for src in $$files; do \
		[ -f "$$src" ] || continue; \
		cp "$$src" $(OUT_DIR)/; \
		src_dir=$$(dirname "$$src")/; \
		for incl in $$(awk 'tolower($$1)=="#include"{print $$2}' "$$src"); do \
			incl=$${incl%\"}; incl=$${incl#\"}; \
			if [ -f "$$src_dir$$incl" ]; then cp "$$src_dir$$incl" $(OUT_DIR)/; \
			elif [ -f "$(MECH_DIR)$$incl" ]; then cp "$(MECH_DIR)$$incl" $(OUT_DIR)/; fi; \
		done; \
	done
	cd $(OUT_DIR) && $(KPP_BIN) $(MECH_NAME).kpp . 2>&1 | tail -1
	for f in $(OUT_DIR)/._*; do \
		base=$${f#$(OUT_DIR)/._}; mv $$f $(OUT_DIR)/$(MECH_NAME)_$$base; \
	done
	sed -i 's/\._/$(MECH_NAME)_/g' $(OUT_DIR)/$(MECH_NAME)_*.c $(OUT_DIR)/$(MECH_NAME)_*.h 2>/dev/null || true
	sed -i 's/\([0-9]\)\.\([0-9]*\)e\([^0-9+-]\)/\1.\2e0\3/g' $(OUT_DIR)/$(MECH_NAME)_*.c 2>/dev/null || true
	printf '\nint kpp_get_nspec(void) { return NSPEC; }\n' >> $(OUT_DIR)/$(MECH_NAME)_Main.c
	printf '\nint kpp_get_nvar(void) { return NVAR; }\n' >> $(OUT_DIR)/$(MECH_NAME)_Main.c
	printf '\nint kpp_get_nreact(void) { return NREACT; }\n' >> $(OUT_DIR)/$(MECH_NAME)_Main.c
	set -e; for src in $(OUT_DIR)/$(MECH_NAME)_*.c; do \
		base=$$(basename $$src); \
		case $$base in *_mex_*) continue;; esac; \
		echo "  CC $$base"; \
		gcc -O2 -fPIC -I$(OUT_DIR) -c $$src -o $${src%.c}.o; \
	done
	echo "  CC kpp_adapter.c"
	gcc -O2 -fPIC -I$(OUT_DIR) $(KPP_ADAPTER_CPPFLAGS) -c $(ADAPTER_SRC) -o $(OUT_DIR)/kpp_adapter.o
	echo "  LINK libkpp_$(MECH_NAME).so"
	gcc -shared -o $(LIB) $(OUT_DIR)/$(MECH_NAME)_*.o $(OUT_DIR)/kpp_adapter.o -lm
	@echo "=== Done: $(LIB) ==="
	@ls -lh $(LIB)
	@echo "  META $(OUT_DIR)/kpp_$(MECH_NAME).json"
	@NSPEC=$$(grep -m1 '^#define NSPEC ' $(OUT_DIR)/$(MECH_NAME)_Parameters.h 2>/dev/null | awk '{print $$3}'); \
	NVAR=$$(grep -m1 '^#define NVAR ' $(OUT_DIR)/$(MECH_NAME)_Parameters.h 2>/dev/null | awk '{print $$3}'); \
	NREACT=$$(grep -m1 '^#define NREACT ' $(OUT_DIR)/$(MECH_NAME)_Parameters.h 2>/dev/null | awk '{print $$3}'); \
	echo "{\"mechanism\":\"$(MECH_NAME)\",\"n_species\":$${NSPEC:-0},\"n_reactions\":$${NREACT:-0},\"n_var\":$${NVAR:-0}}" \
	  > $(OUT_DIR)/kpp_$(MECH_NAME).json

clean:
	rm -rf $(OUT_DIR)
