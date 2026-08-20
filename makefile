#----------------------------------------------------------
# Select Mode (default: serial)
# Usage:
#   make serial   -> Compiles Src_ser with nvfortran (CPU)
#   make cuda     -> Compiles Src_gpu with nvfortran (-cuda -gpu=cc80)
#----------------------------------------------------------
MODE ?= serial

ifeq ($(MODE),cuda)
    FC          := nvfortran
    SRC_DIR     := Src_gpu
    OBJ_DIR     := Build_gpu
    EXE         := wave2D_cuda.exe
    CUDA_FLAGS  := -cuda -gpu=cc90 -Minfo=accel
else
    FC          := nvfortran
    SRC_DIR     := Src_ser
    OBJ_DIR     := Build_ser
    EXE         := wave2D_serial.exe
    CUDA_FLAGS  :=
endif

#----------------------------------------------------------
# Compiler flags
#----------------------------------------------------------
CPPFLAGS    := -cpp
FFLAGS      := -O3 $(CPPFLAGS) $(CUDA_FLAGS)

# Read cppdefs.h, filter for #define, and convert to compiler flags
CPP_OPTIONS := $(shell grep '^\s*#define' cppdefs.h 2>/dev/null | \
                 grep -v 'CPPDEFS_H' | \
                 sed -e 's/^\s*#define\s\+\([^\s]\+\)\s*\(.*\)/-D\1=\2/' -e 's/=\s*$$//')

FFLAGS      += $(CPP_OPTIONS)

#----------------------------------------------------------
# Directories & Setup
#----------------------------------------------------------
MOD_DIR     := $(SRC_DIR)/Modules
NL_DIR      := $(SRC_DIR)
MyNC_DIR    := $(SRC_DIR)/MyNC
MOD_OUTPUT  := $(OBJ_DIR)/mod
PREPROC_DIR := preprocessed_$(MODE)

#----------------------------------------------------------
# NetCDF Configuration
#----------------------------------------------------------
NETCDF ?= /home/Alexander.Kurapov/netcdf-nvhpc
NF_CONFIG   := $(NETCDF)/bin/nf-config
NETCDF_INC  := $(shell $(NF_CONFIG) --fflags 2>/dev/null || echo "-I$(NETCDF)/include")
NETCDF_LIBS := -L$(NETCDF)/lib -lnetcdff -lnetcdf -lm

#----------------------------------------------------------
# Source files
#----------------------------------------------------------
MODULES := \
    $(MOD_DIR)/mync.F90
    # ... Add other module files here in dependency order

NL_SRC  := $(wildcard $(NL_DIR)/*.F90) $(wildcard $(MyNC_DIR)/*.F90) 

# Object files mapped to MODE-specific OBJ_DIR
MODULE_OBJS := $(patsubst $(SRC_DIR)/%.F90,$(OBJ_DIR)/%.o,$(MODULES))
NL_OBJS     := $(patsubst $(SRC_DIR)/%.F90,$(OBJ_DIR)/%.o,$(NL_SRC))
OBJ         := $(MODULE_OBJS) $(NL_OBJS)

#----------------------------------------------------------
# Target Rules
#----------------------------------------------------------
.PHONY: all serial cuda clean debug-preproc

all: serial

serial:
	@$(MAKE) MODE=serial build

cuda:
	@$(MAKE) MODE=cuda build

build: $(EXE)

#----------------------------------------------------------
# Compile modules and sources (.F90 -> .o)
#----------------------------------------------------------
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@) $(MOD_OUTPUT)
	$(FC) $(FFLAGS) $(NETCDF_INC) -module $(MOD_OUTPUT) -I$(MOD_OUTPUT) -I. -c $< -o $@

#----------------------------------------------------------
# Link executable
#----------------------------------------------------------
$(EXE): $(OBJ)
	$(FC) $(CUDA_FLAGS) -o $@ $(OBJ) $(NETCDF_LIBS) -Wl,-rpath,$(NETCDF)/lib

#----------------------------------------------------------
# Generate preprocessed files
#----------------------------------------------------------
PREPROCESSED := $(patsubst $(SRC_DIR)/%.F90,$(PREPROC_DIR)/%.f90,$(MODULES) $(NL_SRC))

$(PREPROC_DIR)/%.f90: $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@)
	$(FC) $(CPPFLAGS) $(CPP_OPTIONS) -E $< -o $@

debug-preproc: $(PREPROCESSED)

#----------------------------------------------------------
# Clean build files
#----------------------------------------------------------
clean:
	rm -rf Build_ser Build_gpu wave2D_serial.exe wave2D_cuda.exe preprocessed_serial preprocessed_cuda

