#----------------------------------------------------------
# Compiler and flags
#----------------------------------------------------------
FC       = nvfortran
CPPFLAGS = -cpp
FFLAGS   = -O2 $(CPPFLAGS)

# --- ADDED SECTION to read cppdefs.h ---

# Read cppdefs.h, filter for #define, and convert to compiler flags (e.g., -DOPTION=1)
# NOTE: Change 'CPPDEFS_H' below if your header guard has a different name.
CPP_OPTIONS := $(shell grep '^\s*#define' cppdefs.h | \
                 grep -v 'CPPDEFS_H' | \
                 sed -e 's/^\s*#define\s\+\([^\s]\+\)\s*\(.*\)/-D\1=\2/' -e 's/=\s*$$//')

# Add the extracted options to your Fortran compiler flags
FFLAGS += $(CPP_OPTIONS)

# --- END of ADDED SECTION ---

#----------------------------------------------------------
# Directories
#----------------------------------------------------------
SRC_DIR      := Src
MOD_DIR      := $(SRC_DIR)/Modules
NL_DIR       := $(SRC_DIR)
MyNC_DIR     := $(SRC_DIR)/MyNC
OBJ_DIR      := Build
MOD_OUTPUT   := $(OBJ_DIR)/mod
PREPROC_DIR  := preprocessed

#----------------------------------------------------------
# NetCDF
#----------------------------------------------------------
NETCDF ?= /home/Alexander.Kurapov/netcdf-nvhpc
NF_CONFIG  := $(NETCDF)/bin/nf-config
NETCDF_INC := $(shell $(NF_CONFIG) --fflags)
NETCDF_LIBS := -L$(NETCDF)/lib -lnetcdff -lnetcdf -lm

#----------------------------------------------------------
# Source files
#----------------------------------------------------------
# List module sources in proper dependency order
MODULES := \
    $(MOD_DIR)/mync.F90      
    #... Add other module files here in dependency order

# Nonlinear source files
NL_SRC := $(wildcard $(NL_DIR)/*.F90) $(wildcard $(MyNC_DIR)/*.F90) 

# Object files
MODULE_OBJS := $(patsubst $(SRC_DIR)/%.F90,$(OBJ_DIR)/%.o,$(MODULES))
NL_OBJS     := $(patsubst $(SRC_DIR)/%.F90,$(OBJ_DIR)/%.o,$(NL_SRC))
OBJ         := $(MODULE_OBJS) $(NL_OBJS)

# Executable
EXE = wave2D.exe

#----------------------------------------------------------
# Default target
#----------------------------------------------------------
all: $(EXE)

#----------------------------------------------------------
# Compile modules and sources
#----------------------------------------------------------
# Pattern rule for .F90 -> .o
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@) $(MOD_OUTPUT)
	$(FC) $(FFLAGS) $(NETCDF_INC) -module $(MOD_OUTPUT) -I$(MOD_OUTPUT) -c $< -o $@

#----------------------------------------------------------
# Link executable
#----------------------------------------------------------
$(EXE): $(OBJ)
	$(FC) -o $@ $(OBJ) $(NETCDF_LIBS) -Wl,-rpath,$(NETCDF)/lib

#----------------------------------------------------------
# Generate preprocessed files
#----------------------------------------------------------
PREPROCESSED := $(patsubst $(SRC_DIR)/%.F90,$(PREPROC_DIR)/%.f90,$(MODULES) $(NL_SRC))

$(PREPROC_DIR)/%.f90: $(SRC_DIR)/%.F90
	@mkdir -p $(dir $@)
	$(FC) $(CPPFLAGS) -E $< -o $@

debug-preproc: $(PREPROCESSED)

#----------------------------------------------------------
# Clean build files
#----------------------------------------------------------
clean:
	rm -rf $(OBJ_DIR) $(EXE) $(PREPROC_DIR)

.PHONY: all clean debug-preproc

