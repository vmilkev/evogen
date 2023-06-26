APPNAME := evogen

CXX := icpc
SRCDIR := src
BDIR_PYLIB := release
BDIR_CLIB := release2
BDIR_APP := release3
TESTDIR := tests
BDIR_TST := $(TESTDIR)/release
TSRCDIR := $(TESTDIR)/src

UNAME_S := $(shell uname -s)
PWDADDR := $(shell pwd)
PYBINDINCL := $(shell python3 -m pybind11 --includes)
PYSUFF := $(shell python3-config --extension-suffix)
VERSUFF := .so

ifeq ($(UNAME_S),Linux)
	CFLAGS_TEST := -O3 -qopenmp -I$(SRCDIR) -DMKL_ILP64 -I${MKLROOT}/include $(PYBINDINCL) -std=c++14
	CFLAGS_APP := -O3 -qopenmp -I$(SRCDIR) -DMKL_ILP64 -I${MKLROOT}/include -std=c++14
	CFLAGS_PYLIB := -O3 -qopenmp -DPYBIND -DMKL_ILP64 -I${MKLROOT}/include $(PYBINDINCL) -fPIC -std=c++14
	CFLAGS_CLIB := -O3 -qopenmp -DMKL_ILP64 -I${MKLROOT}/include -fPIC -std=c++14
	
	EXTRA_CFLAGS := -D UTEST -Wcheck -Wall -w2
#-Wstrict-prototypes -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunused-function -Wunused-variable -Wwrite-strings -Wsign-compare -Wdeprecated -Wextra-tokens -Wic-pointer -Wnon-virtual-dtor -Wp64 -Wpointer-arith -Wport -Wreorder -Wreturn-type -Wshadow
	EXTRA_LIBFLAGS := -L$(BDIR_CLIB) -l$(APPNAME)
	EXTRA_LIBFLAGS_2 :=
	
	LIBFLAGS := -shared -qopenmp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
	LIBFLAGS_TEST := -qopenmp -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm
endif

ifdef perf_num
	EXTRA_CFLAGS += -D PERF_TEST
endif

TARGET_PYLIB := $(BDIR_PYLIB)/$(APPNAME)$(PYSUFF)
TARGET_CLIB := $(BDIR_CLIB)/lib$(APPNAME)$(VERSUFF)
TARGET_APP := $(BDIR_APP)/$(APPNAME)_app.exe
TARGET_TEST := $(BDIR_TST)/$(APPNAME)_test.exe

#icpc: $(TARGET_APP) run
icpc: $(TARGET_PYLIB) $(TARGET_CLIB) $(TARGET_TEST) $(TARGET_APP)

SOURCES_ALL := $(wildcard $(SRCDIR)/*.cpp)
SOURCES = $(filter-out $(SRCDIR)/evo_profiling.cpp, $(wildcard $(SRCDIR)/*.cpp))
OBJECT_PYL := $(patsubst $(SRCDIR)/%,$(BDIR_PYLIB)/%,$(SOURCES:.cpp=.o))
OBJECT_LIB := $(filter-out $(BDIR_CLIB)/wrapper.o, $(patsubst $(SRCDIR)/%,$(BDIR_CLIB)/%,$(SOURCES:.cpp=.o)))
OBJECT_APP := $(filter-out $(BDIR_APP)/wrapper.o, $(patsubst $(SRCDIR)/%,$(BDIR_APP)/%,$(SOURCES_ALL:.cpp=.o)))

SOURCES_TST := $(wildcard $(TESTDIR)/$(SRCDIR)/*.cpp)
OBJECT_TST := $(patsubst $(TSRCDIR)/%,$(BDIR_TST)/%,$(SOURCES_TST:.cpp=.o))

#----------------------------------------------------------------

# build lib for python
$(TARGET_PYLIB): $(OBJECT_PYL)
	@echo " Linking python lib ..."
	$(CXX) $^ $(LIBFLAGS) -o $(TARGET_PYLIB) $(EXTRA_LIBFLAGS_2)

# build catch2 test
$(TARGET_TEST): $(OBJECT_TST)
	@echo " Linking tests ..."
	$(CXX) $^ $(LIBFLAGS_TEST) $(EXTRA_LIBFLAGS) $(EXTRA_LIBFLAGS_2) -o $(TARGET_TEST)

# build shared c++ lib
$(TARGET_CLIB): $(OBJECT_LIB)
	@echo " Linking shared lib ..."
	$(CXX) $^ $(LIBFLAGS) -o $(TARGET_CLIB) $(EXTRA_LIBFLAGS_2)

# build c++ app
$(TARGET_APP): $(OBJECT_APP)
	@echo " Linking c++ app ..."
	$(CXX) $^ $(LIBFLAGS_TEST) $(EXTRA_LIBFLAGS) $(EXTRA_LIBFLAGS_2) -o $(TARGET_APP)

#----------------------------------------------------------------

# Compiling python lib
$(BDIR_PYLIB)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling pylib ..."
	@mkdir -p $(BDIR_PYLIB)
	$(CXX) -c $(CFLAGS_PYLIB) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling c++ lib
$(BDIR_CLIB)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling clib ..."
	@mkdir -p $(BDIR_CLIB)
	$(CXX) -c $(CFLAGS_CLIB) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling catch2 tests
$(BDIR_TST)/%.o: $(TSRCDIR)/%.cpp
	@echo " Compilling catch2 tests ..."
	@mkdir -p $(BDIR_TST)
	$(CXX) -c $(CFLAGS_APP) $(EXTRA_CFLAGS) -o $@ $< 

# Compiling c++ app
$(BDIR_APP)/%.o: $(SRCDIR)/%.cpp
	@echo " Compilling c++ app ..."
	@mkdir -p $(BDIR_APP)
	$(CXX) -c $(CFLAGS_APP) $(EXTRA_CFLAGS) -o $@ $< 

#----------------------------------------------------------------

clean:
	@echo " Cleaning..."
	@rm -fr $(BDIR_PYLIB) $(BDIR_CLIB) $(BDIR_TST) $(BDIR_APP) $(TARGET_PYLIB) $(TARGET_TEST) 2>/dev/null || true
	@rm amatrix_* py_out_solution_model_* cpp_solution_model_* 2>/dev/null || true

test:
	@echo " Testing ..."
	python3 $(TESTDIR)/$(SRCDIR)/*.py
	@export  LD_LIBRARY_PATH=$(PWDADDR)/$(BDIR_CLIB):"$(LIBRARY_PATH)"; $(TARGET_TEST)

run:
	@echo " Running ..."
	@export  LD_LIBRARY_PATH=$(PWDADDR)/$(BDIR_CLIB):"$(LIBRARY_PATH)"; $(TARGET_APP)

cleanlogs:
	@echo " Cleaning logs ..."
	@rm *.log

.PHONY: clean
