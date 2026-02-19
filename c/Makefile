# Possible combinations of compiler flags.
# Serial
# Using Indices: -DUSE_INDICES
# Using Intrinsics: -DUSE_INTRINSICS
# Using Intrinsics and Indices: -DUSE_INTRINSICS -DUSE_INDICES
# Using Intrinsics with Contiguous Loop for Storage: -DUSE_INTRINSICS -DCONTIGUOUS_LOOP
# Using Indices and Intrinsics with Contiguous Loop for Storage: -DUSE_INTRINSICS -DUSE_INDICES -DCONTIGUOUS_LOOP

# The clear winner -> -O2 -march=native -DUSE_INTRINSICS -DUSE_INDICES -DCONTIGUOUS_LOOP

src_dir=src
include_dir=include
object_dir=obj
bin_dir=bin
lib_dir=lib

validate_obj=grid.o intrinsics_utils.o transitionfunction.o validate.o
benchmark_obj=grid.o intrinsics_utils.o transitionfunction.o benchmark.o benchmark_utils.o
library_obj=grid.o intrinsics_utils.o transitionfunction.o

validate_obj:=$(validate_obj:%=$(object_dir)/%)
benchmark_obj:=$(benchmark_obj:%=$(object_dir)/%)
library_obj:=$(library_obj:%=$(object_dir)/%)

cc=gcc
ccflags=-fPIC -O2 -march=native -I$(include_dir) -DUSE_INTRINSICS -DUSE_INDICES -DCONTIGUOUS_LOOP
ldflags=

so_ccflags=-O2 -march=native -I$(include_dir) -DUSE_INTRINSICS -DUSE_INDICES -DCONTIGUOUS_LOOP -fPIC
so_ldflags=-shared

# Debug flags
#ccflags=-O2 -march=native -g -pg -no-pie
#ldflags=-pg -no-pie

validate: $(validate_obj)
	$(cc) $^ -o $(bin_dir)/validate  $(ldflags)

benchmark: $(benchmark_obj)
	$(cc) $^ -o $(bin_dir)/benchmark $(ldflags)

smoothlife: $(library_obj)
	$(cc) $^ -o $(lib_dir)/smoothlife.so $(so_ldflags)

$(object_dir)/benchmark.o: $(src_dir)/benchmark.c
	$(cc) -c $< $(ccflags) -o $@ 

$(object_dir)/validate.o: $(src_dir)/validate.c
	$(cc) -c $< $(ccflags) -o $@ 

$(object_dir)/grid.o: $(src_dir)/grid.c $(include_dir)/grid.h
	$(cc) -c $< $(ccflags) -o $@ 

$(object_dir)/intrinsics_utils.o: $(src_dir)/intrinsics_utils.c $(include_dir)/intrinsics_utils.h
	$(cc) -c $< $(ccflags) -o $@ 

$(object_dir)/benchmark_utils.o: $(src_dir)/benchmark_utils.c $(include_dir)/benchmark_utils.h
	$(cc) -c $< $(ccflags) -o $@ 

$(object_dir)/transitionfunction.o: $(src_dir)/transitionfunction.c $(include_dir)/transitionfunction.h
	$(cc) -c $< $(ccflags) -o $@ 

$(object_dir):
	mkdir -p $(object_dir)

$(bin_dir):
	mkdir -p $(bin_dir)

$(lib_dir):
	mkdir -p $(lib_dir)

clean:
	rm -rf $(object_dir)/* $(bin_dir)/* $(lib_dir)/*

run_validate:
	$(bin_dir)/validate.exe

run_benchmark:
	$(bin_dir)/benchmark.exe

.PHONY: clean, run_validate, run_benchmark
