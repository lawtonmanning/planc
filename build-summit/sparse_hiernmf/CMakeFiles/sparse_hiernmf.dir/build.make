# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /autofs/nccs-svm1_sw/summit/.swci/0-core/opt/spack/20180914/linux-rhel7-ppc64le/gcc-4.8.5/cmake-3.15.2-xit2o3iepxvqbyku77lwcugufilztu7t/bin/cmake

# The command to remove a file.
RM = /autofs/nccs-svm1_sw/summit/.swci/0-core/opt/spack/20180914/linux-rhel7-ppc64le/gcc-4.8.5/cmake-3.15.2-xit2o3iepxvqbyku77lwcugufilztu7t/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /ccs/home/mannlg15/planc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ccs/home/mannlg15/planc/build-summit

# Include any dependencies generated for this target.
include sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/depend.make

# Include the progress variables for this target.
include sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/progress.make

# Include the compile flags for this target's objects.
include sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/flags.make

sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.o: sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/flags.make
sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.o: ../hiernmf/hiernmf.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/ccs/home/mannlg15/planc/build-summit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.o"
	cd /ccs/home/mannlg15/planc/build-summit/sparse_hiernmf && /sw/summit/gcc/6.4.0/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.o -c /ccs/home/mannlg15/planc/hiernmf/hiernmf.cpp

sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.i"
	cd /ccs/home/mannlg15/planc/build-summit/sparse_hiernmf && /sw/summit/gcc/6.4.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /ccs/home/mannlg15/planc/hiernmf/hiernmf.cpp > CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.i

sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.s"
	cd /ccs/home/mannlg15/planc/build-summit/sparse_hiernmf && /sw/summit/gcc/6.4.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /ccs/home/mannlg15/planc/hiernmf/hiernmf.cpp -o CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.s

# Object files for target sparse_hiernmf
sparse_hiernmf_OBJECTS = \
"CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.o"

# External object files for target sparse_hiernmf
sparse_hiernmf_EXTERNAL_OBJECTS =

sparse_hiernmf/sparse_hiernmf: sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/hiernmf.cpp.o
sparse_hiernmf/sparse_hiernmf: sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/build.make
sparse_hiernmf/sparse_hiernmf: /sw/summit/cuda/10.1.243/lib64/libnvblas.so
sparse_hiernmf/sparse_hiernmf: /sw/summit/cuda/10.1.243/lib64/libcublas.so
sparse_hiernmf/sparse_hiernmf: /sw/summit/cuda/10.1.243/lib64/libcudart_static.a
sparse_hiernmf/sparse_hiernmf: /usr/lib64/librt.so
sparse_hiernmf/sparse_hiernmf: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/openblas-0.3.9-g7tkwn4kc7ukzxqqqia4jwpqr27aw3cu/lib/libopenblas.so
sparse_hiernmf/sparse_hiernmf: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/openblas-0.3.9-g7tkwn4kc7ukzxqqqia4jwpqr27aw3cu/lib/libopenblas.so
sparse_hiernmf/sparse_hiernmf: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/spectrum-mpi-10.3.1.2-20200121-awz2q5brde7wgdqqw4ugalrkukeub4eb/lib/libmpiprofilesupport.so
sparse_hiernmf/sparse_hiernmf: /autofs/nccs-svm1_sw/summit/.swci/1-compute/opt/spack/20180914/linux-rhel7-ppc64le/gcc-6.4.0/spectrum-mpi-10.3.1.2-20200121-awz2q5brde7wgdqqw4ugalrkukeub4eb/lib/libmpi_ibm.so
sparse_hiernmf/sparse_hiernmf: sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/ccs/home/mannlg15/planc/build-summit/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable sparse_hiernmf"
	cd /ccs/home/mannlg15/planc/build-summit/sparse_hiernmf && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sparse_hiernmf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/build: sparse_hiernmf/sparse_hiernmf

.PHONY : sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/build

sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/clean:
	cd /ccs/home/mannlg15/planc/build-summit/sparse_hiernmf && $(CMAKE_COMMAND) -P CMakeFiles/sparse_hiernmf.dir/cmake_clean.cmake
.PHONY : sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/clean

sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/depend:
	cd /ccs/home/mannlg15/planc/build-summit && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ccs/home/mannlg15/planc /ccs/home/mannlg15/planc/hiernmf /ccs/home/mannlg15/planc/build-summit /ccs/home/mannlg15/planc/build-summit/sparse_hiernmf /ccs/home/mannlg15/planc/build-summit/sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : sparse_hiernmf/CMakeFiles/sparse_hiernmf.dir/depend

