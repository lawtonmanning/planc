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
CMAKE_COMMAND = /autofs/nccs-svm1_sw/rhea/.swci/0-core/opt/spack/20191017/linux-rhel7-x86_64/gcc-4.8.5/cmake-3.15.2-4xbxfsztfb7z2qemxy2dgrsuabbdgifu/bin/cmake

# The command to remove a file.
RM = /autofs/nccs-svm1_sw/rhea/.swci/0-core/opt/spack/20191017/linux-rhel7-x86_64/gcc-4.8.5/cmake-3.15.2-4xbxfsztfb7z2qemxy2dgrsuabbdgifu/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /ccs/home/mannlg15/planc

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /ccs/home/mannlg15/planc/build-rhea

# Include any dependencies generated for this target.
include dense_distntf/CMakeFiles/dense_distntf.dir/depend.make

# Include the progress variables for this target.
include dense_distntf/CMakeFiles/dense_distntf.dir/progress.make

# Include the compile flags for this target's objects.
include dense_distntf/CMakeFiles/dense_distntf.dir/flags.make

dense_distntf/CMakeFiles/dense_distntf.dir/distntf.cpp.o: dense_distntf/CMakeFiles/dense_distntf.dir/flags.make
dense_distntf/CMakeFiles/dense_distntf.dir/distntf.cpp.o: ../distntf/distntf.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/ccs/home/mannlg15/planc/build-rhea/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object dense_distntf/CMakeFiles/dense_distntf.dir/distntf.cpp.o"
	cd /ccs/home/mannlg15/planc/build-rhea/dense_distntf && /ccs/compilers/gcc/rhel7-x86_64/6.2.0/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/dense_distntf.dir/distntf.cpp.o -c /ccs/home/mannlg15/planc/distntf/distntf.cpp

dense_distntf/CMakeFiles/dense_distntf.dir/distntf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dense_distntf.dir/distntf.cpp.i"
	cd /ccs/home/mannlg15/planc/build-rhea/dense_distntf && /ccs/compilers/gcc/rhel7-x86_64/6.2.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /ccs/home/mannlg15/planc/distntf/distntf.cpp > CMakeFiles/dense_distntf.dir/distntf.cpp.i

dense_distntf/CMakeFiles/dense_distntf.dir/distntf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dense_distntf.dir/distntf.cpp.s"
	cd /ccs/home/mannlg15/planc/build-rhea/dense_distntf && /ccs/compilers/gcc/rhel7-x86_64/6.2.0/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /ccs/home/mannlg15/planc/distntf/distntf.cpp -o CMakeFiles/dense_distntf.dir/distntf.cpp.s

# Object files for target dense_distntf
dense_distntf_OBJECTS = \
"CMakeFiles/dense_distntf.dir/distntf.cpp.o"

# External object files for target dense_distntf
dense_distntf_EXTERNAL_OBJECTS =

dense_distntf/dense_distntf: dense_distntf/CMakeFiles/dense_distntf.dir/distntf.cpp.o
dense_distntf/dense_distntf: dense_distntf/CMakeFiles/dense_distntf.dir/build.make
dense_distntf/dense_distntf: /autofs/nccs-svm1_sw/rhea/.swci/0-core/opt/spack/20191017/linux-rhel7-x86_64/gcc-6.2.0/openblas-0.3.5-pqhmw4pblzsh76ejrzz5hz2o37kl3esm/lib/libopenblas.so
dense_distntf/dense_distntf: /autofs/nccs-svm1_sw/rhea/.swci/0-core/opt/spack/20191017/linux-rhel7-x86_64/gcc-6.2.0/openblas-0.3.5-pqhmw4pblzsh76ejrzz5hz2o37kl3esm/lib/libopenblas.so
dense_distntf/dense_distntf: /autofs/nccs-svm1_sw/rhea/.swci/0-core/opt/spack/20191017/linux-rhel7-x86_64/gcc-6.2.0/openmpi-3.1.4-6i3xptfa3r6joouympddcnja5njc3eda/lib/libmpi_cxx.so
dense_distntf/dense_distntf: /autofs/nccs-svm1_sw/rhea/.swci/0-core/opt/spack/20191017/linux-rhel7-x86_64/gcc-6.2.0/openmpi-3.1.4-6i3xptfa3r6joouympddcnja5njc3eda/lib/libmpi.so
dense_distntf/dense_distntf: dense_distntf/CMakeFiles/dense_distntf.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/ccs/home/mannlg15/planc/build-rhea/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable dense_distntf"
	cd /ccs/home/mannlg15/planc/build-rhea/dense_distntf && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dense_distntf.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dense_distntf/CMakeFiles/dense_distntf.dir/build: dense_distntf/dense_distntf

.PHONY : dense_distntf/CMakeFiles/dense_distntf.dir/build

dense_distntf/CMakeFiles/dense_distntf.dir/clean:
	cd /ccs/home/mannlg15/planc/build-rhea/dense_distntf && $(CMAKE_COMMAND) -P CMakeFiles/dense_distntf.dir/cmake_clean.cmake
.PHONY : dense_distntf/CMakeFiles/dense_distntf.dir/clean

dense_distntf/CMakeFiles/dense_distntf.dir/depend:
	cd /ccs/home/mannlg15/planc/build-rhea && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /ccs/home/mannlg15/planc /ccs/home/mannlg15/planc/distntf /ccs/home/mannlg15/planc/build-rhea /ccs/home/mannlg15/planc/build-rhea/dense_distntf /ccs/home/mannlg15/planc/build-rhea/dense_distntf/CMakeFiles/dense_distntf.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : dense_distntf/CMakeFiles/dense_distntf.dir/depend

