# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework

# Include any dependencies generated for this target.
include CMakeFiles/LB2D.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/LB2D.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LB2D.dir/flags.make

CMakeFiles/LB2D.dir/main.cpp.o: CMakeFiles/LB2D.dir/flags.make
CMakeFiles/LB2D.dir/main.cpp.o: main.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/LB2D.dir/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/LB2D.dir/main.cpp.o -c /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/main.cpp

CMakeFiles/LB2D.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LB2D.dir/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/main.cpp > CMakeFiles/LB2D.dir/main.cpp.i

CMakeFiles/LB2D.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LB2D.dir/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/main.cpp -o CMakeFiles/LB2D.dir/main.cpp.s

CMakeFiles/LB2D.dir/main.cpp.o.requires:
.PHONY : CMakeFiles/LB2D.dir/main.cpp.o.requires

CMakeFiles/LB2D.dir/main.cpp.o.provides: CMakeFiles/LB2D.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/LB2D.dir/build.make CMakeFiles/LB2D.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/LB2D.dir/main.cpp.o.provides

CMakeFiles/LB2D.dir/main.cpp.o.provides.build: CMakeFiles/LB2D.dir/main.cpp.o

CMakeFiles/LB2D.dir/geometry_2D.cpp.o: CMakeFiles/LB2D.dir/flags.make
CMakeFiles/LB2D.dir/geometry_2D.cpp.o: geometry_2D.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/LB2D.dir/geometry_2D.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/LB2D.dir/geometry_2D.cpp.o -c /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/geometry_2D.cpp

CMakeFiles/LB2D.dir/geometry_2D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LB2D.dir/geometry_2D.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/geometry_2D.cpp > CMakeFiles/LB2D.dir/geometry_2D.cpp.i

CMakeFiles/LB2D.dir/geometry_2D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LB2D.dir/geometry_2D.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/geometry_2D.cpp -o CMakeFiles/LB2D.dir/geometry_2D.cpp.s

CMakeFiles/LB2D.dir/geometry_2D.cpp.o.requires:
.PHONY : CMakeFiles/LB2D.dir/geometry_2D.cpp.o.requires

CMakeFiles/LB2D.dir/geometry_2D.cpp.o.provides: CMakeFiles/LB2D.dir/geometry_2D.cpp.o.requires
	$(MAKE) -f CMakeFiles/LB2D.dir/build.make CMakeFiles/LB2D.dir/geometry_2D.cpp.o.provides.build
.PHONY : CMakeFiles/LB2D.dir/geometry_2D.cpp.o.provides

CMakeFiles/LB2D.dir/geometry_2D.cpp.o.provides.build: CMakeFiles/LB2D.dir/geometry_2D.cpp.o

# Object files for target LB2D
LB2D_OBJECTS = \
"CMakeFiles/LB2D.dir/main.cpp.o" \
"CMakeFiles/LB2D.dir/geometry_2D.cpp.o"

# External object files for target LB2D
LB2D_EXTERNAL_OBJECTS =

LB2D: CMakeFiles/LB2D.dir/main.cpp.o
LB2D: CMakeFiles/LB2D.dir/geometry_2D.cpp.o
LB2D: CMakeFiles/LB2D.dir/build.make
LB2D: /usr/lib64/libGLEW.so
LB2D: /usr/lib/x86_64-linux-gnu/libmpfr.so
LB2D: /usr/lib/x86_64-linux-gnu/libgmp.so
LB2D: /usr/lib/libCGAL_Core.so
LB2D: /usr/lib/libCGAL.so
LB2D: /usr/lib/x86_64-linux-gnu/libboost_thread.so
LB2D: /usr/lib/x86_64-linux-gnu/libboost_system.so
LB2D: /usr/lib/x86_64-linux-gnu/libpthread.so
LB2D: /usr/lib/x86_64-linux-gnu/libGLU.so
LB2D: /usr/lib/x86_64-linux-gnu/libGL.so
LB2D: /usr/lib/x86_64-linux-gnu/libSM.so
LB2D: /usr/lib/x86_64-linux-gnu/libICE.so
LB2D: /usr/lib/x86_64-linux-gnu/libX11.so
LB2D: /usr/lib/x86_64-linux-gnu/libXext.so
LB2D: /usr/lib/x86_64-linux-gnu/libglut.so
LB2D: /usr/lib/x86_64-linux-gnu/libXmu.so
LB2D: /usr/lib/x86_64-linux-gnu/libXi.so
LB2D: libBMP.a
LB2D: /usr/lib64/libGLEW.so
LB2D: CMakeFiles/LB2D.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable LB2D"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LB2D.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LB2D.dir/build: LB2D
.PHONY : CMakeFiles/LB2D.dir/build

CMakeFiles/LB2D.dir/requires: CMakeFiles/LB2D.dir/main.cpp.o.requires
CMakeFiles/LB2D.dir/requires: CMakeFiles/LB2D.dir/geometry_2D.cpp.o.requires
.PHONY : CMakeFiles/LB2D.dir/requires

CMakeFiles/LB2D.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/LB2D.dir/cmake_clean.cmake
.PHONY : CMakeFiles/LB2D.dir/clean

CMakeFiles/LB2D.dir/depend:
	cd /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework /home/frank/Desktop/MyStudies/Lattice_Boltzmann/LB2D_Framework/CMakeFiles/LB2D.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LB2D.dir/depend
