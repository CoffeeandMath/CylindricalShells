# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/build

# Include any dependencies generated for this target.
include CMakeFiles/materialclass.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/materialclass.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/materialclass.dir/flags.make

CMakeFiles/materialclass.dir/src/material_class.cc.o: CMakeFiles/materialclass.dir/flags.make
CMakeFiles/materialclass.dir/src/material_class.cc.o: ../src/material_class.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/materialclass.dir/src/material_class.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/materialclass.dir/src/material_class.cc.o -c /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/src/material_class.cc

CMakeFiles/materialclass.dir/src/material_class.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/materialclass.dir/src/material_class.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/src/material_class.cc > CMakeFiles/materialclass.dir/src/material_class.cc.i

CMakeFiles/materialclass.dir/src/material_class.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/materialclass.dir/src/material_class.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/src/material_class.cc -o CMakeFiles/materialclass.dir/src/material_class.cc.s

# Object files for target materialclass
materialclass_OBJECTS = \
"CMakeFiles/materialclass.dir/src/material_class.cc.o"

# External object files for target materialclass
materialclass_EXTERNAL_OBJECTS =

libmaterialclass.a: CMakeFiles/materialclass.dir/src/material_class.cc.o
libmaterialclass.a: CMakeFiles/materialclass.dir/build.make
libmaterialclass.a: CMakeFiles/materialclass.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libmaterialclass.a"
	$(CMAKE_COMMAND) -P CMakeFiles/materialclass.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/materialclass.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/materialclass.dir/build: libmaterialclass.a

.PHONY : CMakeFiles/materialclass.dir/build

CMakeFiles/materialclass.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/materialclass.dir/cmake_clean.cmake
.PHONY : CMakeFiles/materialclass.dir/clean

CMakeFiles/materialclass.dir/depend:
	cd /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/build /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/build /home/kevin/Documents/Research/CylindricalShells/CyclindricalFull3D/build/CMakeFiles/materialclass.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/materialclass.dir/depend

