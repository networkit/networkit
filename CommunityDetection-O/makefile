################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/viz/test/subdir.mk
-include src/viz/subdir.mk
-include src/test/subdir.mk
-include src/spanning/test/subdir.mk
-include src/spanning/subdir.mk
-include src/scoring/subdir.mk
-include src/scd/test/subdir.mk
-include src/scd/subdir.mk
-include src/properties/test/subdir.mk
-include src/properties/subdir.mk
-include src/partitioning/test/subdir.mk
-include src/partitioning/subdir.mk
-include src/overlap/test/subdir.mk
-include src/overlap/subdir.mk
-include src/matching/test/subdir.mk
-include src/matching/subdir.mk
-include src/mapping/test/subdir.mk
-include src/mapping/subdir.mk
-include src/io/test/subdir.mk
-include src/io/subdir.mk
-include src/independentset/test/subdir.mk
-include src/independentset/subdir.mk
-include src/graph/test/subdir.mk
-include src/graph/subdir.mk
-include src/generators/test/subdir.mk
-include src/generators/subdir.mk
-include src/dynamics/subdir.mk
-include src/distmeasures/test/subdir.mk
-include src/distmeasures/subdir.mk
-include src/dcd/test/subdir.mk
-include src/dcd/subdir.mk
-include src/connectedcomponents/test/subdir.mk
-include src/connectedcomponents/subdir.mk
-include src/community/test/subdir.mk
-include src/community/subdir.mk
-include src/coarsening/test/subdir.mk
-include src/coarsening/subdir.mk
-include src/clustering/test/subdir.mk
-include src/clustering/subdir.mk
-include src/base/subdir.mk
-include src/auxiliary/test/subdir.mk
-include src/auxiliary/subdir.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: NetworKit-CommunityDetection-O

# Tool invocations
NetworKit-CommunityDetection-O: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X C++ Linker'
	/usr/local/bin/g++-4.7 -L"/Users/cls/workspace/gtest" -L/usr/local/Cellar/log4cxx/0.10.0/lib -fopenmp -std=c++11 -o "NetworKit-CommunityDetection-O" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) NetworKit-CommunityDetection-O
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets
