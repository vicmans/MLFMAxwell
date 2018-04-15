################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../main.f90 \
../modFMM.f90 \
../modGlobalMetod.f90 \
../modGlobalParam.f90 \
../modIterativo.f90 \
../modMLFMA.f90 \
../modMoM.f90 \
../modRCS.f90 \
../modRWG.f90 \
../modReadStl.f90 

OBJS += \
./main.o \
./modFMM.o \
./modGlobalMetod.o \
./modGlobalParam.o \
./modIterativo.o \
./modMLFMA.o \
./modMoM.o \
./modRCS.o \
./modRWG.o \
./modReadStl.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -g -O0 -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

main.o: ../main.f90 modFMM.o modGlobalMetod.o modGlobalParam.o modIterativo.o modMLFMA.o modMoM.o modRCS.o modRWG.o modReadStl.o

modFMM.o: ../modFMM.f90 modGlobalMetod.o modGlobalParam.o modMoM.o

modGlobalMetod.o: ../modGlobalMetod.f90 modGlobalParam.o

modGlobalParam.o: ../modGlobalParam.f90

modIterativo.o: ../modIterativo.f90 modGlobalMetod.o modGlobalParam.o

modMLFMA.o: ../modMLFMA.f90 modGlobalMetod.o modGlobalParam.o modIterativo.o modMoM.o

modMoM.o: ../modMoM.f90 modGlobalMetod.o modGlobalParam.o

modRCS.o: ../modRCS.f90 modGlobalMetod.o modGlobalParam.o

modRWG.o: ../modRWG.f90 modGlobalParam.o

modReadStl.o: ../modReadStl.f90 modGlobalParam.o


