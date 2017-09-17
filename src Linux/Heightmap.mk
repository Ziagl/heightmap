##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=Heightmap
ConfigurationName      :=Debug
WorkspacePath          :=/home/werner/Cpp_Projekte/Heightmap
ProjectPath            :=/home/werner/Cpp_Projekte/Heightmap
IntermediateDirectory  :=./Debug
OutDir                 := $(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Werner Ziegelwanger
Date                   :=06/09/17
CodeLitePath           :=/home/werner/.codelite
LinkerName             :=/usr/bin/g++
SharedObjectLinkerName :=/usr/bin/g++ -shared -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.i
DebugSwitch            :=-g 
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputFile             :=$(IntermediateDirectory)/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E
ObjectsFileList        :="Heightmap.txt"
PCHCompileFlags        :=
MakeDirCommand         :=mkdir -p
LinkOptions            :=  
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overriden using an environment variables
##
AR       := /usr/bin/ar rcu
CXX      := /usr/bin/g++
CC       := /usr/bin/gcc
CXXFLAGS :=  -g -O0 -Wall $(Preprocessors)
CFLAGS   :=  -g -O0 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := /usr/bin/as


##
## User defined environment variables
##
CodeLiteDir:=/usr/share/codelite
Objects0=$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IntermediateDirectory)/PerlinNoise.cpp$(ObjectSuffix) $(IntermediateDirectory)/TerrainGenerator.cpp$(ObjectSuffix) $(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) $(@D)
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@test -d ./Debug || $(MakeDirCommand) ./Debug


$(IntermediateDirectory)/.d:
	@test -d ./Debug || $(MakeDirCommand) ./Debug

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/main.cpp$(ObjectSuffix): main.cpp $(IntermediateDirectory)/main.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/werner/Cpp_Projekte/Heightmap/main.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/main.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/main.cpp$(DependSuffix): main.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/main.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/main.cpp$(DependSuffix) -MM main.cpp

$(IntermediateDirectory)/main.cpp$(PreprocessSuffix): main.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/main.cpp$(PreprocessSuffix) main.cpp

$(IntermediateDirectory)/PerlinNoise.cpp$(ObjectSuffix): PerlinNoise.cpp $(IntermediateDirectory)/PerlinNoise.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/werner/Cpp_Projekte/Heightmap/PerlinNoise.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/PerlinNoise.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/PerlinNoise.cpp$(DependSuffix): PerlinNoise.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/PerlinNoise.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/PerlinNoise.cpp$(DependSuffix) -MM PerlinNoise.cpp

$(IntermediateDirectory)/PerlinNoise.cpp$(PreprocessSuffix): PerlinNoise.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/PerlinNoise.cpp$(PreprocessSuffix) PerlinNoise.cpp

$(IntermediateDirectory)/TerrainGenerator.cpp$(ObjectSuffix): TerrainGenerator.cpp $(IntermediateDirectory)/TerrainGenerator.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/werner/Cpp_Projekte/Heightmap/TerrainGenerator.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/TerrainGenerator.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/TerrainGenerator.cpp$(DependSuffix): TerrainGenerator.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/TerrainGenerator.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/TerrainGenerator.cpp$(DependSuffix) -MM TerrainGenerator.cpp

$(IntermediateDirectory)/TerrainGenerator.cpp$(PreprocessSuffix): TerrainGenerator.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/TerrainGenerator.cpp$(PreprocessSuffix) TerrainGenerator.cpp

$(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(ObjectSuffix): EasyBMP/EasyBMP.cpp $(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(DependSuffix)
	$(CXX) $(IncludePCH) $(SourceSwitch) "/home/werner/Cpp_Projekte/Heightmap/EasyBMP/EasyBMP.cpp" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(DependSuffix): EasyBMP/EasyBMP.cpp
	@$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) -MG -MP -MT$(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(ObjectSuffix) -MF$(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(DependSuffix) -MM EasyBMP/EasyBMP.cpp

$(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(PreprocessSuffix): EasyBMP/EasyBMP.cpp
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/EasyBMP_EasyBMP.cpp$(PreprocessSuffix) EasyBMP/EasyBMP.cpp


-include $(IntermediateDirectory)/*$(DependSuffix)
##
## Clean
##
clean:
	$(RM) -r ./Debug/


