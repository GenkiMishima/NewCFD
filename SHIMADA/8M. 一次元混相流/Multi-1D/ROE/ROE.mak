# Microsoft Developer Studio Generated NMAKE File, Based on ROE.dsp
!IF "$(CFG)" == ""
CFG=ROE - Win32 Debug
!MESSAGE No configuration specified. Defaulting to ROE - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "ROE - Win32 Release" && "$(CFG)" != "ROE - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ROE.mak" CFG="ROE - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ROE - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "ROE - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 

CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "ROE - Win32 Release"

OUTDIR=.\Release
INTDIR=.\Release
# Begin Custom Macros
OutDir=.\Release
# End Custom Macros

ALL : "$(OUTDIR)\ROE.exe"


CLEAN :
	-@erase "$(INTDIR)\ROE.obj"
	-@erase "$(OUTDIR)\ROE.exe"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/compile_only /nologo /warn:nofileopt /module:"Release/" /object:"Release/" 
F90_OBJS=.\Release/
CPP_PROJ=/nologo /ML /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fp"$(INTDIR)\ROE.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\ROE.bsc" 
BSC32_SBRS= \
	
LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\ROE.pdb" /machine:I386 /out:"$(OUTDIR)\ROE.exe" 
LINK32_OBJS= \
	"$(INTDIR)\ROE.obj"

"$(OUTDIR)\ROE.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "ROE - Win32 Debug"

OUTDIR=.\Debug
INTDIR=.\Debug
# Begin Custom Macros
OutDir=.\Debug
# End Custom Macros

ALL : "$(OUTDIR)\ROE.exe" "$(OUTDIR)\ROE.bsc"


CLEAN :
	-@erase "$(INTDIR)\DF60.PDB"
	-@erase "$(INTDIR)\ROE.obj"
	-@erase "$(INTDIR)\ROE.sbr"
	-@erase "$(OUTDIR)\ROE.bsc"
	-@erase "$(OUTDIR)\ROE.exe"
	-@erase "$(OUTDIR)\ROE.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

F90_PROJ=/browser:"Debug/" /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt /module:"Debug/" /object:"Debug/" /pdbfile:"Debug/DF60.PDB" 
F90_OBJS=.\Debug/
CPP_PROJ=/nologo /MLd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR"$(INTDIR)\\" /Fp"$(INTDIR)\ROE.pch" /YX /Fo"$(INTDIR)\\" /Fd"$(INTDIR)\\" /FD /GZ  /c 
BSC32=bscmake.exe
BSC32_FLAGS=/nologo /o"$(OUTDIR)\ROE.bsc" 
BSC32_SBRS= \
	"$(INTDIR)\ROE.sbr"

"$(OUTDIR)\ROE.bsc" : "$(OUTDIR)" $(BSC32_SBRS)
    $(BSC32) @<<
  $(BSC32_FLAGS) $(BSC32_SBRS)
<<

LINK32=link.exe
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no /pdb:"$(OUTDIR)\ROE.pdb" /debug /machine:I386 /out:"$(OUTDIR)\ROE.exe" /pdbtype:sept 
LINK32_OBJS= \
	"$(INTDIR)\ROE.obj"

"$(OUTDIR)\ROE.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.c{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.obj::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.c{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cpp{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.cxx{$(INTDIR)}.sbr::
   $(CPP) @<<
   $(CPP_PROJ) $< 
<<

.SUFFIXES: .fpp

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.fpp{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  


!IF "$(NO_EXTERNAL_DEPS)" != "1"
!IF EXISTS("ROE.dep")
!INCLUDE "ROE.dep"
!ELSE 
!MESSAGE Warning: cannot find "ROE.dep"
!ENDIF 
!ENDIF 


!IF "$(CFG)" == "ROE - Win32 Release" || "$(CFG)" == "ROE - Win32 Debug"
SOURCE=.\ROE.f90

!IF  "$(CFG)" == "ROE - Win32 Release"


"$(INTDIR)\ROE.obj" : $(SOURCE) "$(INTDIR)"


!ELSEIF  "$(CFG)" == "ROE - Win32 Debug"


"$(INTDIR)\ROE.obj"	"$(INTDIR)\ROE.sbr" : $(SOURCE) "$(INTDIR)"


!ENDIF 


!ENDIF 

