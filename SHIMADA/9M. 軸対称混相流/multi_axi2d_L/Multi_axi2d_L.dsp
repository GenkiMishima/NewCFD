# Microsoft Developer Studio Project File - Name="Multi_axi2d_L" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=Multi_axi2d_L - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "Multi_axi2d_L.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "Multi_axi2d_L.mak" CFG="Multi_axi2d_L - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "Multi_axi2d_L - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "Multi_axi2d_L - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "Multi_axi2d_L - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x411 /d "NDEBUG"
# ADD RSC /l 0x411 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "Multi_axi2d_L - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x411 /d "_DEBUG"
# ADD RSC /l 0x411 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "Multi_axi2d_L - Win32 Release"
# Name "Multi_axi2d_L - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=.\geomet.f90
DEP_F90_GEOME=\
	".\common.h"\
	
# End Source File
# Begin Source File

SOURCE=.\get_data.f90
DEP_F90_GET_D=\
	".\common.h"\
	
# End Source File
# Begin Source File

SOURCE=.\I_boundary.f90
DEP_F90_I_BOU=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# Begin Source File

SOURCE=.\I_Flux.f90
DEP_F90_I_FLU=\
	".\common.h"\
	".\function.h"\
	".\vanLeer2d.h"\
	
# End Source File
# Begin Source File

SOURCE=.\init_field.f90
DEP_F90_INIT_=\
	".\common.h"\
	
# End Source File
# Begin Source File

SOURCE=.\init_io.f90
DEP_F90_INIT_I=\
	".\common.h"\
	
# End Source File
# Begin Source File

SOURCE=.\init_param.f90
DEP_F90_INIT_P=\
	".\common.h"\
	
# End Source File
# Begin Source File

SOURCE=.\J_boundary.f90
DEP_F90_J_BOU=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# Begin Source File

SOURCE=.\J_Flux.f90
DEP_F90_J_FLU=\
	".\common.h"\
	".\function.h"\
	".\vanLeer2d.h"\
	
# End Source File
# Begin Source File

SOURCE=.\LU_SGS.f90
DEP_F90_LU_SG=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# Begin Source File

SOURCE=.\main.f90
DEP_F90_MAIN_=\
	".\common.h"\
	
# End Source File
# Begin Source File

SOURCE=.\out_field.f90
DEP_F90_OUT_F=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# Begin Source File

SOURCE=.\Point_Implicit.f90
DEP_F90_POINT=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# Begin Source File

SOURCE=.\show_time_residual.f90
DEP_F90_SHOW_=\
	".\common.h"\
	
# End Source File
# Begin Source File

SOURCE=.\source.f90
DEP_F90_SOURC=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# Begin Source File

SOURCE=.\time_step.f90
DEP_F90_TIME_=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# Begin Source File

SOURCE=.\update.f90
DEP_F90_UPDAT=\
	".\common.h"\
	".\function.h"\
	
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\common.h
# End Source File
# Begin Source File

SOURCE=.\function.h
# End Source File
# Begin Source File

SOURCE=.\roe2d.h
# End Source File
# Begin Source File

SOURCE=.\vanLeer2d.h
# End Source File
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# Begin Source File

SOURCE=.\grid.data
# End Source File
# Begin Source File

SOURCE=.\input.data
# End Source File
# End Group
# End Target
# End Project
