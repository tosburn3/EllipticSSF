# EllipticSSF

AUTHORS: Thomas Osburn and Nami Nishimura 

This Matehmatica package uses the finite difference method to solve for a massless scalar field in Kerr space time and calculate the associated self-force on the scalar charge. It uses Barry Wardell's effective source code (https://github.com/barrywardell/EffectiveSource) for regularization of the scalar field near the charge. The methods involve separating the "t" and "phi" variables, which leads to an elliptic PDE for each m-mode of the scalar field.


GENERAL COMPILATION INFO

Because Barry's code is written in C, we use the WSTP (formerly known as MathLink) feature of Mathematica to compile an executable with effective source functions that are callable from Mathematica code. Compiling the C code following WSTP requirements is the most challenging part of the setup. Executables have been included for x86-64 Windows and x86-64 Linux, but not for Mac (and it may be that the provided executables don't work on your system and you would need to compile them). To see if a provided executable works for you, rename the executable to remove "_WINDOWS" or "_LINUX" from the filename and skip ahead to the marked discussion about the example notebook.

### If a provided executable works for you skip ahead to the marked discussion about the example notebook.

A dependency for Barry's code is the GNU Scientific Library (https://www.gnu.org/software/gsl/), install that before you proceed.

Wolfram's guide for compiling with WSTP has sections for all major system types: https://reference.wolfram.com/language/tutorial/WSTPAndExternalProgramCommunicationOverview.html


COMPILING ON WINDOWS

Visual Studio (not Visual Studio Code) is required for compilation on Windows. The easiest way to install GSL is through vcpkg (https://solarianprogrammer.com/2020/01/26/getting-started-gsl-gnu-scientific-library-windows-macos-linux/#gsl_installation_windows). To interface WSTP with Visual Studio, you have to copy some files from the Mathematica directories into the Visual Studio directories, see "Installing the WSTP Components" in the WSTP guide (https://reference.wolfram.com/language/tutorial/WSTPDeveloperGuide-Windows.html, and here is a helpful reference for more recent versions of Visual Studio: https://mathematica.stackexchange.com/questions/232846/where-to-place-the-essential-wstp-files-in-visual-studio-2019). Then you need to open the very specific commanld line application "VS20XX x86 Native Tools Command Prompt", which is part of Visual Studio. Using that command prompt, here is an example set of commands to compile on Windows:

SET CL=/nologo /c /DWIN32 /D_WINDOWS /W3 /O2 /DNDEBUG

SET LINK=/NOLOGO /SUBSYSTEM:windows /INCREMENTAL:no /PDB:NONE kernel32.lib user32.lib gdi32.lib

WSPREP kerr-circular.tm -o testtm.cpp

CL -I"C:/DEV/vcpkg/installed/x64-windows/include/" kerr-circular.cpp testtm.cpp

LINK kerr-circular.obj testtm.obj wstp64i4m.lib C:/Users/tosbu/vcpkg/installed/x64-windows/lib/gsl.lib /OUT:kerr-circular.exe


COMPILING ON LINUX/MAC

Idk (I've done it once on each type of system, but I don't have a record of the commands, see WSTP instructions)


USING THE EXAMPLE NOTEBOOK

The example notebook should now be able to link with "kerr-circular.exe" by including the correct path to the executable. See example notebook for usage details.
