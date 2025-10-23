# EllipticSSF

Authors: Thomas Osburn and Nami Nishimura 

This Matehmatica package uses the finite difference method to solve for a massless scalar field in Kerr space time and calculate the associated self-force on the scalar charge. It uses Barry Wardell's effective source code (https://github.com/barrywardell/EffectiveSource) for regularization of the scalar field near the charge. The methods involve separating the "t" and "phi" variables, which leads to an elliptic PDE for each m-mode of the scalar field.

Compiling

Because Barry's code is written in C, we use the WSTP (formerly known as MathLink) feature of Mathematica to compile an executable with effective source functions that are callable from Mathematica code. Compiling the C code following WSTP requirements is the most challenging part of the setup. Executables have been included for x86-64 Windows and x86-64 Linux, but not for Mac (and it may be that the provided executables don't work on your system and you would need to compile them).

### If a provided executable works for you, skip forward to the marked notebook discussion

A dependency for Barry's code is the GNU Scientific Library (https://www.gnu.org/software/gsl/), install that before you proceed.

