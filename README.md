# EllipticSSF

Authors: Nami Nishimura and Thomas Osburn_ 

Development guidelines
1. We should start with only one file in the directory, "EllipticSSF.wl"
2. The Mathematica package serves only to define Modules that are called elsewhere 
3. Eliminate all global variables. Add extra function parameters if necessary.
4. Avoid ".nb" files in the git directory until we have a mature example notebook (Mathematica notebooks are not well suited for git).

Typical git procedure
1. BEFORE making changes, pull the most recent version: 
git pull origin main
2. After you reach a stopping point (small change where the code is functional), push them:
git add .
git commit -m "Type an informative message about your changes here"
git push origin main

