Gauss Computer Programs Chapter 3,
2nd edition Heer/Maussner, Dynamic General Equilibrium Modeling
===============================================================

The files

Ramsey2c.g          
Ramsey3c.g          
SOE.g               
zvec_Fig3_3.fmt              

store the Gauss programs referred to in Chapter 3 of our book. The vector
zvec_Fig3_3.fmt is necessary to solve Problem 3.1.

To run the programs you need the procedure FixvMN2 and other procedures called by FixvMN2. They
are stored in the files NLEQ.src and Derivatives.src (which are in the zip archive of Chapter 11) as well as GraphSettings from
Tools.src (also in the zip archive of Chapter 11).

You can either copy the required code into the respective main file,
use #include statements or add the required procedures to your user.lcg library.

The programs run under Gauss 7.0 and higher.

In case you have any questions please contact me at 

alfred.maussner@wiwi.uni-augsburg.de
