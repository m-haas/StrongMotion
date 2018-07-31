wipe;
source ESDOFparams.tcl;
source ESDOFgravity.tcl;
recorder Drift -file Drift.out -time -iNode 1 -jNode 2 -dof 1 -perpDirn 2;
set GMfile "gm.acc";
set dt 0.005000;
set TotalNumberOfSteps 8768;
set Scalefact 1.0;
source ESDOFdynamicals.tcl;
