set GMtime [expr $dt*$TotalNumberOfSteps];

#Define acceleration series
set accelSeries "Series -dt $dt -filePath $GMfile -factor [expr $Scalefact]";
 
#Create load pattern
# define where and how (pattern tag, dof) acceleration is applied
pattern UniformExcitation 2 1 -accel $accelSeries;
 
#Dynamic analysis
#set dt_analysis 0.0001;
set dt_analysis $dt;
wipeAnalysis;
constraints Plain;
#how to number the dof's Reverse Cuthill-McKee (alternative: Plain)
numberer RCM;
#Generate linear system of equations: General sparse system of equations
system UmfPack;
#convergence test #limit #max iterations
test NormDispIncr 1.0e-6 1000;
algorithm Newton;
#newmark integration #gamma and #beta parameters 
integrator Newmark 0.5 0.25;
#dynamic
analysis Transient;
set NumSteps [expr round(($GMtime + 0.0)/$dt_analysis)];
analyze $NumSteps $dt_analysis;

#Output time at end of analysis	
set currentTime [getTime];
puts "The current time is: $currentTime";

wipe;
