#----------Create Model (2-D 3dof problem)-----------
#----------SAC 3 Storied Building, LA Down Town--------
######################
model BasicBuilder -ndm 2 -ndf 3;

#Unit N,m,pa,kg 

# Create nodes
node 1 0.0 0.0;
node 2 0.0 $H;

# Fix supports at node
fix 1 1 1 1;
mass 2 $Mg 0. 0.;			   # node#, Mx My Mz, Mass=Weight/g.

# Coordinate transformation
#geomTransf PDelta 1;
geomTransf Linear 1;

#Column Elements
# element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag
element elasticBeamColumn 1 1 2 $Aw $Eg $Iplan 1;

#Constant gravity loads
pattern Plain 1 Constant {
load 2 0.0 $Wg 0.0;
}

# Gravity-analysis
constraints Plain;
#constraints Transformation;
numberer Plain;
system BandGeneral;
#system UmfPack;
#test the normal displacement increments #maximum tolerance #maximum iterations
test NormDispIncr 1.0e-6 100;
#Newton-Raphson approach for linearization
algorithm Newton;
#number of steps
set NstepGravity 10;
set DGravity [expr 1.0/$NstepGravity];
integrator LoadControl $DGravity;
#static vs transient: here use static
analysis Static;
analyze $NstepGravity;

# maintain constant gravity loads and reset time to zero
loadConst -time 0.0;

# set damping based on first eigen mode
set freq [expr [eigen -fullGenLapack 1]**0.5]
set dampRatio $DR
rayleigh 0. 0. 0. [expr 2*$dampRatio/$freq]
