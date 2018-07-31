#parameters from ESDOF
set g 9.81;                                   # gravitational acceleration m/s2 
set Aw 0.174400400708;                        # total wall area m2
set H 2.7;                                    # height m
set Mg 91153.0644091;                         # generalized mass kg
set Iplan 0.0025346249806;                    # inertia m4
set Kg 3090534.95577;                         # generalized stiffness N/m
set Eg 8000000000.0;                          # elasticity modulus Pa 
set Wg [expr {-1.*$Mg*$g}];                   # generalized Weight
set DR 0.05;                                  # Damping ratio for material
