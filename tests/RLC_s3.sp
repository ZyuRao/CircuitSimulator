* Example RLC-3 - Flat version
Vin 1 0 sin 0 3 200 0

* First segment (original X1)
R1 1 5 3.5
L1 5 2 1.2m  
C1 2 0 7.3u

* Second segment (original X2)
R2 2 6 3.5
L2 6 3 1.2m
C2 3 0 7.3u

* Third segment (original X3)
R3 3 7 3.5
L3 7 4 1.2m
C3 4 0 7.3u

* Load capacitor
Cl 4 0 10u

.HB 200 20
.TRAN 0.004m 20m
.PLOTNV 4
.PROBE tran V(4)
.print hb V(4)
.end