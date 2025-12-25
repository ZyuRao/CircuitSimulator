* Example RLC-3 - Flat version
Vin 1 0 sin 0 3 200 0

* First segment (original X1)
R1 1 5 3.5
L1 5 2 1.2e-3  
C1 2 0 7.3e-6

* Second segment (original X2)
R2 2 6 3.5
L2 6 3 1.2e-3
C2 3 0 7.3e-6

* Third segment (original X3)
R3 3 7 3.5
L3 7 4 1.2e-3
C3 4 0 7.3e-6

* Load capacitor
Cl 4 0 10e-6

.HB 200 20
.TRAN 0.004e-3 20e-3
.PSS 200 1e-5
.PLOTNV 4
.PROBE tran V(4)
.PROBE pss V(4)
.print hb V(4)
.end