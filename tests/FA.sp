* ==========================================================
* Full Adder (CMOS NAND-only), MOSFET + Capacitor loads
* Nodes:
*   VDD=103, GND=0
*   A=101, B=102, CIN=105
*   SUM=118, COUT=117
*
* Truth table (0/3V):
*   A B Cin | Sum Cout
*   0 0 0   |  0   0
*   0 0 3   |  3   0
*   0 3 0   |  3   0
*   0 3 3   |  0   3
*   3 0 0   |  3   0
*   3 0 3   |  0   3
*   3 3 0   |  0   3
*   3 3 3   |  3   3
* ==========================================================

* -------- supply --------
VDD   103 0  DC 3

* -------- inputs (DC truth table mode: run .op 8 times, edit these 3 lines) --------
VA    101 0  DC 3
VB    102 0  DC 3
VCIN  105 0  DC 0

* -------- (OPTIONAL) inputs as square waves (needs PULSE support in your parser) ----
* 8-combo sequence over 80ns:
* A toggles fastest, then B, then CIN
* VA   101 0  PULSE 0 3  1e-9 0.1e-9 0.1e-9 10e-9 20e-9
* VB   102 0  PULSE 0 3  1e-9 0.1e-9 0.1e-9 20e-9 40e-9
* VCIN 105 0  PULSE 0 3  1e-9 0.1e-9 0.1e-9 40e-9 80e-9

* -------- transistor sizes (buffer.sp style) --------
* PMOS: Wp=60e-6, NMOS: Wn=20e-6, L=0.35e-6

* ==========================================================
* XOR1: XAB = A XOR B via 4 NAND
* N1  = NAND(A,B)
* N2  = NAND(A,N1)
* N3  = NAND(B,N1)
* XAB = NAND(N2,N3)
* ==========================================================

* N1 = ~(A&B) -> node 110, mid 109
MPN1A 110 101 103 p 60e-6 0.35e-6 1
MPN1B 110 102 103 p 60e-6 0.35e-6 1
MNN1A 110 101 109 n 20e-6 0.35e-6 2
MNN1B 109 102 0   n 20e-6 0.35e-6 2

* N2 = ~(A&N1) -> node 112, mid 111
MPN2A 112 101 103 p 60e-6 0.35e-6 1
MPN2B 112 110 103 p 60e-6 0.35e-6 1
MNN2A 112 101 111 n 20e-6 0.35e-6 2
MNN2B 111 110 0   n 20e-6 0.35e-6 2

* N3 = ~(B&N1) -> node 114, mid 113
MPN3A 114 102 103 p 60e-6 0.35e-6 1
MPN3B 114 110 103 p 60e-6 0.35e-6 1
MNN3A 114 102 113 n 20e-6 0.35e-6 2
MNN3B 113 110 0   n 20e-6 0.35e-6 2

* XAB = ~(N2&N3) -> node 116, mid 115
MPX1  116 112 103 p 60e-6 0.35e-6 1
MPX2  116 114 103 p 60e-6 0.35e-6 1
MNX1  116 112 115 n 20e-6 0.35e-6 2
MNX2  115 114 0   n 20e-6 0.35e-6 2

* ==========================================================
* XOR2: SUM = XAB XOR Cin via 4 NAND
* M1  = NAND(XAB,Cin)          (reuse for carry)
* M2  = NAND(XAB,M1)
* M3  = NAND(Cin,M1)
* SUM = NAND(M2,M3)
* ==========================================================

* M1 = ~(XAB&Cin) -> node 120, mid 121
MPM1A 120 116 103 p 60e-6 0.35e-6 1
MPM1B 120 105 103 p 60e-6 0.35e-6 1
MNM1A 120 116 121 n 20e-6 0.35e-6 2
MNM1B 121 105 0   n 20e-6 0.35e-6 2

* M2 = ~(XAB&M1) -> node 122, mid 123
MPM2A 122 116 103 p 60e-6 0.35e-6 1
MPM2B 122 120 103 p 60e-6 0.35e-6 1
MNM2A 122 116 123 n 20e-6 0.35e-6 2
MNM2B 123 120 0   n 20e-6 0.35e-6 2

* M3 = ~(Cin&M1) -> node 124, mid 125
MPM3A 124 105 103 p 60e-6 0.35e-6 1
MPM3B 124 120 103 p 60e-6 0.35e-6 1
MNM3A 124 105 125 n 20e-6 0.35e-6 2
MNM3B 125 120 0   n 20e-6 0.35e-6 2

* SUM = ~(M2&M3) -> node 118, mid 126
MPS1  118 122 103 p 60e-6 0.35e-6 1
MPS2  118 124 103 p 60e-6 0.35e-6 1
MNS1  118 122 126 n 20e-6 0.35e-6 2
MNS2  126 124 0   n 20e-6 0.35e-6 2

* ==========================================================
* Carry:
* COUT = (A&B) OR (Cin&XAB) = NAND( N1 , M1 )
* N1=110, M1=120
* ==========================================================

* COUT = ~(N1&M1) -> node 117, mid 127
MPC1  117 110 103 p 60e-6 0.35e-6 1
MPC2  117 120 103 p 60e-6 0.35e-6 1
MNC1  117 110 127 n 20e-6 0.35e-6 2
MNC2  127 120 0   n 20e-6 0.35e-6 2

* -------- load capacitors (only caps, no resistors) --------
CXAB   116 0 0.05e-12
CSUM   118 0 0.10e-12
CCOUT  117 0 0.10e-12

* -------- models (copy buffer.sp style) --------
.MODEL 1 VT -0.75 MU 5e-2  COX 0.3e-4 LAMBDA 0.05 CJ0 4.0e-14
.MODEL 2 VT 0.83  MU 1.5e-1 COX 0.3e-4 LAMBDA 0.05 CJ0 4.0e-14

* -------- outputs / controls --------
* .PLOTNV: outputs these node voltages for ALL analyses to CSV
.PLOTNV 101
.PLOTNV 102
.PLOTNV 105
.PLOTNV 117
.PLOTNV 118

* If you implemented your new semantics:
* .PRINT hb/tran: put selected into table(csv)
* .PROBE hb/tran: plot selected
.PRINT tran V(117) V(118)
.PRINT hb   V(116) V(117) V(118)

.PROBE tran V(117) V(118)
.PROBE hb   V(116) V(117) V(118)

* -------- analyses --------
* .op
* (optional transient for pulse mode)
.TRAN 0.2e-9 80e-9
* (optional hb for pulse periodic 80ns -> f0=12.5MHz)
* .hb 12.5e6 30
