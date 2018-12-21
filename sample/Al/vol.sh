#!/bin/bash

for a in 1.90 1.95 2.00 2.05 2.10
do
    cat > scf${a}.in <<EOF
&CONTROL
 calculation = 'scf'
/
&SYSTEM
nbnd = 5
         nat = 1
        ntyp = 1
     ecutwfc =30.000000
     ecutrho = 120.000000
/
&ELECTRONS
 mixing_beta = 0.3
conv_thr = 1.000000e-5
electron_maxstep = 100
/
CELL_PARAMETERS
 0.00 ${a} ${a}
 ${a} 0.00 ${a}
 ${a} ${a} 0.00
ATOMIC_SPECIES
 Al al.lda.lps
ATOMIC_POSITIONS 
 Al 0.000000 0.000000 0.000000
K_POINTS
 8 8 8
EOF
    ../../src/pwdft.x < scf${a}.in > scf${a}.out
done
