#!/bin/bash

cat > pdos.gp <<EOF
set terminal pdf color enhanced \
dashed dl 0.5 size 8.0cm, 6.0cm
set output "pdos.pdf"

set border lw 2

set ytics scale 3.0, -0.5 font 'Cmr10,18'
set xtics scale 3.0, -0.5 font 'Cmr10,18'

set yzeroaxis

set ylabel "Partial DOS [eV^{-1}]" offset 0.0, 0.0 font 'Cmr10,18'
set xlabel "Energy from {/Cmmi10 \042}_F [eV]" offset 0.0, 0.0 font 'Cmr10,18'
set key left top font 'Cmr10,18
plot [][0:] \\
EOF

for i in ${@: 2:$#-2}
do
    echo "\"$i\" u (\$1-$1):2 w l tit \"$i\", \\" >> pdos.gp
done
echo "\"${@: $#}\" u (\$1-$1):2 w l tit \"${@: $#}\"" >> pdos.gp
