#!/bin/bash

gnuplot << EOP

set terminal png size 1024, 768
set output 'graphic3_3.png'
set xlabel "X"
set ylabel "F(x)"
set termoption dash
plot 'plotDataPoint' u 1:2 title "Point-Value" w points pt 15 ps 2, 'plotDataApprox1' u 1:2 title "1st power" w lines, 'plotDataApprox2' u 1:2 title "2nd power" w lines

EOP