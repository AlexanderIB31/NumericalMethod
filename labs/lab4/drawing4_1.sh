#!/bin/bash

gnuplot << EOP

set terminal png size 1024, 768
set output 'graphic4_1.png'
set xlabel "X"
set ylabel "F(x)"
set termoption dash
plot 'plotData4_1' u 1:2 title "Runge-Kutta" w lines, 'plotData4_1' u 1:3 title "Euler" w lines, 'plotData4_1' u 1:3 title "Adams" w lines

EOP