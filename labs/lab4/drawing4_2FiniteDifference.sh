#!/bin/bash

gnuplot << EOP

set terminal png size 1024, 768
set output 'graphic4_2FiniteDifference.png'
set xlabel "X"
set ylabel "F(x)"
set termoption dash
plot 'plotDataFiniteDifferenceRR' u 1:2 title "FiniteDifference" w lines, 'plotDataFiniteDifferenceRR' u 1:3 title "Theoretical" w lines, 'plotDataFiniteDifferenceRR' u 1:3 title "RungeRomberg" w lines

EOP