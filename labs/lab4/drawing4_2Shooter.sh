#!/bin/bash

gnuplot << EOP

set terminal png size 1024, 768
set output 'graphic4_2Shooter.png'
set xlabel "X"
set ylabel "F(x)"
set termoption dash
plot 'plotDataShooterRR' u 1:2 title "Shooter" w lines, 'plotDataShooterRR' u 1:3 title "Theoretical" w lines, 'plotDataShooterRR' u 1:3 title "RungeRomberg" w lines

EOP