#!/bin/bash

gnuplot << EOP

set terminal png size 1024, 768
set output 'graphic.png'
set xlabel "Size of matrix"
set ylabel "Time (ms)"
set termoption dash
plot 'plotDataGPU' u 1:2 title "GPU" w lines lc 1, 'plotDataCPU' u 1:2 title "CPU" w lines lc 2

EOP