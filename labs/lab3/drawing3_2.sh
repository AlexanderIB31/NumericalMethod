#!/bin/bash

gnuplot << EOP

set terminal png size 1024, 768
set output 'graphic3_2.png'
plot 'plotData0' u 1:2 w lines lw 3, 'plotData1' u 1:2 w lines lw 3, 'plotData2' u 1:2 w lines lw 3, 'plotData3' u 1:2 w lines lw 3

EOP