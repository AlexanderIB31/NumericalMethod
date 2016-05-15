#!/bin/bash

gnuplot << EOP

set terminal png size 1024, 768
set output 'graphic3_2.png'
plot 'plotData0' u 1:2 w lines, 'plotData1' u 1:2 w lines, 'plotData2' u 1:2 w lines, 'plotData3' u 1:2 w lines

EOP