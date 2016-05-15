#!/bin/bash

val=0

for i in {1..50}
do
	temp=$(echo -e | ./solve)
	val=$((val + temp))
done
echo -e $((val/50))
