#!/bin/sh

for i in m12_gas_new*
do
 suff="${i##*new}"
 mv $i m12_mr$suff
done

