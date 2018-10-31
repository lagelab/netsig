#!/bin/bash

set echo
set verbose

#run benchmark
#max=1000000;
#for ((i=1; i<=max; i*=10 )); do
#    echo -ne "$i\n";
#  ./runNetwork -p input/pvals -j 1 -J 1 -n input/network0.bin -N test -i $((i)) -P 3 -g -o . -M input/InWeb3.HUGO.mapping -V
#done

rm benchmark.log

max=100000;
#max=100000000;
for ((i=1; i<=max; i*=10 )); do
  for ((j=1; j<=5; j++)); do
    echo -ne "\r$i $j";
    ./runNetwork -p input/pvals -j 1 -J 1 -b -n input/network0.bin -N test${j} -i $((i)) -P 3 -g -o . -M input/InWeb3.HUGO.mapping -e >> benchmark.log &
  done
  wait
  echo -ne "\n";
done

#max=10
#for ((i=1; i<=max; i++ )); do
#  echo -ne "$i\n";
#  ./runNetwork -p input/pvals -j 1 -J 1 -n input/network0.bin -N test${i} -i 1000000 -P 3 -g -o . -M input/InWeb3.HUGO.mapping &
#  sleep 5
#done
