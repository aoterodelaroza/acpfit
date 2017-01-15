#! /bin/bash

echo "Testing fit1..."
../acpfit fit1.inp fit1.out
numdiff -q -a 1e-8 fit1.empty ref1.empty
numdiff -q -a 1e-8 fit1.acp ref1.acp
numdiff -q -a 1e-8 fit1.eval ref1.eval
echo "Testing fit2..."
../acpfit fit2.inp fit2.out
numdiff -q -a 1e-8 fit2.empty ref2.empty
numdiff -q -a 1e-8 fit2.acp ref2.acp
numdiff -q -a 1e-8 fit2.eval ref2.eval
echo "Testing fit3..."
../acpfit fit3.inp fit3.out
numdiff -q -a 1e-8 fit3.empty ref1.empty
numdiff -q -a 1e-8 fit3.acp ref1.acp
numdiff -q -a 1e-8 fit3.eval ref1.eval
echo "Testing eval1..."
../acpfit eval1.inp eval1.out
numdiff -q -a 1e-8 eval1.empty ref1.empty
numdiff -q -a 1e-8 eval1.eval ref1.eval
echo "Testing eval2..."
../acpfit eval2.inp eval2.out
numdiff -q -a 1e-8 eval2.empty ref1.empty
numdiff -q -a 1e-8 eval2.eval ref1.eval
