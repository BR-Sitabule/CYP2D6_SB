#!/bin/bas

#explicit minimisation stage 1
echo "
Test
&cntrl
imin = 1,
maxcyc = 1000,
ncyc = 500,
ntb = 1,
igb = 1
cut = 10.0
ntpr = 10,
ntr = 1,
restraint_wt=500.0,
restrainmask='1-466 & !@H='
/
" >explicit_min1.in
sander -O -i explicit_min1.in -o explicit_min1.out -p CYP2D6_min_topology.prmtop -c CYP2D6_min_topology.rst7 -ref CYP2D6_min_topology.rst7 -r relax_CYP2D6_min_explicit.rst
ambpdb -p CYP2D6_min_topology.prmtop -c relax_CYP2D6_min_explicit.rst >explicit_temporary.pdb
reduce -Trim explicit_temporary.pdb >explicit_minimized_CYP2D6.pdb

#explicit minimisation stage 2
echo "
Test
&cntrl
imin = 1,
maxcyc = 2500,
ncyc = 1000,
ntb = 1,
igb = 1
cut = 10.0
ntpr = 10,
ntr = 0
/
" >explicit_min2.in
sander -O -i explicit_min2.in -o explicit_min2.out -p CYP2D6_min_topology.prmtop -c relax_CYP2D6_min_explicit.rst -r relax_CYP2D6_min_explicit2.rst
ambpdb -p CYP2D6_min_topology.prmtop -c relax_CYP2D6_min_explicit2.rst >explicit2_temporary.pdb
reduce -Trim explicit2_temporary.pdb >explicit_minimized2_CYP2D6.pdb
