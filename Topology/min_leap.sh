#!/bin/bash

echo "
s/CYS   412/CYP   412/
" >CYP.sed 
#assigning protonation states of side chains for HIS, LYS, ASP, GLU based on their pka acquired from Propka
sed -i -f script.sed P10635.B99990010.pdb
#renaming BLK (heme and thioridazine) residues as "HEM" and "RTZ" respectively
sed -i -f HEM_RTZ.sed P10635.B99990010.pdb
#grep RTZ from model, grep one drug molecule (residue 467) and redirect it to RTZ1_file.pdb
grep "RTZ   467" P10635.B99990010.pdb > RTZ1_file.pdb
#grep RTZ from model, grep one drug molecule (residue 468) and redirect it to RTZ2_file.pdb
grep "RTZ   468" P10635.B99990010.pdb > RTZ2_file.pdb
#add hydrogen atoms
reduce RTZ1_file.pdb > reduced_RTZ1.pdb
#add hydrogen atoms
reduce RTZ2_file.pdb > reduced_RTZ2.pdb
#convert ligand file to amber compatible pdb
antechamber -i reduced_RTZ1.pdb -fi pdb  -o RTZ1.pdb -fo pdb -rn RTZ
#convert ligand file to amber compatible pdb
antechamber -i reduced_RTZ2.pdb -fi pdb  -o RTZ2.pdb -fo pdb -rn RTZ
#using the bcc method to calculate the atomic point charge and generate mol2 for drug (RTZ)
antechamber -i RTZ1.pdb -fi pdb -o RTZ.mol2 -fo mol2 -c bcc -s 2
#check if there are any missing parameters and generate frcmod file for drug (RTZ)
parmchk2 -i RTZ1.mol2 -f mol2 -o RTZ.frcmod
#Remove ligands, TER and END
grep -v "HETATM\|END\|TER" P10635.B99990010.pdb > CYP2D6_no_ligands.pdb
#add TER after last atom
grep -m 1 "TER" P10635.B99990010.pdb >> CYP2D6_no_ligands.pdb
#add TER after last RTZ atom
grep -m 1 "TER" P10635.B99990010.pdb >> RTZ2.pdb
#add TER after last RTZ atom
grep -m 1 "TER" P10635.B99990010.pdb >> RTZ1.pdb
#add RTZ molecule 1 to CYP2D6_no_ligands.pdb
cat RTZ1.pdb >> CYP2D6_no_ligands.pdb
#add RTZ molecule 2 to CYP2D6_no_ligands.pdb
cat RTZ2.pdb >> CYP2D6_no_ligands.pdb 
#add HEM to CYP2D6_no_ligands.pdb
grep " HEM   469" P10635.B99990010.pdb >> CYP2D6_no_ligands.pdb
#add TER after last HEM atom
grep -m 1 "TER" P10635.B99990010.pdb >> CYP2D6_no_ligands.pdb
#add END at the end of pdb
grep "END" P10635.B99990010.pdb >> CYP2D6_no_ligands.pdb
#final cleaning of CYP2D6 with Heme group prior running leap 
pdb4amber -i CYP2D6_no_ligands.pdb -o CYP2D6.pdb --dry
#running leap
tleap -f topology_CYP2D6.in

#in vacuo minimisation
echo "
Test 
&cntrl
imin = 1,
maxcyc = 500,
ncyc = 250,
ntb = 0,
igb = 1,
cut = 12, 
ntpr = 10, 
ntr= 1, 
restraint_wt=2.0,
restraintmask=':1-469'
/ 
" >relax.in
sander -O -i relax.in -o relax.out -p Clean_P10635.B99990010.prmtop -c Clean_P10635.B99990010.rst7 -ref Clean_P10635.B99990010.rst7 -r relax_CYP2D6_min.rst
ambpdb -p Clean_P10635.B99990010.prmtop -c relax_CYP2D6_min.rst >temporary.pdb 
reduce -Trim temporary.pdb >minimized_CYP2D6.pdb
#in vacuo minimisation

echo "
#Source leaprc file for ff14SB protein force field
source /opt/exp_soft/bioinf/amber20/dat/leap/cmd/oldff/leaprc.ff14SB
#Source leaprc file for gaff
source /opt/exp_soft/bioinf/amber20/dat/leap/cmd/leaprc.gaff
#source leaprc file for TIP3 water model
source /opt/exp_soft/bioinf/amber20/dat/leap/cmd/leaprc.water.tip3p
#Load the mol2 file for HEM
HEM =loadmol2 HEM.mol2
#Load the mol2 file for RTZ
RTZ =loadmol2 RTZ.mol2
#Load the mol2 file for CYP
CYP =loadmol2 CYP.mol2
#Load the ionsjc_tip3p frcmod file for water
loadamberparams frcmod.ionsjc_tip3p
#Load the IC6 frcmod file for the heme group
loadamberparams IC6.frcmod
#Load the RTZ frcmod file for the drug molecules
loadamberparams RTZ.frcmod
#load minimized CYP2D6 pdb
mol=loadpdb minimized_CYP2D6.pdb

#generate bond between proximal cysteine and Fe
bond mol.412.SG mol.469.FE
bond mol.412.N mol.411.C
bond mol.412.C mol.413.N

#add counterions to neutralise system and solvate box
addions mol Cl- 0
solvateoct mol TIP3PBOX 15.0

#Save AMBER topology and coordinate files
saveamberparm mol CYP2D6_min_topology.prmtop CYP2D6_min_topology.rst7
#Save pdb output file
savepdb mol CYP2D6_min_topology.pdb
#Quit tleap program
quit
" >CYP2D6_min_topology.in
 
tleap -f CYP2D6_min_topology.in
#running leap
