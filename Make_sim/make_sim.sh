echo $'###################################################\n
#############MEMBRANE PARAMS#######################\n
###################################################\n'

BOX_X=0
read -p $'What are your membrane x,y dimensions?\n' BOX_X

BOX_Z=0
read -p $'What are your membrane z dimensions?\n' BOX_Z

MEMB_FILE=" "
read -p $'What is the name of your membrane\n' MEMB_FILE

LIPIDS=0
read -p $'What is the lipid composition of your bilayer? \n 1. POPG \n 2. POPC \n 3. POPE:POPG(3:1) \n' LIPIDS

LIPNUM=0
read -p $'How many lipids in bilayer?\n' LIPNUM


PEP_FILE=" "
read -p $'What is the name of your peptide pdb file?\n' PEP_FILE



PEPNUM=0
read -p $'How many peptides in simulation? \n' PEPNUM


###################################################
#############INSERT PEPTIDES#######################
###################################################

#Edit vdwradii
cp radii_init.dat vdwradii.dat
sed -i '10s/.*/???  C     1.0/' vdwradii.dat
#Insert molecules -  edit box for different bilayer
gmx insert-molecules -ci $PEP_FILE -o peptides.pdb -nmol $PEPNUM -box $BOX_X $BOX_X 1.0 -try 1000
#seperate chains
tail -n +5 peptides.pdb | head -n -2 > ter.pdb
#Insert ter every n lines dependign on petpide length
PEPLEN=$(wc -l ter.pdb | awk '{ print $1 }')
RESLEN=$((PEPLEN/PEPNUM))
sed "0~${RESLEN} s/$/\nTER/g" < ter.pdb > ter_2.pdb


###################################################
###################################################
###################################################

echo '$###################################################\n
#############PROCESS PDB FILE######################\n
###################################################\n'

HISPRO=" "
HISCHOICE=0
read -p $'Are histidine residues protonated? \n 1. YES \n 2. NO \n' HISCHOICE
if [ $HISCHOICE == 1 ]
then
HISPRO="-his"
fi

ASPCHOICE=0
ASPRO=" "
read -p $'Are Aspartate residues protonated? \n 1. YES \n 2. NO \n' ASPCHOICE
if [ $ASPCHOICE == 1 ]
then
ASPRO="-asp"
fi

GLUCHOICE=0
GLUPRO=" "
read -p $'Are Glutamate residues protonated? \n 1. YES \n 2. NO \n' GLUCHOICE
if [ $GLUCHOICE == 1 ]
then
GLUPRO="-glu"
fi

gmx pdb2gmx -ff charmm36 -f ter_2.pdb -o pro.pdb -chainsep ter -ter $HISPRO $ASPRO $GLUPRO

###################################################
###################################################
###################################################
dist=8.8
box_prep_z=$(bc <<< "$BOX_Z+$dist")
gmx editconf -f pro.pdb -c -box $BOX_X $BOX_X $box_prep_z -o newbox.pdb
head -n -2 newbox.pdb > temp_1.pdb
sed '/^\(ATOM\|-t\(h\|o\)\)/!d' $MEMB_FILE > temp_3.pdb
cat temp_1.pdb temp_3.pdb > system.pdb
gmx editconf -f system.pdb -o system_newbox.pdb -c -box $BOX_X $BOX_X 14

###################################################
############'Solvate and add ions##################
###################################################

###WORK ON THE BUILDING TOPOLOGY PART!!!
cp base_topology.top topol.top

if [ $LIPIDS == 1 ]
then
echo "POPG      $LIPNUM" >> topol.top
fi

if [ $LIPIDS == 2 ]
then
echo "POPC      $LIPNUM" >> topol.top
fi

if [ $LIPIDS == 3 ]
then
PE_MIX=$(((LIPNUM/4)*3))
PG_MIX=$((LIPNUM/4))
echo "POPE      $PE_MIX" >> topol.top
echo "POPG      $PG_MIX" >> topol.top

fi




sed -i '10s/.*/???  C     0.4/' vdwradii.dat
gmx solvate -cp system_newbox.pdb -p topol.top -o system_solv.pdb
gmx grompp -f ions.mdp -c system_solv.pdb -p topol.top -o ions.tpr -maxwarn 1
gmx genion -s ions.tpr -o system_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
rm vdwradii.dat

