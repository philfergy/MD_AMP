
#####ENTER STANDARDIZED TRAJECTORY#######
###############NAME######################

#enter trajectory name
read -p "Enter name of trajectory (.xtc) file without extension  " TRAJ

#enter topology name
read -p "Enter name of toplogy (.tpr) file without extension  " TOPOL

#enter Number for membrane type
read -p "1. PEPG 2.PG " LIPIDS

#enter Peptide name
read -p "Enter peptide name " FILENAME


dir=$(pwd)

#enter residue range for calculation (with -)
read -p "enter residue start" START
read -p "enter residue end" END

######FIND NUMBER OF INDEX GROUPS########
#########################################


echo "0" | gmx trjconv -f $dir/traj_comp.xtc -s $dir/md_0_1.tpr -b 0 -e 0 -o $dir/fullsys.gro 


echo "q" | gmx make_ndx -f $dir/fullsys.gro > $dir/temp1 
grep atoms $dir/temp1 > $dir/temp2
tot_lines=$(wc -l ${dir}/temp2 | awk '{print $1}')
group_num_low=$((tot_lines-1))
group_num_high=$((tot_lines))


echo "#########$tot_lines################"

echo "#########$group_num_high################"






###output files

if [ $LIPIDS == 1 ]
then

########CREATE INDEX FILES###############
#########################################



####                       MAKE IN4
echo "resname POPG and not same residue as within 0.4 of resid $START to $END" | gmx select -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr  -n $dir/index.ndx -on $dir/PG_notin4.ndx -b 200000 -e 200000
cat $dir/PG_notin4.ndx $dir/index.ndx > $dir/temp.ndx
##make sn1 index file
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPG_sn1_notin4.ndx -n $dir/temp.ndx << EOF
0 & a C31
0 & a C32
0 & a C33
0 & a C34
0 & a C35
0 & a C36
0 & a C37
0 & a C38
0 & a C39
0 & a C310
0 & a C311
0 & a C312
0 & a C313
0 & a C314
0 & a C315
0 & a C316
del 0-${group_num_high}
q
EOF
echo "resname POPE and not same residue as within 0.4 of resid $START to $END" | gmx select -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -n $dir/index.ndx -on $dir/PE_notin4.ndx -b 200000 -e 200000
cat $dir/PE_notin4.ndx $dir/index.ndx > $dir/temp.ndx
##make sn1 index file
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPE_sn1_notin4.ndx -n $dir/temp.ndx << EOF
0 & a C31
0 & a C32
0 & a C33
0 & a C34
0 & a C35
0 & a C36
0 & a C37
0 & a C38
0 & a C39
0 & a C310
0 & a C311
0 & a C312
0 & a C313
0 & a C314
0 & a C315
0 & a C316
del 0-${group_num_high}
q
EOF


####                       MAKE NOT IN4

echo "resname POPG and same residue as within 0.4 of resid $START to $END" | gmx select -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -n $dir/index.ndx -on $dir/PG_in4.ndx -b 200000 -e 200000
cat $dir/PG_in4.ndx $dir/index.ndx > $dir/temp.ndx
##make sn1 index file
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPG_sn1_in4.ndx -n $dir/temp.ndx << EOF
0 & a C31
0 & a C32
0 & a C33
0 & a C34
0 & a C35
0 & a C36
0 & a C37
0 & a C38
0 & a C39
0 & a C310
0 & a C311
0 & a C312
0 & a C313
0 & a C314
0 & a C315
0 & a C316
del 0-${group_num_high}
q
EOF
echo "resname POPE and same residue as within 0.4 of resid $START to $END" | gmx select -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -n $dir/index.ndx -on $dir/PE_in4.ndx -b 200000 -e 200000
cat $dir/PE_in4.ndx index.ndx > $dir/temp.ndx
##make sn1 index file
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPE_sn1_in4.ndx -n $dir/temp.ndx << EOF
0 & a C31
0 & a C32
0 & a C33
0 & a C34
0 & a C35
0 & a C36
0 & a C37
0 & a C38
0 & a C39
0 & a C310
0 & a C311
0 & a C312
0 & a C313
0 & a C314
0 & a C315
0 & a C316
del 0-${group_num_high}
q
EOF


###     MAKE WHOLE
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPG_sn1.ndx -n $dir/index.ndx << EOF
13 & a C31
13 & a C32
13 & a C33
13 & a C34
13 & a C35
13 & a C36
13 & a C37
13 & a C38
13 & a C39
13 & a C310
13 & a C311
13 & a C312
13 & a C313
13 & a C314
13 & a C315
13 & a C316
del 0-${group_num_low}
q
EOF
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPE_sn1.ndx -n $dir/index.ndx << EOF
14 & a C31
14 & a C32
14 & a C33
14 & a C34
14 & a C35
14 & a C36
14 & a C37
14 & a C38
14 & a C39
14 & a C310
14 & a C311
14 & a C312
14 & a C313
14 & a C314
14 & a C315
14 & a C316
del 0-${group_num_low}
q
EOF


########RUN ANALYSIS#####################
#########################################

gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPG_sn1.ndx -n $dir/POPG_sn1.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPG.xvg
rm order.xvg
gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPE_sn1.ndx -n $dir/POPE_sn1.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPE.xvg
rm order.xvg
gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPG_sn1_in4.ndx -n $dir/POPG_sn1_in4.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPG_in4.xvg
rm order.xvg
gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPE_sn1_in4.ndx -n $dir/POPE_sn1_in4.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPE_in4.xvg
rm order.xvg
gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPG_sn1_notin4.ndx -n $dir/POPG_sn1_notin4.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPG_notin4.xvg
rm order.xvg
gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPE_sn1_notin4.ndx -n $dir/POPE_sn1_notin4.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPE_notin4.xvg
fi

if [ $LIPIDS == 2 ]
then

####                       MAKE IN4
echo "resname POPG and not same residue as within 0.4 of resid $START to $END" | gmx select -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -n $dir/index.ndx -on $dir/PG_notin4.ndx -b 200000 -e 200000
cat $dir/PG_notin4.ndx $dir/index.ndx > $dir/temp.ndx
##make sn1 index file
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPG_sn1_notin4.ndx -n $dir/temp.ndx << EOF
0 & a C31
0 & a C32
0 & a C33
0 & a C34
0 & a C35
0 & a C36
0 & a C37
0 & a C38
0 & a C39
0 & a C310
0 & a C311
0 & a C312
0 & a C313
0 & a C314
0 & a C315
0 & a C316
del 0-${group_num_high}
q
EOF



####                       MAKE NOT IN4

echo "resname POPG and same residue as within 0.4 of resid $START to $END" | gmx select -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -n $dir/index.ndx -on $dir/PG_in4.ndx -b 200000 -e 200000
cat $dir/PG_in4.ndx $dir/index.ndx > $dir/temp.ndx
##make sn1 index file
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPG_sn1_in4.ndx -n $dir/temp.ndx << EOF
0 & a C31
0 & a C32
0 & a C33
0 & a C34
0 & a C35
0 & a C36
0 & a C37
0 & a C38
0 & a C39
0 & a C310
0 & a C311
0 & a C312
0 & a C313
0 & a C314
0 & a C315
0 & a C316
del 0-${group_num_high}
q
EOF



###     MAKE WHOLE
gmx make_ndx -f $dir/fullsys.gro -o $dir/POPG_sn1.ndx << EOF
13 & a C31
13 & a C32
13 & a C33
13 & a C34
13 & a C35
13 & a C36
13 & a C37
13 & a C38
13 & a C39
13 & a C310
13 & a C311
13 & a C312
13 & a C313
13 & a C314
13 & a C315
13 & a C316
del 0-${group_num_low}
q
EOF

########RUN ANALYSIS#####################
#########################################

gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPG_sn1.ndx -n $dir/POPG_sn1.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPG.xvg
gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPG_sn1_in4.ndx -n $dir/POPG_sn1_in4.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPG_in4.xvg
gmx order -f $dir/$TRAJ.xtc -s $dir/$TOPOL.tpr -nr $dir/POPG_sn1_notin4.ndx -n $dir/POPG_sn1_notin4.ndx -b 100000 -e 200000 -xvg none -od $dir/${FILENAME}_POPG_notin4.xvg

fi

