#!/usr/bin/tclsh
mol new ./visual.gro waitfor all

set sel [atomselect 0 all]
set numframes [molinfo 0 get numframes]
set fp [open "./visual_index.txt" r]

for {set i 0} {$i < $numframes} {incr i} {
        gets $fp line
        $sel frame $i
        puts "Setting color for frame $i ..."
        $sel set user $line
}
close $fp
$sel delete


#H2O
mol modstyle 0 0 DynamicBonds 3.6 0.1
mol modselect 0 0 "resname H2O and user != 10.000000"
mol modcolor 0 0 User
mol addrep 0
mol modstyle 1 0 VDW 0.1
mol modselect 1 0 "resname H2O and user = 10.000000"
mol addrep 0
mol modstyle 2 0 VDW 0.2
mol modselect 2 0 "resname H2O and user != 10.000000"
mol modcolor 2 0 User
mol addrep 0
#additive guest
mol modstyle 3 0 DynamicBonds 1.4 0.1
mol modselect 3 0 "resname add"
mol addrep 0
mol modstyle 4 0 VDW 0.2
mol modselect 4 0 "resname add guest"

display resetview
rotate y by 90
display projection orthographic
display depthcue 0
pbc box
