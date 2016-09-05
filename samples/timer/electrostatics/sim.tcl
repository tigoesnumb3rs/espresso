#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#

set backup_timer          20
set n_part                100
set int_n_times           100
set int_steps             10
set number_of_timers      22
set output_file           [open "simulation.output" "w"]
set backup_file           [open "simulation.backup" "w"]
set sim_nodes             1
set sim_procs             1

puts  $output_file "\{"
puts  $output_file "\"Nodes\":\"$sim_nodes\","
puts  $output_file "\"Processes\":\"$sim_procs\","
puts  $output_file "\"Backuprate\":\"$backup_timer\","
puts  $output_file "\"Particles\":\"$n_part\","
puts  $output_file "\"Runs\":\"$int_n_times\","
puts  $output_file "\"IntPerRun\":\"$int_steps\","
flush $output_file

puts  $backup_file "\{"
puts $backup_file "\"Nodes\":\"$sim_nodes\","
puts $backup_file "\"Processes\":\"$sim_procs\","
puts  $backup_file "\"Backuprate\":\"$backup_timer\","
puts  $backup_file "\"Particles\":\"$n_part\","
puts  $backup_file "\"Runs\":\"$int_n_times\","
puts  $backup_file "\"IntPerRun\":\"$int_steps\","
flush $backup_file


# Set system parameters
#set n_part 200
set density 0.7
set box_l [expr pow($n_part/$density,1./3.)]

# Select electrostatics method
set method "p3m"
#set method "memd"

# Setup system geometry in Espresso
setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1

# Place particles
set q 1; set type 0
for {set i 0} { $i < $n_part } {incr i} {
  set posx [expr $box_l*[t_random]]
  set posy [expr $box_l*[t_random]]
  set posz [expr $box_l*[t_random]]
  set q [expr -$q]; set type [expr 1-$type]
  part $i pos $posx $posy $posz q $q type $type
}

# Simulation parameters
setmd time_step 0.01; setmd skin 0.3
# Thermostat
set temp 1; set gamma 1
thermostat langevin $temp $gamma


# Lennard-Jones interactions
set sig 1.0; set cut [expr 1.12246*$sig]
set eps 1.0; set shift [expr 0.25*$eps]
inter 0 0 lennard-jones $eps $sig $cut $shift 0
inter 1 0 lennard-jones $eps $sig $cut $shift 0
inter 1 1 lennard-jones $eps $sig $cut $shift 0


# Check if electrostatics method is clear...
if { ![ info exists method ] } {
	puts "Please select an electrostatics method in the script."
	exit
}

# Distinguish between different methods
if { $method == "p3m" } {
	puts [inter coulomb 10.0 p3m tunev2 accuracy 1e-3]
} elseif { $method == "memd" } {
	# MEMD need no Verlet lists!
	cellsystem domain_decomposition -no_verlet_list
	set memd_mesh 12
	set f_mass [expr 100.0*pow( ([setmd time_step]*$memd_mesh/$box_l) , 2.0 ) ]
	puts "memd parameters: mesh=$memd_mesh, f_mass=$f_mass"
	puts [inter coulomb 10.0 memd $f_mass $memd_mesh]
} else {
	puts "Electrostatics method must be one of 'memd' or 'p3m'."
	exit
}


set p3m_params [inter coulomb]
foreach f $p3m_params { eval inter $f }

if { [regexp "ROTATION" [code_info]] } {
    set deg_free 6
} {
    set deg_free 3
}

set integ_steps 200
for {set cap 20} {$cap < 200} {incr cap 20} {
  puts "t=[setmd time] E=[analyze energy total]"
  inter forcecap $cap
  integrate $integ_steps
}
inter forcecap 0

timer reset

for {set i 0} { $i < $int_n_times } { incr i} {
    set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
    puts "t=[setmd time] E=[analyze energy total], T=$temp"
    integrate $int_steps

    if {$i%$backup_timer==0} {
        puts $backup_file "\"$i\": \["
        set n 1
        foreach t [timer] {
            puts $backup_file "  \["
            puts $backup_file "    \"[lindex $t 0]\","
            puts $backup_file "    \"[lindex $t 1]\","
            puts $backup_file "    \"[lindex $t 2]\","
            puts $backup_file "    \"[lindex $t 3]\","
            puts $backup_file "    \"[lindex $t 4]\","
            puts $backup_file "    \"[lindex $t 5]\","
            puts $backup_file "    \"[lindex $t 6]\","
            puts $backup_file "    \"[lindex $t 7]\","
            puts $backup_file "    \"[lindex $t 8]\""
            if {$n%[expr {$number_of_timers * $sim_procs}]} {
                puts $backup_file "  \],"
            } else {
                puts $backup_file "  \]"
            }
            incr n
        }
        puts $backup_file "  \],"
        flush $backup_file
    }



#    set f [open "config_$i" "w"]
#    blockfile $f write tclvariable {box_l density}
#    blockfile $f write variable box_l
#    blockfile $f write particles {id pos type}
#    close $f
}


puts $output_file "\"Data\": \["
set n 1
foreach t [timer] {
    puts $output_file "  \["
    puts $output_file "    \"[lindex $t 0]\","
    puts $output_file "    \"[lindex $t 1]\","
    puts $output_file "    \"[lindex $t 2]\","
    puts $output_file "    \"[lindex $t 3]\","
    puts $output_file "    \"[lindex $t 4]\","
    puts $output_file "    \"[lindex $t 5]\","
    puts $output_file "    \"[lindex $t 6]\","
    puts $output_file "    \"[lindex $t 7]\""
    if {$n%[expr {$number_of_timers * $sim_procs}]} {
        puts $output_file "  \],"
    } else {
        puts $output_file "  \]"
    }
    incr n
}
puts $output_file "\],"


puts $output_file "\"RunsTaken\":\"$i\","
puts $output_file "\"State\":\"SUCCESS\""
puts $output_file "\}"
flush $output_file

puts $backup_file "\}"
flush $backup_file



puts "MPIPID FUNCTIONNAME MEAN STANDARD-DEVIATION VARIANCE MIN MAX #CALLS"
foreach t [timer] {
    puts $t
}
