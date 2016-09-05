set monomer_density 0.85
set n_polymers  10
set n_beads 15
set polymer_bond_length 1.0
set monomer_type 0
set kT 1.0
set fene_cut 2.0
set fene_k 7.0

set box [expr (($n_polymers * $n_beads) / $monomer_density)**(1./3.)]

puts "box $box"

setmd time_step 0.01
setmd skin 0.4
setmd box_l $box $box $box
thermostat langevin $kT 0.5

inter 0 fene $fene_k $fene_cut

inter 0 0 lennard-jones 1.0 1.0 [expr 2.**(1./6.)] auto

puts [inter]

polymer $n_polymers $n_beads $polymer_bond_length types $monomer_type bond 0

puts [analyze energy]

set energy_fd [open "energy.dat" "w"]

for { set i 0  }  { $i < 500 } { incr i } {
inter forcecap [expr 5.0 +$i * 0.5]
integrate 1000
puts [analyze energy]

puts $energy_fd "[setmd time] [analyze energy total]"
flush $energy_fd
}

puts [timer]

inter forcecap 0

for { set i 0 } { $i < 1000 } { incr i } {
integrate 1000

#correlation $cor_sf autoupdate
# set c_times [correlation $cor_sf print correlation_time]
# 
# set max 0.0
# foreach c_time $c_times {
#   if { $c_time > $max } {
# 	  set max $c_time
#   }
# }
# 
# puts "c_time $max"

}

close $energy_fd

puts [timer]
