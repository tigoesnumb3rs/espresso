require_feature EXTERNAL_FORCES
require_feature COLLISION_DETECTION
require_feature LENNARD_JONES
require_feature BOND_ANGLE

# Get environment variable or default
proc getenv {envvar defaultval} {
    if {[info exists ::env($envvar)]} {
        return $::env($envvar)
    } else {
        return $defaultval
    }
}

proc is_envset {envvar} {
    if {[info exists ::env($envvar)]} {
        return "(env)"
    } else {
        return ""
    }
}

# Get the hostname as hostname(1) outputs it. [info hostname] somehow includes
# domain information.
proc get_hostname {} {
    set f [open "/etc/hostname"]
    gets $f hn
    close $f
    return $hn
}

# Dump all simulation parameters
proc dump_config {} {
    global outdir
    set f [open "$outdir/CONFIGURATION" "w"]
    blockfile $f write tclvariable {rcut eps sigma temp gamma n_part box_l skin dt n_cycle n_steps ckp_output vmd_output outdir checkpoint ckp_type hostnam capincr cap warmup_us nwarmupiter mindist bind_radius}
    blockfile $f write variable all
    blockfile $f write interactions
    blockfile $f write thermostat
    close $f
}

# Formats microseconds to seconds with two fractional digits and appended " s"
proc format_sec {us} {
    return [format "%.2f s" [expr $us / 1e6]]
}

# Write a blockfile
proc write_checkpoint_blockfile {iter} {
    global outdir checkpoint
    set f [open "|gzip -c - > $outdir/$checkpoint$iter.dat.gz" "w"]
    blockfile $f write random
    blockfile $f write seed
    blockfile $f write particles {id pos v f}
    blockfile $f write bonds all
    close $f
}

# Dump via MPI-IO
proc write_checkpoint_mpiio {iter} {
    global outdir checkpoint
    set f [open "$outdir/$checkpoint$iter.blockfile" "w"]
    blockfile $f write random
    blockfile $f write seed
    close $f

    mpiio "$outdir/$checkpoint$iter" write type pos v bond
}

# Checkpoint dispatching function
# Used method depends on global ckp_type variable
# Aborts the simulation if ckp_type is not meaningful
proc write_checkpoint {iter} {
    global ckp_type
    set ckpnt_us [lindex [time {
        if { $ckp_type == "blockfile" } {
            write_checkpoint_blockfile $iter
        } elseif { $ckp_type == "mpiio" } {
            write_checkpoint_mpiio $iter
        } else {
            puts stderr "Unknown checkpoint type: $ckp_type"
            puts stderr "Exiting."
            exit 1
        }
    }] 0]
    puts "\[Checkpointing\] Wrote iteration $iter (Took [format_sec $ckpnt_us])"
}

## Variables section ##

# Particles
set n_part [getenv SIM_NPART 10000]

# Box
set box_l [getenv SIM_BOXL 640]

# Simulation time
set dt [getenv SIM_DT 0.01]
set n_cycle [getenv SIM_NCYCLE 10000]
set n_steps [getenv SIM_NSTEPS 1000]


# Thermostat parameters
set skin [getenv SIM_SKIN 0.4]
set temp [getenv SIM_TEMP 1.0]
set gamma [getenv SIM_GAMMA 1.0]

# Lennard-Jones parameters
set eps [getenv SIM_EPS 1.0]
set sigma [getenv SIM_SIGMA 1.0]
set rcut [getenv SIM_RCUT 2.5]

# Bonding
set bind_radius [getenv SIM_BIND_RADIUS 1.1]

# Checkpointing and visualization options
set vmd_output [getenv SIM_VMDOUTPUT "yes"]
set ckp_output [getenv SIM_CKPOUTPUT "no"]
set ckp_type [getenv SIM_CKPTYPE "mpiio"]

set outdir [getenv SIM_OUTDIR ""]
set checkpoint [getenv SIM_CKPNAME "para_part-"]

if {$ckp_output == "yes"} {
    puts "Creating dir $outdir"
    file mkdir $outdir
}

# Output all parameters
puts "========================="
puts "= Simulation Parameters ="
puts "========================="
puts "Size"
puts "├─> N_part:  $n_part [is_envset SIM_NPART]"
puts "└─> Box_l:   $box_l [is_envset SIM_BOXL]"
puts "    (rho:    [expr double($n_part) / $box_l / $box_l])"
puts "Time"
puts "├─> Dt:      $dt [is_envset SIM_DT]"
puts "├─> N_Cycle: $n_cycle [is_envset SIM_NCYCLE]"
puts "└─> N_Steps: $n_steps [is_envset SIM_NSTEPS]"
puts "Thermostat"
puts "├─> Gamma:   $gamma [is_envset SIM_GAMMA]"
puts "└─> Temp:    $temp [is_envset SIM_TEMP]"
puts "LJ Interaction"
puts "├─> Eps:     $eps [is_envset SIM_EPS]"
puts "├─> Sigma:   $sigma [is_envset SIM_SIGMA]"
puts "└─> Rcut:    $rcut [is_envset SIM_RCUT]"
puts "Bonding"
puts "└─> Bind rad:$bind_radius [is_envset SIM_BIND_RADIUS]"
puts "Misc"
puts "└─> Skin:    $skin [is_envset SIM_SKIN]"
puts "Output/Visualization"
puts "├─> VMD Visu: $vmd_output [is_envset SIM_VMDOUTPUT]"
puts "├─> Ckpoints: $ckp_output [is_envset SIM_CKPOUTPUT]"
puts "├─> Ckp-type: $ckp_type [is_envset SIM_CKPTYPE]"
puts "├─> Outdir:   $outdir [is_envset SIM_OUTDIR]"
puts "└─> Ckp-name: $checkpoint [is_envset SIM_CKPNAME]"
puts ""

# Simulation Environment
setmd box_l $box_l $box_l $box_l
setmd time_step $dt
setmd skin $skin
setmd periodic 1 1 1

# Interactions
inter 0 0 lennard-jones $eps $sigma $rcut auto
inter 2 angle 1.0 [PI]
#inter 0 harmonic 959.6 1.0
inter 1 harmonic 2.0 1.0

# set particles
puts "\[Setup\] Setting $n_part particles"
for {set i 0} {$i < $n_part} {incr i} {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz 0.0
    set vx 0.0
    set vy 0.0
    set vz 0.0
    part $i pos $posx $posy $posz type 0 v $vx $vy $vz fix 0 0 1
}


if { $vmd_output=="yes" } {
    puts "\[VMD Visu\] Run `vmd -e para_part.vmd_start.script`."
    prepare_vmd_connection "para_part" wait 30000
}

# limit forces to allow smooth adaption
puts "\[Setup\] Warming up"
set nwarmupiter 0
set capincr 10
set cap 10
set mindist 0

set warmup_us [lindex [time {while { $mindist < 0.9 } {
    puts "  \[Warmup\] Cap: $cap, Mindist: $mindist"
    inter forcecap $cap
    integrate 40
    if { $vmd_output=="yes" } { imd positions }
    set mindist [analyze mindist 0 0]
    incr nwarmupiter
    incr cap $capincr
}}] 0]
puts "  \[Warmup\] Final mindist: $mindist"
puts "\[Setup\] Warmup took [format_sec $warmup_us] ($nwarmupiter iterations)"
# turn off the ljforcecap
inter forcecap 0

# Do not turn on bonding until after warmup.
thermostat langevin $temp $gamma
on_collision bind_three_particles $bind_radius 1 2 0

if {$ckp_output == "yes"} {
    dump_config
}

puts "\[Simulation\] Start"
for {set iter 0} { $iter < $n_cycle } {incr iter} {
    set int_us [lindex [time {
        integrate $n_steps
    }] 0]
    puts "Cycle [expr $iter + 1]/$n_cycle done (Took: [format_sec $int_us])"
    if { $vmd_output == "yes" } { imd positions }
    if { $ckp_output == "yes" } { write_checkpoint $iter }
}
puts "\[Simulation\] End."

