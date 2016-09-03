#!/bin/bash

NODES="1"
PARTICLES="100 1000 10000 100000"
DENSITY="0.2 0.5 0.7"
PROCS="1 2 3 4 5 6 7 8 12 16"
ID=$(date +"%Y-%m-%d-%H-%M")
SIM_TEMPLATE="/tikhome/tiwa/sim/simulation_template.sh"
HEADER="/tikhome/tiwa/sim/HEADER"
FOOTER="/tikhome/tiwa/sim/FOOTER"

mkdir $ID
pushd $ID


for part in $PARTICLES 
do
  for den in $DENSITY
  do
    for process in $PROCS
    do
        mkdir "N$Nodes-n$process-p$part-d$den"
        pushd "N$Nodes-n$process-p$part-d$den"

        cp $SIM_TEMPLATE script.sh
        sed -i "s/#SBATCH -n 4/#SBATCH -n $process/" script.sh
        sed -i "s/#SBATCH -J espressomd-testsim/#SBATCH -J sim-part-$part-den-$den-procs-$process/" script.sh
        NAME="$(pwd)/lj_liquid.tcl"
        sed -i "s!TCLscriptNAME!$NAME!" script.sh

        cat $HEADER > lj_liquid.tcl
        echo 'set max_time              7200' >> lj_liquid.tcl
        echo 'set backup_timer          100' >> lj_liquid.tcl
        echo "set density               $den" >> lj_liquid.tcl
        echo "set n_part                $part" >> lj_liquid.tcl
        echo 'set int_n_times           10' >> lj_liquid.tcl
        echo 'set int_steps             1000' >> lj_liquid.tcl
        echo 'set sim_type              "Exclusive"' >> lj_liquid.tcl
        echo "set sim_nodes             $NODES" >> lj_liquid.tcl
        echo "set sim_procs             $process" >> lj_liquid.tcl
        echo 'set espresso_version      "ESPResSo-unknown"' >> lj_liquid.tcl
        echo 'set boost_version         "1.60.0"' >> lj_liquid.tcl
        echo 'set mpi_version           "1.10.1"' >> lj_liquid.tcl
        echo "set simulation_date       $ID" >> lj_liquid.tcl
        echo 'set modules_loaded        "look into script.sh"' >> lj_liquid.tcl
        echo 'set simulation_name       "look at dirname"' >> lj_liquid.tcl
        echo 'set custom_compileflags   "none"' >> lj_liquid.tcl
        echo 'set measure_sdterr        "yes"' >> lj_liquid.tcl
        echo 'set number_of_timers      11' >> lj_liquid.tcl
        echo 'set output_file           [open "simulation.output" "w"]' >> lj_liquid.tcl
        echo 'set backup_file           [open "simulation.backup" "w"]' >> lj_liquid.tcl
        cat $FOOTER >> lj_liquid.tcl
        sbatch script.sh
        popd
    done
  done
done
popd 

