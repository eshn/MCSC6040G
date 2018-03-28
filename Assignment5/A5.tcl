#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
exec $ESPRESSO_SOURCE/Espresso $0 $*

set N1 40
set N2 5
set rseed 1234

# Seed
t_random seed $rseed
set first [expr srand($rseed)]

set PI 3.14159265359

set tmax 1000

# Boundaries
set boxx 100
set boxy 10
set boxz 10

# File IO
# First 2 data files are for time analysis
set cardat [open [join [list car.dat] {}] "w"]
set truckdat [open [join [list truck.dat] {}] "w"]
set ftrj [open [join [list test.xyz] {}] "w"]

# Center of the domain
set cx [expr $boxx/2.0]
set cy [expr $boxy/2.0]

# Lennard Jones Parameters
set eps 1.0
set sig 1.0
set rc [expr pow(2.0,1.0/6.0)*$sig]
set shi 0.25

thermostat langevin 0.00 1.0

# Initialization of domain
setmd box_l $boxx $boxy $boxz
setmd periodic 1 1 1
setmd time_step 0.01
setmd skin 0.4

## Potentials
# Lennard Jones (inter <type> <type> <potential> <params>)
inter 1 1 lennard-jones $eps $sig $rc $shi 0.0
inter 1 0 lennard-jones $eps $sig $rc $shi 0.0
inter 2 0 lennard-jones $eps $sig $rc $shi 0.0

inter 1 2 lennard-jones $eps $sig $rc $shi 0.0

inter 2 2 lennard-jones $eps $sig $rc $shi 0.0

inter 0 fene 30 2
inter 3 angle 1 $PI

# Boundaries
for { set i 1000 } { $i < [expr $boxx + 1000.0] } { incr i} {
	set x [expr $i - 1000.0]
	set y [expr $boxy]
	part $i pos $x $y 0.0 type 0 fix 1 1 1 
	part [expr $i+1000] pos $x [expr $boxy] 0.0 type 0
}
# Obstacles
for {set j 1} { $j < 9 } { incr j} {
	for {set i 0 } { $i < 4} {incr i } {
		set x [expr $boxx / 2 + 4*$j + $i]
		if {$j%2 == 0} {
			set y [expr $i + 1]
		}
		if {$j%2 == 1} {
			set y [expr $boxy - $i -1]
		}
		part [expr 3000 + 4*($j - 1) + $i] pos $x $y 0.0 type 0 fix 1 1 1
	}
}

# Particles
# Regular cars
set yoff 2
for { set i 0 } { $i < $N1 } { incr i } {
	set x [expr $i%10 *1.5 + 1]
	if {$i%10 == 0} {
		set yoff [expr $yoff + 1.5]
	}
	set y $yoff
	set z 0.0
	part $i pos $x $y $z type 1 ext_force 0.5 0 0
}

#Longer cars (length 5)
for {set j 0} {$j < $N2} {incr j} {
	set x [expr ($j+1) * 5 + 25]
	set y [expr ($j+1) * $boxy / ($N2 + 1)]
	for { set i 0 } { $i < 5} {incr i} {
		set ind [expr 100 + $j * $N2 + $i]
		if {$i == 0} {
			part $ind pos $x $y 0.0 type 2 ext_force 0.5 0.0 0.0
		}
		if {$i != 0} {
			part $ind pos [expr $x - $i] $y 0.0 type 2 ext_force 0.0 0.0 0.0
			part $ind bond 0 [expr $ind - 1]
		}
	}
	for { set i 0 } { $i < 5} {incr i} {
		set ind [expr 100 + $j * $N2 + $i]
		if {$i != 0 && $i != 4} {
			part $ind bond 3 [expr $ind + 1] [expr $ind - 1]
		}
	}
}


# Write to file
for { set t 0 } { $t < $tmax } { incr t } {
	if { $t%2 == 0 } {
		puts $ftrj "[expr $N1 + 5*$N2 + 2*$boxx + 32]"
		puts $ftrj "title"

		puts $cardat "$N1"
		puts $truckdat "$N2"
		# Write cars
		for { set i 0 } { $i < $N1 } { incr i } {
			set x [lindex [part $i print folded_position] 0]
			set y [lindex [part $i print folded_position] 1]
			puts $ftrj "a$i $x $y 0.0"
			# Car info only
			puts $cardat "a$i $x"
		}
		# Write Longer Cars
		for { set j 0 } { $j < $N2 } { incr j} {
		    for {set i 0 } { $i < 5 } { incr i} {
		        set ind [expr 100 + $j * $N2 + $i]
		        set x [lindex [part $ind print folded_position] 0]
		        set y [lindex [part $ind print folded_position] 1]
			if {$i == 0} {
		        	puts $ftrj "a$ind $x $y 0.2"
				puts $truckdat "a$ind $x"
			}
			if {$i != 0 } {
				puts $ftrj "a$ind $x $y 0.3"
			}
		    }
		}
		# Write Walls
		for {set i 1000 } { $i < [expr $boxx + 1000] } { incr i } {
			set x [lindex [part $i print pos] 0]
			set y [lindex [part $i print pos] 1]
			puts $ftrj "a$i $x $y 0.1"
			puts $ftrj "a[expr $i+1000] $x 0.0 0.1"
		}
		# Write Obstacles
		for {set j 1} { $j < 9 } { incr j} {
	            for {set i 0 } { $i < 4} {incr i } {
        	        set ind [expr 3000 + 4*($j - 1) + $i]
	                set x [lindex [part $ind print pos] 0]
        	        set y [lindex [part $ind print pos] 1]
	                puts $ftrj "a$ind $x $y 0.1"
        	    }
	        }
	}
	integrate 100
}
