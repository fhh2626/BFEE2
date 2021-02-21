##############################################################
# FEP SCRIPT
# Jerome Henin <henin@ibpc.fr>
#
# Changes:
# 2020-05-03: added prev_lambda option to runFEP, improved doc
# 2018-04-18: added interleaved double-wide sampling (IDWS)
# 2010-04-24: added runFEPmin
# 2009-11-17: changed for NAMD 2.7 keywords
# 2008-06-25: added TI routines
# 2007-11-01: fixed runFEP to handle backwards transformations
#             (i.e. dLambda < 0)
##############################################################

##############################################################
## Example NAMD input: calculation with smaller windows at
## the ends
#
# source fep.tcl
#
# alch                  on
# alchFile              system.fep
# alchCol               B
# alchOutFreq           10
# alchOutFile           system.fepout
# alchEquilSteps        500
#
# set nSteps      5000
#
## A) Simple schedule: 20 windows from 0 to 1, in a single run
#
# runFEP 0.0 1.0 0.05 $nSteps
#
## B) Same thing, in two NAMD separate runs with a restart
#
## First run
# runFEP 0.0 0.5 0.05 $nSteps
#
## Restart
# runFEP 0.5 1.0 0.05 $nSteps
#
## C) Lambda schedule with narrower windows at the end points
## Using two explicit lists of lambda points
#
# set init {0.0 0.05 0.1}
# set end {0.9 0.95 1.0}
#
# runFEPlist $init $nSteps
# runFEP 0.1 0.9 0.1 $nSteps
# runFEPlist $end $nSteps
#
## Alternately, in one step:
#
# runFEPlist [concat $init [FEPlist 0.1 0.9 0.1] $end] $nSteps
#
##############################################################

##############################################################
## Special usage for Interleaved Double-Wide sampling
## A) Simple schedule: 20 windows from 0 to 1, in a single run
#
# runFEP 0.0 1.0 0.05 $nSteps true
#
## B) Same thing, in two NAMD separate runs with a restart
#
## First run
# runFEP 0.0 0.5 0.05 $nSteps true
#
## Restart - need to tell the script the previous lambda point: 0.45
# runFEP 0.5 1.0 0.05 $nSteps true 0.45
#
## C) Example of a piecewise calculation with restarts
# and a nonlinear lambda schedule
#
## Run individual points 0, 0.05 then the series from 0.1 to 0.5
#
# runFEPlist [concat {0. 0.05} [FEPlist 0.1 0.5 0.1]] $numSteps true
#
## Continue series from 0.5 to 0.9, sampling backward dE from 0.4
#
# runFEPlist [FEPlist 0.5 0.9 0.1] $numSteps true 0.4
#
## Add two values 0.95 and 1, sampling backward dE from 0.9
## (automatically adds final backward window from 1. to 0.95)
#
# runFEPlist {0.95 1.} $numSteps true 0.9
#
##############################################################

##############################################################
# proc runFEPlist { lambdaList nSteps {IDWS} {prev_lambda} }
#
# Run n FEP windows joining (n + 1) lambda-points
# Provide prev_lambda value if continuing a sequential
# transformation and using IDWS
##############################################################

proc runFEPlist { lambdaList nSteps { IDWS false } { prev_lambda -1 } } {
    set epsilon 1e-15

    # Keep track of window number
    global win
    if {![info exists win]} {
      set win 1
    }

    set l1 [lindex $lambdaList 0]
    foreach l2 [lrange $lambdaList 1 end] {
      print [format "Running FEP window %3s: Lambda1 %-6s Lambda2 %-6s \[dLambda %-6s\]"\
        $win $l1 $l2 [expr $l2 - $l1]]
      firsttimestep 0
      alchLambda       $l1
      alchLambda2      $l2

      if { $IDWS && ($prev_lambda >= 0.) } {
        alchLambdaIDWS $prev_lambda
      }
      run $nSteps

      # Keep track of previous value to set is as target for backward calculation in IDWS
      set prev_lambda $l1
      set l1 $l2
      incr win
    }

    if { $IDWS && ($l1 > [expr {1. - $epsilon}] || $l1 < $epsilon)} {
      # If the list ends at 1 or zero, we add a final window, which is backward from the end point
      # to complete double-wide sampling
      # this will be look like "forward" sampling in the fepout file ("FepEnergy:" keyword)
      print [format "Running FEP window %3s: Lambda1 %-6s Lambda2 %-6s \[dLambda %-6s\]"\
        $win $l1 $l2 [expr $l2 - $l1]]
      firsttimestep 0
      alchLambda       $l1
      alchLambda2      $prev_lambda
      alchLambdaIDWS   -1
      run $nSteps
    }
}


##############################################################
# proc runFEP { start stop dLambda nSteps {IDWS} {prev_lambda} }
#
# run FEP windows of width dLambda between values start and stop
##############################################################

proc runFEP { start stop dLambda nSteps { IDWS false } { prev_lambda -1 } } {

    runFEPlist [FEPlist $start $stop $dLambda] $nSteps $IDWS $prev_lambda
}


##############################################################
# proc FEPlist { start stop dLambda nSteps }
#
# Create list of FEP windows
##############################################################

proc FEPlist { start stop dLambda } {
    set epsilon 1e-15

    if { ($stop < $start) && ($dLambda > 0) } {
      set dLambda [expr {-$dLambda}]
    }

    if { $start == $stop } {
      set ll [list $start $start]
    } else {
      set ll [list $start]
      set l2 [increment $start $dLambda]

      if { $dLambda > 0} {
        # A small workaround for numerical rounding errors
        while { [expr {$l2 <= ($stop + $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      } else {
        while { [expr {$l2 >= ($stop - $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      }
    }

    return $ll
}


##############################################################
##############################################################

proc runFEPmin { start stop dLambda nSteps nMinSteps temp} {
    set epsilon 1e-15

    if { ($stop < $start) && ($dLambda > 0) } {
      set dLambda [expr {-$dLambda}]
    }

    if { $start == $stop } {
      set ll [list $start $start]
    } else {
      set ll [list $start]
      set l2 [increment $start $dLambda]

      if { $dLambda > 0} {
        # A small workaround for numerical rounding errors
        while { [expr {$l2 <= ($stop + $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      } else {
        while { [expr {$l2 >= ($stop - $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      }
    }

    if { $nMinSteps > 0 } {
      alchLambda       $start
      alchLambda2      $start
      minimize $nMinSteps
      reinitvels $temp
    }

    runFEPlist $ll $nSteps
}

##############################################################
##############################################################

proc runTIlist { lambdaList nSteps } {
    # Keep track of window number
    global win
    if {![info exists win]} {
	    set win 1
    }

    foreach l $lambdaList {
	    print [format "Running TI window %3s: Lambda %-6s "	$win $l ]
	    firsttimestep 0
	    alchLambda       $l
	    run $nSteps
	    incr win
    }
}


##############################################################
##############################################################

proc runTI { start stop dLambda nSteps } {
    set epsilon 1e-15

    if { ($stop < $start) && ($dLambda > 0) } {
      set dLambda [expr {-$dLambda}]
    }

    if { $start == $stop } {
      set ll [list $start $start]
    } else {
      set ll [list $start]
      set l2 [increment $start $dLambda]

      if { $dLambda > 0} {
        # A small workaround for numerical rounding errors
        while { [expr {$l2 <= ($stop + $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      } else {
        while { [expr {$l2 >= ($stop - $epsilon) } ] } {
          lappend ll $l2
          set l2 [increment $l2 $dLambda]
        }
      }
    }

    runTIlist $ll $nSteps
}

##############################################################
# Increment lambda and try to correct truncation errors around
# 0 and 1
##############################################################

proc increment { lambda dLambda } {
    set epsilon 1e-15
    set new [expr { $lambda + $dLambda }]

    if { [expr $new > - $epsilon && $new < $epsilon] } {
      return 0.0
    }
    if { [expr ($new - 1) > - $epsilon && ($new - 1) < $epsilon] } {
      return 1.0
    }
    return $new
}
