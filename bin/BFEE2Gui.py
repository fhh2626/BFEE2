#!/usr/bin/env python3

## binding free energy estimator (BFEE) v2.0 alpha
## by Haohao Fu (fhh2626_at_gmail.com)
##
## this is a completely rewritten version of BFEE
##
##
## Changelog:
##    2018.2.6 v0.4 beta:
##        first public version
##    2018.11.8 v0.5 beta:
##        H-mass repartitioning supported and numstep reduced
##        added water box based on the protein-ligand distance vector
##        reduced frequency of writing colvars outputs
##        now default pressure is 1.01325 bar
##    2018.11.12 v0.51 beta:
##        more strict check for polar_theta and polar_phi to guarantee the rotational invariance
##        by default, WTM-eABF instead classical eABF is used to accelerate sampling
##    2018.11.13 v0.52 beta:
##        minor bug fixes
##    2018.11.20 v0.53 beta:
##        removed H-mass repartitioning, which is not stable in free-energy calculation
##        by default, colvarsRestartFrequency=outputFreq
##        reduced stepsPerCycle for stability
##    2020.6.12 v0.6 beta:
##        adapted the latest version of Colvars
##        uses stochastic velocity rescaling instead of Langevin barostat to improve efficiency
##    2020.6.15 v1.0:
##        alchemical free-energy calculation now supported
##    2020.6.16 v1.01:
##        minor bug fixes
##    2020.6.18 v1.02:
##        switched to bidirectional FEP by default
##    2020.7.16 v1.1:
##        now users can set reference and protein seperately
##        making tuning the pulling direction easily
##    2020.11.15 v1.11:
##        minor bug fixes
##

import sys
from PySide2.QtWidgets import QApplication
import BFEE2.gui as gui
        
if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = gui.mainUI()
    sys.exit(app.exec_())
