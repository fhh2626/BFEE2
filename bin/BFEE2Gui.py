#!/usr/bin/env python3
###################################################
## binding free energy estimator (BFEE) v2.x
## by Haohao Fu (fhh2626_at_gmail.com)
##
## this is a completely rewritten version of BFEE
##
###################################################

import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QAction
import BFEE2.gui as gui
        
if __name__ == '__main__':
    
    app = QApplication(sys.argv)
    ex = gui.mainUI()
    sys.exit(app.exec())
