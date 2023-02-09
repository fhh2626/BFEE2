# the GUI of new BFEE

import os
import shutil
import sys
import webbrowser

import numpy as np
# use appdirs to manage persistent configuration
from appdirs import user_config_dir
from PySide2 import QtCore
from PySide2.QtGui import QFont, QIcon
from PySide2.QtWidgets import (QAction, QApplication, QCheckBox, QComboBox,
                               QFileDialog, QGridLayout, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QListWidget,
                               QMainWindow, QMessageBox, QPushButton,
                               QSplitter, QTabWidget, QToolBar, QVBoxLayout,
                               QWidget)

import BFEE2.inputGenerator as inputGenerator
import BFEE2.postTreatment as postTreatment
import BFEE2.templates_gromacs.BFEEGromacs as BFEEGromacs
from BFEE2.commonTools import commonSlots, fileParser, ploter

try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

import BFEE2.version
from BFEE2 import doc

__PROGRAM_NAME__ = f'BFEEstimator v{BFEE2.version.__VERSION__}'

class mainSettings(QWidget):
    """settings in the menubar
       set pathes of third-party softwares (VMD, gmx and tleap)
    """
    
    def __init__(self):
        super().__init__()
        self.config_dir = user_config_dir('BFEE2', 'chinfo')
        # test if config directory exists
        if not os.path.exists(self.config_dir):
            # create it if not exists
            os.makedirs(self.config_dir)
        self._initUI()
        self._initSingalsSlots()
        self.setWindowTitle('Settings')
        self._readConfig()
        #self.setGeometry(0,0,0,0)
        #self.show()

    def _initUI(self):
        """settings GUI
        """

        self.mainLayout = QVBoxLayout()

        self.thirdPartySoftware = QGroupBox('Third party software:')
        self.thirdPartySoftwareLayout = QVBoxLayout()

        # settings grid
        self.settingsGridLayout = QGridLayout()

        # vmd
        self.vmdLabel = QLabel('VMD:')
        self.vmdLineEdit = QLineEdit()
        self.vmdButton = QPushButton('Browse')
        self.settingsGridLayout.addWidget(self.vmdLabel, 0, 0)
        self.settingsGridLayout.addWidget(self.vmdLineEdit, 0, 1)
        self.settingsGridLayout.addWidget(self.vmdButton, 0, 2)

        # gmx
        #self.gromacsLayout = QHBoxLayout()
        #self.gromacsLabel = QLabel('Gromacs:')
        #self.gromacsLineEdit = QLineEdit()
        #self.gromacsButton = QPushButton('Browse')
        #self.settingsGridLayout.addWidget(self.gromacsLabel, 1, 0)
        #self.settingsGridLayout.addWidget(self.gromacsLineEdit,1 ,1)
        #self.settingsGridLayout.addWidget(self.gromacsButton, 1, 2)

        # tleap
        #self.tleapLayout = QHBoxLayout()
        #self.tleapLabel = QLabel('tleap:')
        #self.tleapLineEdit = QLineEdit()
        #self.tleapButton = QPushButton('Browse')
        #self.settingsGridLayout.addWidget(self.tleapLabel, 2, 0)
        #self.settingsGridLayout.addWidget(self.tleapLineEdit, 2, 1)
        #self.settingsGridLayout.addWidget(self.tleapButton, 2, 2)

        # OK and Cancel
        self.settingsButtonLayout = QHBoxLayout()
        self.settingsOKButton = QPushButton('OK')
        #self.settingsCancelButton = QPushButton('Cancel')
        self.settingsButtonLayout.addWidget(QSplitter())
        self.settingsButtonLayout.addWidget(self.settingsOKButton)
        #self.settingsButtonLayout.addWidget(self.settingsCancelButton)

        self.thirdPartySoftwareLayout.addLayout(self.settingsGridLayout)
        self.thirdPartySoftwareLayout.addLayout(self.settingsButtonLayout)

        self.thirdPartySoftware.setLayout(self.thirdPartySoftwareLayout)
        self.mainLayout.addWidget(self.thirdPartySoftware)
        self.setLayout(self.mainLayout)

    def _initSingalsSlots(self):
        """initialize singals and slots
        """
        self.vmdButton.clicked.connect(commonSlots.openFileDialog('exe', self.vmdLineEdit))
        #self.gromacsButton.clicked.connect(commonSlots.openFileDialog('exe', self.gromacsLineEdit))
        #self.tleapButton.clicked.connect(commonSlots.openFileDialog('exe', self.tleapLineEdit))
        self.settingsOKButton.clicked.connect(self._OKSlot())

    def _readConfig(self):
        """read the config saving paths for third-party softwares
        """
        if not os.path.exists(f'{self.config_dir}/3rdSoft.ini'):
            return

        with open(f'{self.config_dir}/3rdSoft.ini', 'r') as cFile:
            line = cFile.readline()
            self.vmdLineEdit.setText(line.strip())
            #line = cFile.readline()
            #self.gromacsLineEdit.setText(line.strip())
            #line = cFile.readline()
            #self.tleapLineEdit.setText(line.strip())

    def _writeConfig(self):
        """write the config saving paths for third-party softwares
        """

        with open(f'{self.config_dir}/3rdSoft.ini', 'w') as cFile:
            cFile.write(self.vmdLineEdit.text() + '\n')
            #cFile.write(self.gromacsLineEdit.text() + '\n')
            #cFile.write(self.tleapLineEdit.text() + '\n')

    def _OKSlot(self):
        """the slot corresponding the OK button
        """
        def f(): 
            self._writeConfig()
            self.close()
        return f

class geometricAdvancedSettings(QWidget):
    """advanced settings for the geometric route
       set pulling direction, non-standard solvent
       and number of stratification windows
    """

    def __init__(self):
        super().__init__()
        self._initUI()
        self._initSingalsSlots()
        self.setWindowTitle('Geometric advanced settings')
        self.setGeometry(0,0,0,0)
        #self.show()

    def _initUI(self):
        """initialize UI for the geometric advanced settings
        """

        self.mainLayout = QVBoxLayout()

        # user-defined pulling direction
        self.userDefinedDirection = QGroupBox('User-defined separation direction')
        self.userDefinedDirectionLayout = QHBoxLayout()

        self.userDefinedDirectionLabel = QLabel('Reference:')
        self.userDefinedDirectionLineEdit = QLineEdit()
        self.userDefinedDirectionLayout.addWidget(self.userDefinedDirectionLabel)
        self.userDefinedDirectionLayout.addWidget(self.userDefinedDirectionLineEdit)

        self.userDefinedDirection.setLayout(self.userDefinedDirectionLayout)

        # non-standard solvent
        self.nonStandardSolvent = QGroupBox('User-provided large box')
        self.nonStandardSolventLayout = QGridLayout()

        self.nonStandardSolventPsfLabel = QLabel('psf/parm file:')
        self.nonStandardSolventPsfLineEdit = QLineEdit()
        self.nonStandardSolventPsfButton = QPushButton('Browse')
        self.nonStandardSolventPdbLabel = QLabel('pdb/rst7 file:')
        self.nonStandardSolventPdbLineEdit = QLineEdit()
        self.nonStandardSolventPdbButton = QPushButton('Browse')
        self.nonStandardSolventLayout.addWidget(self.nonStandardSolventPsfLabel, 0, 0)
        self.nonStandardSolventLayout.addWidget(self.nonStandardSolventPsfLineEdit, 0, 1)
        self.nonStandardSolventLayout.addWidget(self.nonStandardSolventPsfButton, 0, 2)
        self.nonStandardSolventLayout.addWidget(self.nonStandardSolventPdbLabel, 1, 0)
        self.nonStandardSolventLayout.addWidget(self.nonStandardSolventPdbLineEdit, 1, 1)
        self.nonStandardSolventLayout.addWidget(self.nonStandardSolventPdbButton, 1, 2)

        self.nonStandardSolvent.setLayout(self.nonStandardSolventLayout)

        # stratification
        self.stratification = QGroupBox('Stratification (number of strata)')
        self.stratificationLayout = QGridLayout()

        self.stratificationRMSDBoundLabel = QLabel('RMSD(Bound):')
        self.stratificationRMSDBoundLineEdit = QLineEdit('1')
        self.stratificationTheta = QLabel('Theta:')
        self.stratificationThetaLineEdit = QLineEdit('1')
        self.stratificationPhiLabel = QLabel('Phi:')
        self.stratificationPhiLineEdit = QLineEdit('1')
        self.stratificationPsiLabel = QLabel('Psi:')
        self.stratificationPsiLineEdit = QLineEdit('1')
        self.stratificationthetaLabel = QLabel('theta:')
        self.stratificationthetaLineEdit = QLineEdit('1')
        self.stratificationphiLabel = QLabel('phi:')
        self.stratificationphiLineEdit = QLineEdit('1')
        self.stratificationRLabel = QLabel('r:')
        self.stratificationRLineEdit = QLineEdit('1')
        self.stratificationRMSDUnboundLabel = QLabel('RMSD(Unbound):')
        self.stratificationRMSDUnboundLineEdit = QLineEdit('1')

        self.stratificationLayout.addWidget(self.stratificationRMSDBoundLabel, 0, 0)
        self.stratificationLayout.addWidget(self.stratificationRMSDBoundLineEdit, 0, 1)
        self.stratificationLayout.addWidget(self.stratificationTheta, 0, 2)
        self.stratificationLayout.addWidget(self.stratificationThetaLineEdit, 0, 3)
        self.stratificationLayout.addWidget(self.stratificationPhiLabel, 0, 4)
        self.stratificationLayout.addWidget(self.stratificationPhiLineEdit, 0, 5)
        self.stratificationLayout.addWidget(self.stratificationPsiLabel, 0, 6)
        self.stratificationLayout.addWidget(self.stratificationPsiLineEdit, 0, 7)
        self.stratificationLayout.addWidget(self.stratificationthetaLabel, 1, 0)
        self.stratificationLayout.addWidget(self.stratificationthetaLineEdit, 1, 1)
        self.stratificationLayout.addWidget(self.stratificationphiLabel, 1, 2)
        self.stratificationLayout.addWidget(self.stratificationphiLineEdit, 1, 3)
        self.stratificationLayout.addWidget(self.stratificationRLabel, 1, 4)
        self.stratificationLayout.addWidget(self.stratificationRLineEdit, 1, 5)
        self.stratificationLayout.addWidget(self.stratificationRMSDUnboundLabel, 1, 6)
        self.stratificationLayout.addWidget(self.stratificationRMSDUnboundLineEdit, 1, 7)

        self.stratification.setLayout(self.stratificationLayout)
        
        # compatibility
        self.compatibility = QGroupBox('Compatibility')
        self.compatibilityLayout = QGridLayout()
        
        self.pinDownProCheckbox = QCheckBox('Pinning down the protein')
        self.pinDownProCheckbox.setChecked(True)
        
        self.useOldCvCheckbox = QCheckBox('Use quaternion-based CVs')
        self.useOldCvCheckbox.setChecked(False)
        
        self.reflectingBoundaryCheckbox = QCheckBox('Use reflecting boundary')
        self.reflectingBoundaryCheckbox.setChecked(True)
        
        self.compatibilityLayout.addWidget(self.pinDownProCheckbox, 0, 0)
        self.compatibilityLayout.addWidget(self.useOldCvCheckbox, 0, 1)
        self.compatibilityLayout.addWidget(self.reflectingBoundaryCheckbox, 1, 0)
        self.compatibility.setLayout(self.compatibilityLayout)
        
        # force field settings
        self.FFSettings = QGroupBox('Force field settings')
        self.FFSettingsLayout = QHBoxLayout()
        
        self.OPLSMixingRuleCheckbox = QCheckBox('OPLS mixing rules')
        self.OPLSMixingRuleCheckbox.setChecked(False)
        
        self.FFSettingsLayout.addWidget(self.OPLSMixingRuleCheckbox)
        self.FFSettings.setLayout(self.FFSettingsLayout)

        # strategy settings
        self.strategy = QGroupBox('Strategy settings')
        self.strategyLayout = QHBoxLayout()

        self.considerRMSDCVCheckbox = QCheckBox('Take into account RMSD CV')
        self.considerRMSDCVCheckbox.setChecked(True)

        self.useGaWTMCheckbox = QCheckBox('Use GaWTM-eABF')
        self.useGaWTMCheckbox.setChecked(False)

        self.strategyLayout.addWidget(self.considerRMSDCVCheckbox)
        self.strategyLayout.addWidget(self.useGaWTMCheckbox)
        self.strategy.setLayout(self.strategyLayout)
        
        # membrane protein
        self.modeling = QGroupBox('Modeling (avaiable for CHARMM FF)')
        self.modelingLayout = QVBoxLayout()

        self.memProCheckbox = QCheckBox('Membrane protein')
        self.memProCheckbox.setChecked(False)
        
        self.neutralizeLigOnlyLayout = QHBoxLayout()
        self.neutralizeLigOnlyLabel = QLabel('Auto-neutralize ligand-only system by:')
        self.neutralizeLigOnlyCombobox = QComboBox()
        self.neutralizeLigOnlyCombobox.addItem('NaCl')
        self.neutralizeLigOnlyCombobox.addItem('KCl')
        self.neutralizeLigOnlyCombobox.addItem('CaCl2')
        self.neutralizeLigOnlyCombobox.addItem('None')
        self.neutralizeLigOnlyLayout.addWidget(self.neutralizeLigOnlyLabel)
        self.neutralizeLigOnlyLayout.addWidget(self.neutralizeLigOnlyCombobox)

        self.modelingLayout.addWidget(self.memProCheckbox)
        self.modelingLayout.addLayout(self.neutralizeLigOnlyLayout)
        self.modeling.setLayout(self.modelingLayout)

        # parallel runs for error estimation
        self.parallelRuns = QGroupBox('Parallel runs')
        self.parallelRunsLayout = QHBoxLayout()

        self.parallelRunsLabel = QLabel('Number of parallel runs: ')
        self.parallelRunsLineEdit = QLineEdit('1')

        self.parallelRunsLayout.addWidget(self.parallelRunsLabel)
        self.parallelRunsLayout.addWidget(self.parallelRunsLineEdit)
        self.parallelRuns.setLayout(self.parallelRunsLayout)
        

        self.geometricAdvancedSettingsButtonLayout = QHBoxLayout()
        self.geometricAdvancedSettingsOKButton = QPushButton('OK')
        #self.geometricAdvancedSettingsCancelButton = QPushButton('Cancel')
        self.geometricAdvancedSettingsButtonLayout.addWidget(QSplitter())
        self.geometricAdvancedSettingsButtonLayout.addWidget(self.geometricAdvancedSettingsOKButton)
        #self.geometricAdvancedSettingsButtonLayout.addWidget(self.geometricAdvancedSettingsCancelButton)

        self.mainLayout.addWidget(self.userDefinedDirection)
        self.mainLayout.addWidget(self.nonStandardSolvent)
        self.mainLayout.addWidget(self.stratification)
        self.mainLayout.addWidget(self.compatibility)
        self.mainLayout.addWidget(self.FFSettings)
        self.mainLayout.addWidget(self.strategy)
        self.mainLayout.addWidget(self.modeling)
        self.mainLayout.addWidget(self.parallelRuns)
        self.mainLayout.addLayout(self.geometricAdvancedSettingsButtonLayout)
        self.setLayout(self.mainLayout)

    def _initSingalsSlots(self):
        """initialize (connect) signial and slots for geometric advanced settings
        """

        self.nonStandardSolventPsfButton.clicked.connect(
            commonSlots.openFileDialog(
                'psf/parm/top', 
                self.nonStandardSolventPsfLineEdit
            )
        )
        self.nonStandardSolventPdbButton.clicked.connect(
            commonSlots.openFileDialog(
                'pdb/gro', 
                self.nonStandardSolventPdbLineEdit
            )
        )
        self.geometricAdvancedSettingsOKButton.clicked.connect(self.close)


class alchemicalAdvancedSettings(QWidget):
    """advanced settings for the alchemical route
       set the number of stratification windows
    """

    def __init__(self):
        super().__init__()
        self._initUI()
        self._initSingalsSlots()
        self.setWindowTitle('Alchemical advanced settings')
        self.setGeometry(0,0,0,0)
        #self.show()

    def _initUI(self):
        """initialize UI for the alchemical advanced settings
        """

        self.mainLayout = QVBoxLayout()

        # stratification windows
        self.stratification = QGroupBox('Stratification (number of strata)')
        self.stratificationLayout = QGridLayout()

        self.boundLigandLabel = QLabel('Ligand/Bound state:')
        self.boundLigandLineEdit = QLineEdit('50')
        self.unboundLigandLabel = QLabel('Ligand/Unbound state:')
        self.unboundLigandLineEdit = QLineEdit('20')
        self.boundRestraintsLabel = QLabel('Restraints/Bound state:')
        self.boundRestraintsLineEdit = QLineEdit('50')
        self.unboundRestraintsLabel = QLabel('Restraints/Unbound state:')
        self.unboundRestraintsLineEdit = QLineEdit('20')

        self.stratificationLayout.addWidget(self.boundLigandLabel, 0, 0)
        self.stratificationLayout.addWidget(self.boundLigandLineEdit, 0, 1)
        self.stratificationLayout.addWidget(self.unboundLigandLabel, 0, 2)
        self.stratificationLayout.addWidget(self.unboundLigandLineEdit, 0, 3)
        self.stratificationLayout.addWidget(self.boundRestraintsLabel, 1, 0)
        self.stratificationLayout.addWidget(self.boundRestraintsLineEdit, 1, 1)
        self.stratificationLayout.addWidget(self.unboundRestraintsLabel, 1, 2)
        self.stratificationLayout.addWidget(self.unboundRestraintsLineEdit, 1, 3)

        self.stratification.setLayout(self.stratificationLayout)

        # double-wide simulation
        self.doubleWide = QGroupBox('Double-wide sampling simulation')
        self.doubleWideLayout = QGridLayout()

        self.doubleWideCheckbox = QCheckBox('Generate input files for double-wide sampling')
        self.doubleWideCheckbox.setChecked(False)
        self.doubleWideLayout.addWidget(self.doubleWideCheckbox)
        self.doubleWide.setLayout(self.doubleWideLayout)

        # minimize before sampling in each window
        self.minBeforeSample = QGroupBox('Minimization before sampling')
        self.minBeforeSampleLayout = QVBoxLayout()

        self.minBeforeSampleCheckbox = QCheckBox('Minimize before sampling in each window')
        self.minBeforeSampleCheckbox.setChecked(False)
        self.minBeforeSampleLayout.addWidget(self.minBeforeSampleCheckbox)
        self.minBeforeSample.setLayout(self.minBeforeSampleLayout)
        
        # compatibility
        self.compatibility = QGroupBox('Compatibility')
        self.compatibilityLayout = QHBoxLayout()
        
        self.pinDownProCheckbox = QCheckBox('Pinning down the protein')
        self.pinDownProCheckbox.setChecked(True)
        
        self.useOldCvCheckbox = QCheckBox('Use quaternion-based CVs')
        self.useOldCvCheckbox.setChecked(False)
        
        self.compatibilityLayout.addWidget(self.pinDownProCheckbox)
        self.compatibilityLayout.addWidget(self.useOldCvCheckbox)
        self.compatibility.setLayout(self.compatibilityLayout)
        
        # force field settings
        self.FFSettings = QGroupBox('Force field settings')
        self.FFSettingsLayout = QHBoxLayout()
        
        self.OPLSMixingRuleCheckbox = QCheckBox('OPLS mixing rules')
        self.OPLSMixingRuleCheckbox.setChecked(False)
        
        self.FFSettingsLayout.addWidget(self.OPLSMixingRuleCheckbox)
        self.FFSettings.setLayout(self.FFSettingsLayout)

        # membrane protein
        self.modeling = QGroupBox('Modeling (avaiable for CHARMM FF)')
        self.modelingLayout = QVBoxLayout()

        self.memProCheckbox = QCheckBox('Membrane protein')
        self.memProCheckbox.setChecked(False)
        
        self.neutralizeLigOnlyLayout = QHBoxLayout()
        self.neutralizeLigOnlyLabel = QLabel('Auto-neutralize ligand-only system by:')
        self.neutralizeLigOnlyCombobox = QComboBox()
        self.neutralizeLigOnlyCombobox.addItem('NaCl')
        self.neutralizeLigOnlyCombobox.addItem('KCl')
        self.neutralizeLigOnlyCombobox.addItem('CaCl2')
        self.neutralizeLigOnlyCombobox.addItem('None')
        self.neutralizeLigOnlyLayout.addWidget(self.neutralizeLigOnlyLabel)
        self.neutralizeLigOnlyLayout.addWidget(self.neutralizeLigOnlyCombobox)

        self.modelingLayout.addWidget(self.memProCheckbox)
        self.modelingLayout.addLayout(self.neutralizeLigOnlyLayout)
        self.modeling.setLayout(self.modelingLayout)

        self.alchemicalAdvancedSettingsButtonLayout = QHBoxLayout()
        self.alchemicalAdvancedSettingsOKButton = QPushButton('OK')
        #self.alchemicalAdvancedSettingsCancelButton = QPushButton('Cancel')
        self.alchemicalAdvancedSettingsButtonLayout.addWidget(QSplitter())
        self.alchemicalAdvancedSettingsButtonLayout.addWidget(self.alchemicalAdvancedSettingsOKButton)
        #self.alchemicalAdvancedSettingsButtonLayout.addWidget(self.alchemicalAdvancedSettingsCancelButton)
        
        self.mainLayout.addWidget(self.stratification)
        self.mainLayout.addWidget(self.doubleWide)
        self.mainLayout.addWidget(self.compatibility)
        self.mainLayout.addWidget(self.FFSettings)
        self.mainLayout.addWidget(self.minBeforeSample)
        self.mainLayout.addWidget(self.modeling)
        self.mainLayout.addLayout(self.alchemicalAdvancedSettingsButtonLayout)
        self.setLayout(self.mainLayout)

    def _initSingalsSlots(self):
        """initialize (connect) signals and slots for the alchemical advanced settings
        """

        self.alchemicalAdvancedSettingsOKButton.clicked.connect(self.close)

class mainUI(QMainWindow):
    """the main window UI 
       include the preTreatment, postTreatment and QuickPlot tab
       the preTreatment tab include NAMD and Gromacs tab
       the postTreatment tab include geometric and alchemical tab
    """
    
    def __init__(self):
        super().__init__()
        
        self._initActions()

        self._initNAMDTab()
        self._initGromacsTab()
        self._initPreTreatmentTab()
        self._initQuickPlotTab()

        self._initGeometricTab()
        self._initAlchemicalTab()
        self._initPostTreatmentTab()

        self._initSingalsSlots()

        self._initMainUI()

        # other dialogs
        self.mainSettings = mainSettings()
        self.alchemicalAdvancedSettings = alchemicalAdvancedSettings()
        self.geometricAdvancedSettings = geometricAdvancedSettings()

        self.setGeometry(0,0,0,0)
        self.setWindowTitle(__PROGRAM_NAME__)    
        self.setWindowIcon(QIcon("BFEE2/icon/icon.png"))
        self.show()

    def _initActions(self):
        ''' initialize actions for menubar '''

        # settings
        self.settingsAction = QAction('&Settings', self)
        self.settingsAction.setStatusTip('Set pathes for third-party softwares')
        self.settingsAction.triggered.connect(self._mainSettings())

        # exit
        self.exitAction = QAction('&Exit', self)        
        self.exitAction.setStatusTip('Exit application')
        self.exitAction.triggered.connect(QApplication.quit)

        # help
        self.helpAction = QAction('&Help', self)
        self.helpAction.setStatusTip('Open user manual')
        self.helpAction.triggered.connect(self._openDocFile)

        # python API
        self.pythonAPIAction = QAction('&Python API', self)
        self.pythonAPIAction.setStatusTip('Open Python API Documentation')
        self.pythonAPIAction.triggered.connect(self._openPythonAPIFile)

        # about
        self.aboutAction = QAction('&About', self)
        self.aboutAction.setStatusTip('About BFEEstimator')
        self.aboutAction.triggered.connect(self._showAboutBox())

        
    def _initMainUI(self):
        """initialize main window
        """
        
        # status bar
        #self.statusBar()

        # menu bar
        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        self.fileMenu = menubar.addMenu('&File')
        self.fileMenu.addAction(self.settingsAction)
        self.fileMenu.addSeparator()
        self.fileMenu.addAction(self.exitAction)

        self.helpMenu = menubar.addMenu('&Help')
        self.helpMenu.addAction(self.helpAction)
        self.helpMenu.addAction(self.pythonAPIAction)
        self.helpMenu.addSeparator()
        self.helpMenu.addAction(self.aboutAction)

        # main layout
        self.mainLayout = QVBoxLayout()

        # title
        self.title = QLabel('Binding Free Energy Estimator')
        titleFont = QFont()
        titleFont.setBold(True)
        self.title.setFont(titleFont)
        self.titleBox = QGroupBox()
        self.titleBoxLayout = QVBoxLayout()
        self.titleBoxLayout.addWidget(self.title, alignment=QtCore.Qt.AlignCenter)
        self.titleBox.setLayout(self.titleBoxLayout)

        # tabs
        self.mainTabs = QTabWidget()
        #self.preTreatmentTab = QWidget()
        #self.postTreatmentTab = QWidget()
        #self.quickPlot = QWidget()

        self.mainTabs.addTab(self.preTreatmentTab, 'Pre-treatment')
        self.mainTabs.addTab(self.postTreatmentTab, 'Post-treatment')
        self.mainTabs.addTab(self.quickPlot, 'Quick-Plot')

        # main layout
        #self.mainLayout.addWidget(self.titleBox)
        self.mainLayout.addWidget(self.mainTabs)
        self.mainWidgit = QWidget()
        self.mainWidgit.setLayout(self.mainLayout)
        self.setCentralWidget(self.mainWidgit)

    def _initPreTreatmentTab(self):
        """initialize pre-treatment tab
        """

        self.preTreatmentTab = QWidget()
        
        # pre-treatment tabs
        # NAMD and gromacs
        self.preTreatmentMainTabs = QTabWidget()

        self.preTreatmentMainTabs.addTab(self.NAMDTab, 'NAMD/Gromacs(CHARMM/Amber files)')
        self.preTreatmentMainTabs.addTab(self.GromacsTab, 'Gromacs(Gromacs files)')

        self.preTreatmentMainLayout = QVBoxLayout()
        self.preTreatmentMainLayout.addWidget(self.preTreatmentMainTabs)

        # other parameters
        self.otherParameters = QGroupBox('Other parameters')
        self.otherParametersLayout = QVBoxLayout()

        # temperature, selection protein and ligand layout
        self.otherParametersChildLayout = QGridLayout()

        # temperature
        self.temperatureLabel = QLabel('Temperature:      ')
        self.temperatureLineEdit = QLineEdit('300')
        self.otherParametersChildLayout.addWidget(self.temperatureLabel, 0, 0)
        self.otherParametersChildLayout.addWidget(self.temperatureLineEdit, 0, 1)

        # select protein
        self.selectProteinLabel = QLabel('Select protein:   ')
        self.selectProteineLineEdit = QLineEdit('segid SH3D')
        self.otherParametersChildLayout.addWidget(self.selectProteinLabel, 1, 0)
        self.otherParametersChildLayout.addWidget(self.selectProteineLineEdit, 1, 1)

        # select ligand
        self.selectLigandLabel = QLabel('Select ligand:    ')
        self.selectLigandLineEdit = QLineEdit('segid PPRO')
        self.otherParametersChildLayout.addWidget(self.selectLigandLabel, 2, 0)
        self.otherParametersChildLayout.addWidget(self.selectLigandLineEdit, 2, 1)

        # select strategy
        self.selectStrategyLayout = QHBoxLayout()
        
        self.selectStrategyLabel = QLabel('Select MD engine and strategy: ')
        self.selectStrategyCombobox = QComboBox()
        self.selectStrategyCombobox.addItem('Geometric')
        self.selectStrategyCombobox.addItem('Alchemical')
        self.selectStrategyAdvancedButton = QPushButton('Advanced settings')
        
        self.selectMDEngineCombobox = QComboBox()
        self.selectMDEngineCombobox.addItem('NAMD')
        self.selectMDEngineCombobox.addItem('Gromacs')
        
        self.selectStrategyChildLayout = QHBoxLayout()
        self.selectStrategyChildLayout.addWidget(self.selectMDEngineCombobox)
        self.selectStrategyChildLayout.addWidget(self.selectStrategyCombobox)
        self.selectStrategyChildLayout.addWidget(self.selectStrategyAdvancedButton)

        self.selectStrategyLayout.addWidget(self.selectStrategyLabel)
        self.selectStrategyLayout.addLayout(self.selectStrategyChildLayout)
        

        # generate input button
        self.generateInputButton = QPushButton('Generate Inputs')
 
        self.otherParametersLayout.addLayout(self.otherParametersChildLayout)
        self.otherParametersLayout.addLayout(self.selectStrategyLayout)
        self.otherParameters.setLayout(self.otherParametersLayout)

        self.preTreatmentMainLayout.addWidget(self.otherParameters)
        self.preTreatmentMainLayout.addWidget(self.generateInputButton)

        self.preTreatmentTab.setLayout(self.preTreatmentMainLayout)

    def _initNAMDTab(self):
        """initialize NAMD Tab in pre-treatment Tab
        """

        self.NAMDTab = QWidget()
        self.NAMDTabMainLayout = QVBoxLayout()

        # inputs for the complex
        self.inputsForComplex = QGroupBox('Inputs for complex')
        self.inputsForComplexLayout = QGridLayout()

        # psf/parm
        self.psfLabel = QLabel('psf/parm file:')
        self.psfLineEdit = QLineEdit()
        self.psfButton = QPushButton('Browse')
        self.inputsForComplexLayout.addWidget(self.psfLabel, 0, 0)
        self.inputsForComplexLayout.addWidget(self.psfLineEdit, 0, 1)
        self.inputsForComplexLayout.addWidget(self.psfButton, 0, 2)

        # coor
        self.coorLabel = QLabel('pdb/rst file:')
        self.coorLineEdit = QLineEdit()
        self.coorButton = QPushButton('Browse')
        self.inputsForComplexLayout.addWidget(self.coorLabel, 1, 0)
        self.inputsForComplexLayout.addWidget(self.coorLineEdit, 1, 1)
        self.inputsForComplexLayout.addWidget(self.coorButton, 1, 2)

        # force fields
        self.forceFields = QGroupBox('Force fields')
        self.forceFieldsLayout = QVBoxLayout()

        # force field type
        self.forceFieldTypeLayout = QHBoxLayout()
        self.forceFieldTypeLabel = QLabel('Force field type:')
        self.forceFieldCombobox = QComboBox()
        self.forceFieldCombobox.addItem('CHARMM')
        self.forceFieldCombobox.addItem('Amber')

        self.forceFieldTypeLayout.addWidget(self.forceFieldTypeLabel)
        self.forceFieldTypeLayout.addWidget(self.forceFieldCombobox)

        # CHARMM force field files
        self.forceFieldFilesLayout = QVBoxLayout()
        self.forceFieldFilesLabel = QLabel('Force field files:')
        self.forceFieldFilesBox = QListWidget()
        self.forceFieldFilesChildLayout = QHBoxLayout()
        self.forceFieldAddButton = QPushButton('Add')
        self.forceFieldClearButton = QPushButton('Clear')
        self.forceFieldFilesChildLayout.addWidget(self.forceFieldAddButton)
        self.forceFieldFilesChildLayout.addWidget(self.forceFieldClearButton)
        self.forceFieldFilesLayout.addWidget(self.forceFieldFilesLabel)
        self.forceFieldFilesLayout.addWidget(self.forceFieldFilesBox)
        self.forceFieldFilesLayout.addLayout(self.forceFieldFilesChildLayout)

        self.forceFieldsLayout.addLayout(self.forceFieldTypeLayout)
        self.forceFieldsLayout.addLayout(self.forceFieldFilesLayout)

        self.inputsForComplex.setLayout(self.inputsForComplexLayout)
        self.forceFields.setLayout(self.forceFieldsLayout)
        self.NAMDTabMainLayout.addWidget(self.inputsForComplex)
        self.NAMDTabMainLayout.addWidget(self.forceFields)
        self.NAMDTab.setLayout(self.NAMDTabMainLayout)

    def _initGromacsTab(self):
        """initialize GMX Tab in pre-treatment Tab
        """

        self.GromacsTab = QWidget()
        self.GromacsTabMainLayout = QVBoxLayout()

        # inputs for the complex
        self.GromacsInputsForComplex = QGroupBox('Inputs for complex')
        self.GromacsInputsForComplexLayout = QGridLayout()

        # top
        self.topLabel = QLabel('Topology file: ')
        self.topLineEdit = QLineEdit()
        self.topButton = QPushButton('Browse')
        self.GromacsInputsForComplexLayout.addWidget(self.topLabel, 0, 0)
        self.GromacsInputsForComplexLayout.addWidget(self.topLineEdit, 0, 1)
        self.GromacsInputsForComplexLayout.addWidget(self.topButton, 0, 2)

        # gro
        self.gromacsPdbLabel = QLabel('Structure file: ')
        self.gromacsPdbLineEdit = QLineEdit()
        self.gromacsPdbButton = QPushButton('Browse')
        self.GromacsInputsForComplexLayout.addWidget(self.gromacsPdbLabel, 1, 0)
        self.GromacsInputsForComplexLayout.addWidget(self.gromacsPdbLineEdit, 1, 1)
        self.GromacsInputsForComplexLayout.addWidget(self.gromacsPdbButton, 1, 2)

        # structure file format
        self.gromacsStructureFileFormatLayout = QHBoxLayout()
        self.gromacsStructureFileFormatLabel = QLabel('Structure file format:')
        self.gromacsStructureFileFormatCombobox = QComboBox()
        self.gromacsStructureFileFormatCombobox.addItem('pdb')
        self.gromacsStructureFileFormatCombobox.addItem('xpdb')
        self.gromacsStructureFileFormatCombobox.setToolTip('Select "<b>xpdb</b>" if your PDB file has more than 9,999 number of residues')
        self.gromacsStructureFileFormatLayout.addWidget(self.gromacsStructureFileFormatLabel)
        self.gromacsStructureFileFormatLayout.addWidget(self.gromacsStructureFileFormatCombobox)
        self.GromacsInputsForComplexLayout.addLayout(self.gromacsStructureFileFormatLayout, 2, 0, 1, 3)

        self.GromacsInputsForComplexLayout.addWidget(QSplitter(), 3, 1)
        self.GromacsInputsForComplex.setLayout(self.GromacsInputsForComplexLayout)

        # inputs for the ligand-only system
        self.GromacsInputsLigandOnly = QGroupBox('Inputs for ligand-only system')
        self.GromacsInputsLigandOnlyLayout = QGridLayout()

        # top
        self.gromacsLigandOnlyTopLabel = QLabel('Topology file: ')
        self.gromacsLigandOnlyTopLineEdit = QLineEdit()
        self.gromacsLigandOnlyTopButton = QPushButton('Browse')
        self.GromacsInputsLigandOnlyLayout.addWidget(self.gromacsLigandOnlyTopLabel, 0, 0)
        self.GromacsInputsLigandOnlyLayout.addWidget(self.gromacsLigandOnlyTopLineEdit, 0, 1)
        self.GromacsInputsLigandOnlyLayout.addWidget(self.gromacsLigandOnlyTopButton, 0, 2)

        # gro
        self.gromacsLigandOnlyPdbLabel = QLabel('Structure file: ')
        self.gromacsLigandOnlyPdbLineEdit = QLineEdit()
        self.gromacsLigandOnlyPdbButton = QPushButton('Browse')
        self.GromacsInputsLigandOnlyLayout.addWidget(self.gromacsLigandOnlyPdbLabel, 1, 0)
        self.GromacsInputsLigandOnlyLayout.addWidget(self.gromacsLigandOnlyPdbLineEdit, 1, 1)
        self.GromacsInputsLigandOnlyLayout.addWidget(self.gromacsLigandOnlyPdbButton, 1, 2)

        # structure file format
        self.gromacsLigandOnlyStructureFileFormatLayout = QHBoxLayout()
        self.gromacsLigandOnlyStructureFileFormatLabel = QLabel('Structure file format:')
        self.gromacsLigandOnlyStructureFileFormatCombobox = QComboBox()
        self.gromacsLigandOnlyStructureFileFormatCombobox.addItem('pdb')
        self.gromacsLigandOnlyStructureFileFormatCombobox.addItem('xpdb')
        self.gromacsLigandOnlyStructureFileFormatCombobox.setToolTip('Select "<b>xpdb</b>" if your PDB file has more than 9,999 number of residues')
        self.gromacsLigandOnlyStructureFileFormatLayout.addWidget(self.gromacsLigandOnlyStructureFileFormatLabel)
        self.gromacsLigandOnlyStructureFileFormatLayout.addWidget(self.gromacsLigandOnlyStructureFileFormatCombobox)
        self.GromacsInputsLigandOnlyLayout.addLayout(self.gromacsLigandOnlyStructureFileFormatLayout, 2, 0, 1, 3)

        self.GromacsInputsLigandOnlyLayout.addWidget(QSplitter(), 3, 1)
        self.GromacsInputsLigandOnly.setLayout(self.GromacsInputsLigandOnlyLayout)

        self.GromacsTabMainLayout.addWidget(self.GromacsInputsForComplex)
        self.GromacsTabMainLayout.addWidget(self.GromacsInputsLigandOnly)

        self.GromacsTab.setLayout(self.GromacsTabMainLayout)

    def _initPostTreatmentTab(self):
        """initialize pre-treatment tab
        """

        self.postTreatmentTab = QWidget()

        # post-treatment tabs
        # Geometric and alchemical
        self.postTreatmentMainTabs = QTabWidget()

        self.postTreatmentMainTabs.addTab(self.geometricTab, 'Geometric')
        self.postTreatmentMainTabs.addTab(self.alchemicalTab, 'Alchemical')

        self.postTreatmentMainLayout = QVBoxLayout()
        self.postTreatmentMainLayout.addWidget(self.postTreatmentMainTabs)

        self.calculateButton = QPushButton('Calculate binding free energy')
        self.postTreatmentMainLayout.addWidget(self.calculateButton)

        self.postTreatmentTab.setLayout(self.postTreatmentMainLayout)

    def _initGeometricTab(self):
        """initialize geometric tab of post-treatment
        """

        self.geometricTab = QWidget()
        self.geometricTabLayout = QVBoxLayout()

        # pmf inputs
        self.pmfInputs = QGroupBox('PMF inputs (.czar.pmf/.UI.pmf):')
        self.pmfInputsLayout = QVBoxLayout()

        # bound stats
        self.boundStateLabel = QLabel('Bound state:')
        self.boundStateLayout = QGridLayout()

        # RMSD
        self.rmsdBoundLabel = QLabel('RMSD: ')
        self.rmsdBoundLineEdit = QLineEdit()
        self.rmsdBoundButton = QPushButton('Browse')
        self.boundStateLayout.addWidget(self.rmsdBoundLabel, 0, 0)
        self.boundStateLayout.addWidget(self.rmsdBoundLineEdit, 0, 1)
        self.boundStateLayout.addWidget(self.rmsdBoundButton, 0, 2)

        # Theta
        self.ThetaLabel = QLabel('Theta:')
        self.ThetaLineEdit = QLineEdit()
        self.ThetaButton = QPushButton('Browse')
        self.boundStateLayout.addWidget(self.ThetaLabel, 1, 0)
        self.boundStateLayout.addWidget(self.ThetaLineEdit, 1, 1)
        self.boundStateLayout.addWidget(self.ThetaButton, 1, 2)

        # Phi
        self.PhiLabel = QLabel('Phi:  ')
        self.PhiLineEdit = QLineEdit()
        self.PhiButton = QPushButton('Browse')
        self.boundStateLayout.addWidget(self.PhiLabel, 2, 0)
        self.boundStateLayout.addWidget(self.PhiLineEdit, 2, 1)
        self.boundStateLayout.addWidget(self.PhiButton, 2, 2)

        # Psi
        self.PsiLabel = QLabel('Psi:  ')
        self.PsiLineEdit = QLineEdit()
        self.PsiButton = QPushButton('Browse')
        self.boundStateLayout.addWidget(self.PsiLabel, 3, 0)
        self.boundStateLayout.addWidget(self.PsiLineEdit, 3, 1)
        self.boundStateLayout.addWidget(self.PsiButton, 3, 2)

        # theta
        self.thetaLayout = QHBoxLayout()
        self.thetaLabel = QLabel('theta:')
        self.thetaLineEdit = QLineEdit()
        self.thetaButton = QPushButton('Browse')
        self.boundStateLayout.addWidget(self.thetaLabel, 4, 0)
        self.boundStateLayout.addWidget(self.thetaLineEdit, 4, 1)
        self.boundStateLayout.addWidget(self.thetaButton, 4, 2)

        # phi
        self.phiLayout = QHBoxLayout()
        self.phiLabel = QLabel('phi:  ')
        self.phiLineEdit = QLineEdit()
        self.phiButton = QPushButton('Browse')
        self.boundStateLayout.addWidget(self.phiLabel, 5, 0)
        self.boundStateLayout.addWidget(self.phiLineEdit, 5, 1)
        self.boundStateLayout.addWidget(self.phiButton, 5, 2)

        # r
        self.rLabel = QLabel('r:    ')
        self.rLineEdit = QLineEdit()
        self.rButton = QPushButton('Browse')
        self.boundStateLayout.addWidget(self.rLabel, 6, 0)
        self.boundStateLayout.addWidget(self.rLineEdit, 6, 1)
        self.boundStateLayout.addWidget(self.rButton, 6, 2)

        self.unboundStateLabel = QLabel('Unbound state:')

        # RMSD unbound
        self.rmsdUnboundLayout = QHBoxLayout()
        self.rmsdUnboundLabel = QLabel('RMSD: ')
        self.rmsdUnboundLineEdit = QLineEdit()
        self.rmsdUnboundButton = QPushButton('Browse')
        self.rmsdUnboundLayout.addWidget(self.rmsdUnboundLabel)
        self.rmsdUnboundLayout.addWidget(self.rmsdUnboundLineEdit)
        self.rmsdUnboundLayout.addWidget(self.rmsdUnboundButton)

        self.pmfInputsLayout.addWidget(self.boundStateLabel)
        self.pmfInputsLayout.addLayout(self.boundStateLayout)
        self.pmfInputsLayout.addWidget(self.unboundStateLabel)
        self.pmfInputsLayout.addLayout(self.rmsdUnboundLayout)

        self.pmfInputs.setLayout(self.pmfInputsLayout)

        # force constants
        self.forceConstants = QGroupBox('Force constants (in Colvars unit):')
        self.forceConstantsLayout = QGridLayout()

        self.fcBoundStateLabel = QLabel('Bound state:')

        # all widgets
        self.fcRMSDLabel = QLabel('RMSD:')
        self.fcRMSDLineEdit = QLineEdit('10')
        self.fcThetaLabel = QLabel('Theta:')
        self.fcThetaLineEdit = QLineEdit('0.1')
        self.fcPhiLabel = QLabel('Phi:')
        self.fcPhiLineEdit = QLineEdit('0.1')
        self.fcPsiLabel = QLabel('Psi:')
        self.fcPsiLineEdit = QLineEdit('0.1')
        self.fcthetaLabel = QLabel('theta:')
        self.fcthetaLineEdit = QLineEdit('0.1')
        self.fcphiLabel = QLabel('phi:')
        self.fcphiLineEdit = QLineEdit('0.1')

        self.forceConstantsLayout.addWidget(self.fcRMSDLabel, 0, 0)
        self.forceConstantsLayout.addWidget(self.fcRMSDLineEdit, 0, 1)
        self.forceConstantsLayout.addWidget(self.fcThetaLabel, 0, 2)
        self.forceConstantsLayout.addWidget(self.fcThetaLineEdit, 0, 3)
        self.forceConstantsLayout.addWidget(self.fcPhiLabel, 0, 4)
        self.forceConstantsLayout.addWidget(self.fcPhiLineEdit, 0, 5)
        self.forceConstantsLayout.addWidget(self.fcPsiLabel, 1, 0)
        self.forceConstantsLayout.addWidget(self.fcPsiLineEdit, 1, 1)
        self.forceConstantsLayout.addWidget(self.fcthetaLabel, 1, 2)
        self.forceConstantsLayout.addWidget(self.fcthetaLineEdit, 1, 3)
        self.forceConstantsLayout.addWidget(self.fcphiLabel, 1, 4)
        self.forceConstantsLayout.addWidget(self.fcphiLineEdit, 1, 5)

        self.forceConstants.setLayout(self.forceConstantsLayout)

        # other parameters
        self.postOtherParams = QGroupBox('Other parameters:')
        self.postOtherParamsLayout = QHBoxLayout()

        self.postTemperatureLabel = QLabel('temperature:')
        self.postTemperatureLineEdit = QLineEdit('300')
        self.postRstarLabel = QLabel(' r*:')
        self.postRstarLineEdit = QLineEdit('30')
        self.postPMFTypeLabel = QLabel(' PMF type:')
        self.postPMFTypeBox = QComboBox()
        self.postPMFTypeBox.addItem('NAMD')
        self.postPMFTypeBox.addItem('Gromacs')

        self.postOtherParamsLayout.addWidget(self.postTemperatureLabel)
        self.postOtherParamsLayout.addWidget(self.postTemperatureLineEdit)
        self.postOtherParamsLayout.addWidget(self.postRstarLabel)
        self.postOtherParamsLayout.addWidget(self.postRstarLineEdit)
        self.postOtherParamsLayout.addWidget(self.postPMFTypeLabel)
        self.postOtherParamsLayout.addWidget(self.postPMFTypeBox)

        self.postOtherParams.setLayout(self.postOtherParamsLayout)

        self.fcBoundStateLabel = QLabel('Bound state:')

        self.geometricTabLayout.addWidget(self.pmfInputs)
        self.geometricTabLayout.addWidget(self.forceConstants)
        self.geometricTabLayout.addWidget(self.postOtherParams)
        self.geometricTab.setLayout(self.geometricTabLayout)

    def _initAlchemicalTab(self):
        """initialize alchemical tab of post-treatment
        """

        self.alchemicalTab = QWidget()
        self.alchemicalTabLayout = QVBoxLayout()

        
        self.restraintInputs = QGroupBox('Inputs for alchemical simulations (.fepout/.log):')
        self.restraintInputsLayout = QVBoxLayout()

        # bound state
        self.alchemicalBoundStateLabel = QLabel('Atoms/Bound state (.fepout):')
        self.alchemicalBoundStateLayout = QGridLayout()
        
        self.alchemicalForwardLabel1 = QLabel('Forward:')
        self.alchemicalForwardLineEdit1 = QLineEdit()
        self.alchemicalForwardButton1 = QPushButton('Browse')
        self.alchemicalBoundStateLayout.addWidget(self.alchemicalForwardLabel1, 0, 0)
        self.alchemicalBoundStateLayout.addWidget(self.alchemicalForwardLineEdit1, 0, 1)
        self.alchemicalBoundStateLayout.addWidget(self.alchemicalForwardButton1, 0, 2)

        self.alchemicalBackwardLabel1 = QLabel('Backward:')
        self.alchemicalBackwardLineEdit1 = QLineEdit()
        self.alchemicalBackwardButton1 = QPushButton('Browse')
        self.alchemicalBoundStateLayout.addWidget(self.alchemicalBackwardLabel1, 1, 0)
        self.alchemicalBoundStateLayout.addWidget(self.alchemicalBackwardLineEdit1, 1, 1)
        self.alchemicalBoundStateLayout.addWidget(self.alchemicalBackwardButton1, 1, 2)

        self.alchemicalBoundStateLabel2 = QLabel('Restraints/Bound state (.log):')
        self.alchemicalBoundStateLayout2 = QGridLayout()
        
        self.alchemicalForwardLabel2 = QLabel('Forward:')
        self.alchemicalForwardLineEdit2 = QLineEdit()
        self.alchemicalForwardButton2 = QPushButton('Browse')
        self.alchemicalBoundStateLayout2.addWidget(self.alchemicalForwardLabel2, 0, 0)
        self.alchemicalBoundStateLayout2.addWidget(self.alchemicalForwardLineEdit2, 0, 1)
        self.alchemicalBoundStateLayout2.addWidget(self.alchemicalForwardButton2, 0, 2)

        self.alchemicalBackwardLabel2 = QLabel('Backward:')
        self.alchemicalBackwardLineEdit2 = QLineEdit()
        self.alchemicalBackwardButton2 = QPushButton('Browse')
        self.alchemicalBoundStateLayout2.addWidget(self.alchemicalBackwardLabel2, 1, 0)
        self.alchemicalBoundStateLayout2.addWidget(self.alchemicalBackwardLineEdit2, 1, 1)
        self.alchemicalBoundStateLayout2.addWidget(self.alchemicalBackwardButton2, 1, 2)

        # unbound state

        self.alchemicalUnboundStateLabel = QLabel('Atoms/Unbound state (.fepout):')
        self.alchemicalUnboundStateLayout = QGridLayout()
        
        self.alchemicalForwardLabel3 = QLabel('Forward:')
        self.alchemicalForwardLineEdit3 = QLineEdit()
        self.alchemicalForwardButton3 = QPushButton('Browse')
        self.alchemicalUnboundStateLayout.addWidget(self.alchemicalForwardLabel3, 0, 0)
        self.alchemicalUnboundStateLayout.addWidget(self.alchemicalForwardLineEdit3, 0, 1)
        self.alchemicalUnboundStateLayout.addWidget(self.alchemicalForwardButton3, 0, 2)

        self.alchemicalBackwardLabel3 = QLabel('Backward:')
        self.alchemicalBackwardLineEdit3 = QLineEdit()
        self.alchemicalBackwardButton3 = QPushButton('Browse')
        self.alchemicalUnboundStateLayout.addWidget(self.alchemicalBackwardLabel3, 1, 0)
        self.alchemicalUnboundStateLayout.addWidget(self.alchemicalBackwardLineEdit3, 1, 1)
        self.alchemicalUnboundStateLayout.addWidget(self.alchemicalBackwardButton3, 1, 2)

        self.alchemicalUnboundStateLabel2 = QLabel('Restraints/Unbound state (.log):')
        self.alchemicalUnboundStateLayout2 = QGridLayout()
        
        self.alchemicalForwardLabel4 = QLabel('Forward:')
        self.alchemicalForwardLineEdit4 = QLineEdit()
        self.alchemicalForwardButton4 = QPushButton('Browse')
        self.alchemicalUnboundStateLayout2.addWidget(self.alchemicalForwardLabel4, 0, 0)
        self.alchemicalUnboundStateLayout2.addWidget(self.alchemicalForwardLineEdit4, 0, 1)
        self.alchemicalUnboundStateLayout2.addWidget(self.alchemicalForwardButton4, 0, 2)

        self.alchemicalBackwardLabel4 = QLabel('Backward:')
        self.alchemicalBackwardLineEdit4 = QLineEdit()
        self.alchemicalBackwardButton4 = QPushButton('Browse')
        self.alchemicalUnboundStateLayout2.addWidget(self.alchemicalBackwardLabel4, 1, 0)
        self.alchemicalUnboundStateLayout2.addWidget(self.alchemicalBackwardLineEdit4, 1, 1)
        self.alchemicalUnboundStateLayout2.addWidget(self.alchemicalBackwardButton4, 1, 2)

        self.restraintInputsLayout.addWidget(self.alchemicalBoundStateLabel)
        self.restraintInputsLayout.addLayout(self.alchemicalBoundStateLayout)
        self.restraintInputsLayout.addWidget(self.alchemicalBoundStateLabel2)
        self.restraintInputsLayout.addLayout(self.alchemicalBoundStateLayout2)
        self.restraintInputsLayout.addWidget(self.alchemicalUnboundStateLabel)
        self.restraintInputsLayout.addLayout(self.alchemicalUnboundStateLayout)
        self.restraintInputsLayout.addWidget(self.alchemicalUnboundStateLabel2)
        self.restraintInputsLayout.addLayout(self.alchemicalUnboundStateLayout2)
        self.restraintInputs.setLayout(self.restraintInputsLayout)

        # alchemical force constants
        self.alchemicalForceConstants = QGroupBox('Force constants (in Colvars unit):')

        # all widgets
        self.alchemicalFCLayout = QHBoxLayout()
        self.alchemicalfcThetaLabel = QLabel('Theta:')
        self.alchemicalfcThetaLineEdit = QLineEdit('0.1')
        self.alchemicalfcPhiLabel = QLabel(' Phi: ')
        self.alchemicalfcPhiLineEdit = QLineEdit('0.1')
        self.alchemicalfcPsiLabel = QLabel('Psi:  ')
        self.alchemicalfcPsiLineEdit = QLineEdit('0.1')
        self.alchemicalfcthetaLabel = QLabel('theta:')
        self.alchemicalfcthetaLineEdit = QLineEdit('0.1')
        self.alchemicalfcphiLabel = QLabel(' phi: ')
        self.alchemicalfcphiLineEdit = QLineEdit('0.1')
        self.alchemicalfcRLabel = QLabel('r: ')
        self.alchemicalfcRLineEdit = QLineEdit('10')

        self.alchemicalFCLayout.addWidget(self.alchemicalfcThetaLabel)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcThetaLineEdit)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcPhiLabel)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcPhiLineEdit)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcPsiLabel)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcPsiLineEdit)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcthetaLabel)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcthetaLineEdit)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcphiLabel)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcphiLineEdit)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcRLabel)
        self.alchemicalFCLayout.addWidget(self.alchemicalfcRLineEdit)

        self.alchemicalForceConstants.setLayout(self.alchemicalFCLayout)

        # alchemical restraint centers
        self.alchemicalRestraintCenters = QGroupBox('Restraint centers (in Colvars unit), temperature and others:')

        # all widgets
        self.alchemicalRCLayout = QHBoxLayout()
        self.alchemicalRCThetaLabel = QLabel('Theta:')
        self.alchemicalRCThetaLineEdit = QLineEdit('0')
        self.alchemicalRCthetaLabel = QLabel('theta:')
        self.alchemicalRCthetaLineEdit = QLineEdit('90')
        self.alchemicalRCRLabel = QLabel('r: ')
        self.alchemicalRCRLineEdit = QLineEdit('8')
        self.alchemicalPostTemperatureLabel = QLabel('temperature:')
        self.alchemicalPostTemperatureLineEdit = QLineEdit('300')
        self.alchemicalPostTypeLabel = QLabel('Post-treatment type:')
        self.alchemicalPostTypeBox = QComboBox()
        self.alchemicalPostTypeBox.addItem('FEP')
        self.alchemicalPostTypeBox.addItem('BAR')

        self.alchemicalRCLayout.addWidget(self.alchemicalRCThetaLabel)
        self.alchemicalRCLayout.addWidget(self.alchemicalRCThetaLineEdit)
        self.alchemicalRCLayout.addWidget(self.alchemicalRCthetaLabel)
        self.alchemicalRCLayout.addWidget(self.alchemicalRCthetaLineEdit)
        self.alchemicalRCLayout.addWidget(self.alchemicalRCRLabel)
        self.alchemicalRCLayout.addWidget(self.alchemicalRCRLineEdit)
        self.alchemicalRCLayout.addWidget(self.alchemicalPostTemperatureLabel)
        self.alchemicalRCLayout.addWidget(self.alchemicalPostTemperatureLineEdit)
        self.alchemicalRCLayout.addWidget(self.alchemicalPostTypeLabel)
        self.alchemicalRCLayout.addWidget(self.alchemicalPostTypeBox)

        self.alchemicalRestraintCenters.setLayout(self.alchemicalRCLayout)

        self.alchemicalTabLayout.addWidget(self.restraintInputs)
        self.alchemicalTabLayout.addWidget(self.alchemicalForceConstants)
        self.alchemicalTabLayout.addWidget(self.alchemicalRestraintCenters)
        self.alchemicalTab.setLayout(self.alchemicalTabLayout)

    def _initQuickPlotTab(self):
        """initialize quick-plot tab
        """

        self.quickPlot = QWidget()
        self.quickPlotLayout = QVBoxLayout()

        # plot a (stratified) pmf
        self.plotPmf = QGroupBox('Plot (stratified) PMFs:')
        self.plotPmfLayout = QVBoxLayout()

        self.plotPmfLabel = QLabel('PMF files:')
        self.plotPmfBox = QListWidget()
        self.plotPmfChildLayout = QHBoxLayout()
        self.plotPmfAddButton = QPushButton('Add')
        self.plotPmfClearButton = QPushButton('Clear')
        self.plotPmfPlotButton = QPushButton('Plot')
        self.plotPmfChildLayout.addWidget(self.plotPmfAddButton)
        self.plotPmfChildLayout.addWidget(self.plotPmfClearButton)
        self.plotPmfChildLayout.addWidget(self.plotPmfPlotButton)

        self.plotPmfLayout.addWidget(self.plotPmfLabel)
        self.plotPmfLayout.addWidget(self.plotPmfBox)
        self.plotPmfLayout.addLayout(self.plotPmfChildLayout)
        self.plotPmf.setLayout(self.plotPmfLayout)

        # calculate pmf RMSD convergence
        self.plotPmfConvergence = QGroupBox('Calculate PMF RMSD convergence:')
        self.plotPmfConvergenceLayout = QVBoxLayout()

        self.plotPmfConvergenceLabel = QLabel('history file:')
        self.plotPmfConvergenceBox = QLineEdit()
        self.plotPmfConvergenceChildLayout = QHBoxLayout()
        self.plotPmfConvergenceBrowseButton = QPushButton('Browse')
        self.plotPmfConvergencePlotButton = QPushButton('Plot')
        self.plotPmfConvergenceChildLayout.addWidget(self.plotPmfConvergenceBrowseButton)
        self.plotPmfConvergenceChildLayout.addWidget(self.plotPmfConvergencePlotButton)

        self.plotPmfConvergenceLayout.addWidget(self.plotPmfConvergenceLabel)
        self.plotPmfConvergenceLayout.addWidget(self.plotPmfConvergenceBox)
        self.plotPmfConvergenceLayout.addLayout(self.plotPmfConvergenceChildLayout)
        self.plotPmfConvergence.setLayout(self.plotPmfConvergenceLayout)

        # merge a (stratified) pmf
        self.mergePmf = QGroupBox('Merge (stratified) PMFs:')
        self.mergePmfLayout = QVBoxLayout()

        self.mergePmfLabel = QLabel('PMF files:')
        self.mergePmfBox = QListWidget()
        self.mergePmfChildLayout = QHBoxLayout()
        self.mergePmfAddButton = QPushButton('Add')
        self.mergePmfClearButton = QPushButton('Clear')
        self.mergePmfmergeButton = QPushButton('Merge')
        self.mergePmfChildLayout.addWidget(self.mergePmfAddButton)
        self.mergePmfChildLayout.addWidget(self.mergePmfClearButton)
        self.mergePmfChildLayout.addWidget(self.mergePmfmergeButton)

        self.mergePmfLayout.addWidget(self.mergePmfLabel)
        self.mergePmfLayout.addWidget(self.mergePmfBox)
        self.mergePmfLayout.addLayout(self.mergePmfChildLayout)
        self.mergePmf.setLayout(self.mergePmfLayout)

        # plot hysteresis between forward and backward simulations
        self.plotHysteresis = QGroupBox('Plot hysteresis between bidirectional simulations:')
        self.plotHysteresisLayout = QVBoxLayout()

        self.plotHysteresisForwardLayout = QHBoxLayout()
        self.plotHysteresisForwardLabel = QLabel('Forward (fepout/log): ')
        self.plotHysteresisForwardLineEdit = QLineEdit()
        self.plotHysteresisForwardButton = QPushButton('Browse')
        self.plotHysteresisForwardLayout.addWidget(self.plotHysteresisForwardLabel)
        self.plotHysteresisForwardLayout.addWidget(self.plotHysteresisForwardLineEdit)
        self.plotHysteresisForwardLayout.addWidget(self.plotHysteresisForwardButton)

        self.plotHysteresisBackwardLayout = QHBoxLayout()
        self.plotHysteresisBackwardLabel = QLabel('Backward (fepout/log):')
        self.plotHysteresisBackwardLineEdit = QLineEdit()
        self.plotHysteresisBackwardButton = QPushButton('Browse')
        self.plotHysteresisBackwardLayout.addWidget(self.plotHysteresisBackwardLabel)
        self.plotHysteresisBackwardLayout.addWidget(self.plotHysteresisBackwardLineEdit)
        self.plotHysteresisBackwardLayout.addWidget(self.plotHysteresisBackwardButton)

        self.plotHysteresisPlotButton = QPushButton('Plot')

        self.plotHysteresisLayout.addLayout(self.plotHysteresisForwardLayout)
        self.plotHysteresisLayout.addLayout(self.plotHysteresisBackwardLayout)
        self.plotHysteresisLayout.addWidget(self.plotHysteresisPlotButton)

        self.plotHysteresis.setLayout(self.plotHysteresisLayout)

        self.quickPlotLayout.addWidget(self.plotPmf)
        self.quickPlotLayout.addWidget(self.mergePmf)
        self.quickPlotLayout.addWidget(self.plotPmfConvergence)
        self.quickPlotLayout.addWidget(self.plotHysteresis)
        self.quickPlot.setLayout(self.quickPlotLayout)

    # slots are defined below
    # otherwise they are defined in slots.py
    def _mainSettings(self):
        """call main settings
        """
        def f():
            self.mainSettings.show()
        return f

    def _advancedSettings(self, comboBox):
        """call advanced settings in pre-treatment
           the returned function is depended on the comboBox(strategy type)
        """

        def f():
            if comboBox.currentText() == 'Geometric':
                self.geometricAdvancedSettings.show()
            elif comboBox.currentText() == 'Alchemical':
                self.alchemicalAdvancedSettings.show()

        return f

    def _changeFFButtonState(self):
        """enable/disable the add and clear button of force field section
        """

        if self.forceFieldCombobox.currentText() == 'CHARMM':
            self.forceFieldAddButton.setEnabled(True)
            self.forceFieldClearButton.setEnabled(True)
            self.forceFieldFilesBox.setEnabled(True)
        elif self.forceFieldCombobox.currentText() == 'Amber':
            self.forceFieldAddButton.setEnabled(False)
            self.forceFieldClearButton.setEnabled(False)
            self.forceFieldFilesBox.setEnabled(False)
            
    def _changeStrategySettingState(self):
        """enable/disable the alchemical route and some advanced settings based on the MD engine
        """
        
        if self.selectMDEngineCombobox.currentText() == 'NAMD':
            self.selectStrategyCombobox.setEnabled(True)
            self.geometricAdvancedSettings.useOldCvCheckbox.setEnabled(True)
            self.geometricAdvancedSettings.OPLSMixingRuleCheckbox.setEnabled(True)
            self.geometricAdvancedSettings.useGaWTMCheckbox.setEnabled(True)
            
        elif self.selectMDEngineCombobox.currentText() == 'Gromacs':
            index = self.selectStrategyCombobox.findText('Geometric', QtCore.Qt.MatchFixedString)
            if index >= 0:
                self.selectStrategyCombobox.setCurrentIndex(index)
            self.selectStrategyCombobox.setEnabled(False)

            self.geometricAdvancedSettings.useOldCvCheckbox.setChecked(False)
            self.geometricAdvancedSettings.useOldCvCheckbox.setEnabled(False)

            self.geometricAdvancedSettings.OPLSMixingRuleCheckbox.setChecked(False)
            self.geometricAdvancedSettings.OPLSMixingRuleCheckbox.setEnabled(False)

            self.geometricAdvancedSettings.useGaWTMCheckbox.setEnabled(False)
            self.geometricAdvancedSettings.useGaWTMCheckbox.setChecked(False)
            
    def _changeStrategySettingStateForOldGromacs(self):
        """enable/disable a lot of options for the old Gromacs tab
        """
        
        if self.preTreatmentMainTabs.currentWidget() == self.NAMDTab:
            self.selectStrategyAdvancedButton.setEnabled(True)
            self.selectMDEngineCombobox.setEnabled(True)
            
        elif self.preTreatmentMainTabs.currentWidget() == self.GromacsTab:
            index = self.selectMDEngineCombobox.findText('Gromacs', QtCore.Qt.MatchFixedString)
            if index >= 0:
                self.selectMDEngineCombobox.setCurrentIndex(index)
            
            index = self.selectStrategyCombobox.findText('Geometric', QtCore.Qt.MatchFixedString)
            if index >= 0:
                self.selectStrategyCombobox.setCurrentIndex(index)
            
            self.selectMDEngineCombobox.setEnabled(False)
            self.selectStrategyCombobox.setEnabled(False)
            self.selectStrategyAdvancedButton.setEnabled(False)

    def _openDocFile(self):
        """open Documentation file
        """

        webbrowser.open('https://www.nature.com/articles/s41596-021-00676-1')

    def _openPythonAPIFile(self):
        """open Python API Documentation file
        """

        webbrowser.open('https://fhh2626.github.io/BFEE2APIDocs/')

    def _showAboutBox(self):
        """the about message box
        
        Returns:
            function obj: a slot function that shows that about message box
        """

        def f():
            QMessageBox.about(
                self,
                'About',
                f'<center><b>{__PROGRAM_NAME__}</b></center><br>'+r'''
                Binding free energy estimator (BFEE) is a python-based software
                that automates absolute binding free energy calculations through
                either the alchemical or geometric route by molecular dynamics
                simulations.<br>
                <b>Authors:</b><br>
                Haohao Fu (<a href="mailto:fhh2626@gmail.com">fhh2626@gmail.com</a>)<br>
                Haochuan Chen (<a href="mailto:summersnow9403@gmail.com">summersnow9403@gmail.com</a>)<br>
                <b>License:</b><br>
                BFEE2 is free software: you can redistribute it and/or modify it
                under the terms of the GNU General Public License as published by
                the Free Software Foundation, either version 3 of the License, or
                (at your option) any later version.<br>
                <b>Reference:</b><br>
                <b>BFEE2:</b> Fu et al. Nat. Protoc. 2022, 17, 1114-1141 and Fu et al. J. Chem. Inf. Model. 2021, 61, 2116–2123<br>
                <b>Alchemical and geometric routes:</b> Gumbart et al. J. Chem. Theory Comput. 2013, 9, 794–802<br>
                <b>WTM-eABF:</b> Fu et al. Acc. Chem. Res. 2019, 52, 3254–3264 and Fu et al. J. Phys. Chem. Lett. 2018, 9, 4738–4745<br>
                <b>Collective variables:</b> Fu et al. J. Chem. Theory Comput. 2017, 13, 5173–5178<br>
                <b>Contact</b> Wensheng Cai (<a href="mailto:wscai@nankai.edu.cn">wscai@nankai.edu.cn</a>)
                and Chris Chipot (<a href="mailto:chipot@ks.uiuc.edu">chipot@ks.uiuc.edu</a>)
                for further copyright information.
                ''')
        return f

    def _showGeometricResults(self, unit):
        ''' calculate binding from the geometric route,
            parameters in the Geometric tab will be read.
            Show a QMessageBox for the result
            Inputs:
                unit (string): 'namd' or 'gromacs' '''
        
        pTreat = postTreatment.postTreatment(
            float(self.postTemperatureLineEdit.text()), unit, 'geometric')
            
        pmfs = [
                    self.rmsdBoundLineEdit.text(), 
                    self.ThetaLineEdit.text(), 
                    self.PhiLineEdit.text(), 
                    self.PsiLineEdit.text(), 
                    self.thetaLineEdit.text(), 
                    self.phiLineEdit.text(), 
                    self.rLineEdit.text(), 
                    self.rmsdUnboundLineEdit.text()
        ]
            
        try:
            parameters = [
                    float(self.fcRMSDLineEdit.text()), 
                    float(self.fcThetaLineEdit.text()), 
                    float(self.fcPhiLineEdit.text()), 
                    float(self.fcPsiLineEdit.text()),
                    float(self.fcthetaLineEdit.text()), 
                    float(self.fcphiLineEdit.text()), 
                    float(self.postRstarLineEdit.text()), 
                    float(self.fcRMSDLineEdit.text())
            ]
        except:
            QMessageBox.warning(self, 'Error', f'Force constant or r* input error!')
            return

        # check inputs
        for index, item in enumerate(pmfs):
            if (index != 0) and (index != 7) and (not os.path.exists(item)):
                QMessageBox.warning(self, 'Error', f'file {item} does not exist!')
                return
            if (index == 7) and (not os.path.exists(item)):
                QMessageBox.warning(self, 'Warning', f'Supposing dealing with a rigid ligand!')

        # calculate free energies
        try:
            result = pTreat.geometricBindingFreeEnergy(pmfs, parameters)
        except postTreatment.RStarTooLargeError:
            QMessageBox.warning(
                self, 
                'Error', 
                f'\
r_star cannot be larger than r_max of step 7!\n'
            )
            return
        except Exception as e:
            print(e)
            QMessageBox.warning(
                self, 
                'Error', 
                f'\
Unknown error! The error message is: \n\
{e}\n'
            )
            return

        QMessageBox.about(
            self,
            'Result',
            f'\
Results:\n\
ΔG(site,c)            = {result[0]:.2f} kcal/mol\n\
ΔG(site,eulerTheta)   = {result[1]:.2f} kcal/mol\n\
ΔG(site,eulerPhi)     = {result[2]:.2f} kcal/mol\n\
ΔG(site,eulerPsi)     = {result[3]:.2f} kcal/mol\n\
ΔG(site,polarTheta)   = {result[4]:.2f} kcal/mol\n\
ΔG(site,polarPhi)     = {result[5]:.2f} kcal/mol\n\
(1/beta)*ln(S*I*C0)   = {result[6]:.2f} kcal/mol\n\
ΔG(bulk,c)            = {result[7]:.2f} kcal/mol\n\
ΔG(bulk,o)            = {result[8]:.2f} kcal/mol\n\
\n\
Standard Binding Free Energy:\n\
ΔG(total)             = {result[9]:.2f} kcal/mol\n'
        )

    def _showAlchemicalResults(self, unit):
        """calculate binding from the alchemical route,
            parameters in the alchemical tab will be read.
            Show a QMessageBox for the result

        Args:
            unit (str): 'namd' or 'gromacs'
        """
        
        pTreat = postTreatment.postTreatment(
            float(self.alchemicalPostTemperatureLineEdit.text()), unit, 'geometric')
        
        # alchemical outputs
        files = [
                    self.alchemicalForwardLineEdit1.text(), 
                    self.alchemicalBackwardLineEdit1.text(), 
                    self.alchemicalForwardLineEdit2.text(), 
                    self.alchemicalBackwardLineEdit2.text(), 
                    self.alchemicalForwardLineEdit3.text(), 
                    self.alchemicalBackwardLineEdit3.text(), 
                    self.alchemicalForwardLineEdit4.text(), 
                    self.alchemicalBackwardLineEdit4.text()
        ]
            
        try:
            parameters = [
                    float(self.alchemicalRCThetaLineEdit.text()), 
                    float(self.alchemicalRCthetaLineEdit.text()), 
                    float(self.alchemicalRCRLineEdit.text()), 
                    float(self.alchemicalfcThetaLineEdit.text()),
                    float(self.alchemicalfcPhiLineEdit.text()), 
                    float(self.alchemicalfcPsiLineEdit.text()), 
                    float(self.alchemicalfcthetaLineEdit.text()), 
                    float(self.alchemicalfcphiLineEdit.text()),
                    float(self.alchemicalfcRLineEdit.text())
            ]
            temperature = float(self.alchemicalPostTemperatureLineEdit.text())
        except:
            QMessageBox.warning(self, 'Error', f'Force constant or restraint center or temperature input error!')
            return
        if self.alchemicalPostTypeBox.currentText() == 'FEP':
                jobType = 'fep'
        elif self.alchemicalPostTypeBox.currentText() == 'BAR':
                jobType = 'bar'

        # check inputs
        for item in [
                    self.alchemicalBackwardLineEdit1.text(), 
                    self.alchemicalBackwardLineEdit2.text(), 
                    self.alchemicalBackwardLineEdit3.text(), 
                    self.alchemicalBackwardLineEdit4.text()
        ]:
            if (not os.path.exists(item)) and item != '':
                QMessageBox.warning(self, 'Error', f'backward file {item} does not exist and is not empty!')
                return

        for item in [
                    self.alchemicalForwardLineEdit1.text(), 
                    self.alchemicalForwardLineEdit2.text(), 
                    self.alchemicalForwardLineEdit3.text(), 
                    self.alchemicalForwardLineEdit4.text(), 
        ]:
            if not os.path.exists(item):
                QMessageBox.warning(self, 'Error', f'file {item} does not exist!')
                return

        # calculate free energies
        result, errors = pTreat.alchemicalBindingFreeEnergy(files, parameters, temperature, jobType)

        QMessageBox.about(
            self,
            'Result',
            f'\
Results:\n\
ΔG(site,couple)   = {result[0]:.2f} ± {errors[0]:.2f} kcal/mol\n\
ΔG(site,c+o+a+r)  = {result[1]:.2f} ± {errors[1]:.2f} kcal/mol\n\
ΔG(bulk,decouple) = {result[2]:.2f} ± {errors[2]:.2f} kcal/mol\n\
ΔG(bulk,c)        = {result[3]:.2f} ± {errors[3]:.2f} kcal/mol\n\
ΔG(bulk,o+a+r)    = {result[4]:.2f} kcal/mol\n\
\n\
Standard Binding Free Energy:\n\
ΔG(total)         = {result[5]:.2f} ± {errors[5]:.2f} kcal/mol\n'
        )

    def _showFinalResults(self):
        """calculate binding free energy and show the final results
        
        Returns:
            function obj: a slot function that calculates binding free energy and show the final results
        """

        def f():
            if self.postPMFTypeBox.currentText() == 'NAMD':
                unit = 'namd'
            elif self.postPMFTypeBox.currentText() == 'Gromacs':
                unit = 'gromacs'

            if self.postTreatmentMainTabs.currentIndex() == 0:
                jobType = 'geometric'
            elif self.postTreatmentMainTabs.currentIndex() == 1:
                jobType = 'alchemical'

            if jobType == 'geometric':
                self._showGeometricResults(unit)
            elif jobType == 'alchemical':
                self._showAlchemicalResults(unit)
            
        return f

    def _generateInputFiles(self):
        """generate input files for binding free energy simulation
        
        Returens:
            function obj: a slot function that generates input files for binding free energy simulation
        """
        
        def f():
            path = QFileDialog.getExistingDirectory(None, 'Select a directory')
            # cancel
            if path == '':
                return

            iGenerator = inputGenerator.inputGenerator()

            # third-party softwares and user-provided solvation boxes
            for item in [
                self.mainSettings.vmdLineEdit.text(),
                #self.mainSettings.gromacsLineEdit.text(),
                #self.mainSettings.tleapLineEdit.text(),
                self.geometricAdvancedSettings.nonStandardSolventPdbLineEdit.text(),
                self.geometricAdvancedSettings.nonStandardSolventPsfLineEdit.text()
            ]:
                if ((not os.path.exists(item)) and item != ''):
                    QMessageBox.warning(self, 'Error', f'file {item} does not exist!')
                    return

            if self.preTreatmentMainTabs.currentIndex() == 0:
                # force fields
                forceFieldFiles = [self.forceFieldFilesBox.item(i).text() for i in range(self.forceFieldFilesBox.count())]

                # NAMD files
                for item in [self.psfLineEdit.text(), self.coorLineEdit.text()] + forceFieldFiles:
                    if not os.path.exists(item):
                        QMessageBox.warning(self, 'Error', f'file {item} does not exist!')
                        return

                # check inputs
                try:
                    float(self.temperatureLineEdit.text())
                    stratification = [
                        int(self.geometricAdvancedSettings.stratificationRMSDBoundLineEdit.text()),
                        int(self.geometricAdvancedSettings.stratificationThetaLineEdit.text()),
                        int(self.geometricAdvancedSettings.stratificationPhiLineEdit.text()),
                        int(self.geometricAdvancedSettings.stratificationPsiLineEdit.text()),
                        int(self.geometricAdvancedSettings.stratificationthetaLineEdit.text()),
                        int(self.geometricAdvancedSettings.stratificationphiLineEdit.text()),
                        int(self.geometricAdvancedSettings.stratificationRLineEdit.text()),
                        int(self.geometricAdvancedSettings.stratificationRMSDUnboundLineEdit.text())
                    ]
                    alchemicalStratification = [
                        int(self.alchemicalAdvancedSettings.boundLigandLineEdit.text()),
                        int(self.alchemicalAdvancedSettings.boundRestraintsLineEdit.text()),
                        int(self.alchemicalAdvancedSettings.unboundLigandLineEdit.text()),
                        int(self.alchemicalAdvancedSettings.unboundRestraintsLineEdit.text())
                    ]
                except:
                    QMessageBox.warning(self, 'Error', f'Force constant or r* input error!')
                    return

                # job type
                if self.forceFieldCombobox.currentText() == 'CHARMM':
                    forceFieldType = 'charmm'
                elif self.forceFieldCombobox.currentText() == 'Amber':
                    forceFieldType = 'amber'

                # make sure there are CHARMM FF files
                if forceFieldFiles == [] and forceFieldType == 'charmm':
                    QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
CHARMM force field files must be specified!'
                        )
                    return
                
                if self.selectStrategyCombobox.currentText() == 'Geometric':
                    
                    if self.selectMDEngineCombobox.currentText().lower() == 'gromacs':
                        QMessageBox.warning(
                            self, 'Warning', f'Please pay attention if you use the new Gromacs interface: \n\
1. Occationally, the atom number in complex.ndx do not correct. If so, just modify it manually; \n\
2. Make sure in complex.gro, the center of box is approximately at (x/2, y/2, z/2) rather than (0, 0, 0);  \n\
3. Gromacs patched by the latest (master branch) version of Colvars is needed.'
                )

                    # for the amber force field, files of large box must be provided
                    if forceFieldType == 'amber':

                        if self.geometricAdvancedSettings.nonStandardSolventPdbLineEdit.text() == '' or \
                            self.geometricAdvancedSettings.nonStandardSolventPsfLineEdit.text() == '':
                            QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Coordinate and topology files of large box must be \
provided in "Advanced Settings"when using the Amber \
force fields!'
                            )
                            return
                    
                    if self.geometricAdvancedSettings.useGaWTMCheckbox.isChecked():
                        QMessageBox.warning(self, 'Error', 
                                f'\
The feature of using GaWTM-eABF as the workhorse engine is \
experimental! Known issues:\n \
The config files for extending GaWTM-eABF simulations \
will do pre-equilibration again. This may cause some errors. \
One can revise the config file manually to avoid it'
                        )

                    try:
                        iGenerator.generateNAMDGeometricFiles(
                            path,
                            self.psfLineEdit.text(),
                            self.coorLineEdit.text(),
                            forceFieldType,
                            forceFieldFiles,
                            float(self.temperatureLineEdit.text()),
                            self.selectProteineLineEdit.text(),
                            self.selectLigandLineEdit.text(),
                            self.geometricAdvancedSettings.userDefinedDirectionLineEdit.text(),
                            self.geometricAdvancedSettings.nonStandardSolventPsfLineEdit.text(),
                            self.geometricAdvancedSettings.nonStandardSolventPdbLineEdit.text(),
                            stratification,
                            self.geometricAdvancedSettings.memProCheckbox.isChecked(),
                            self.geometricAdvancedSettings.neutralizeLigOnlyCombobox.currentText(),
                            self.geometricAdvancedSettings.pinDownProCheckbox.isChecked(),
                            self.geometricAdvancedSettings.useOldCvCheckbox.isChecked(),
                            int(self.geometricAdvancedSettings.parallelRunsLineEdit.text()),
                            self.mainSettings.vmdLineEdit.text(),
                            self.geometricAdvancedSettings.reflectingBoundaryCheckbox.isChecked(),
                            self.selectMDEngineCombobox.currentText().lower(),
                            self.geometricAdvancedSettings.OPLSMixingRuleCheckbox.isChecked(),
                            self.geometricAdvancedSettings.considerRMSDCVCheckbox.isChecked(),
                            self.geometricAdvancedSettings.useGaWTMCheckbox.isChecked()
                        )
                    except fileParser.SelectionError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Selection corresponding to nothing!\n\
Check your selection again!'
                        )
                        if os.path.exists(f'{path}/BFEE'):
                            shutil.rmtree(f'{path}/BFEE')
                        return
                    except inputGenerator.DirectoryExistError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
./BFEE* directory already exists!'
                        )
                        return
                    except inputGenerator.FileTypeUnknownError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Unkwn input file types! The following file types are supported:\n\
psf, parm, prm7, parm7, prmtop, pdb, coor, rst, rst7, inpcrd '
                        )
                        return
                    except PermissionError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Cannot read input files due to the permission reason!\n\
Restart the program or check the authority of the files!'
                        )
                        return
                    except Exception as e:
                        print(e)
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Unknown error! The error message is: \n\
{e}\n'
                        )
                        return

                elif self.selectStrategyCombobox.currentText() == 'Alchemical':

                    try:
                        iGenerator.generateNAMDAlchemicalFiles(
                            path,
                            self.psfLineEdit.text(),
                            self.coorLineEdit.text(),
                            forceFieldType,
                            forceFieldFiles,
                            float(self.temperatureLineEdit.text()),
                            self.selectProteineLineEdit.text(),
                            self.selectLigandLineEdit.text(),
                            alchemicalStratification,
                            self.alchemicalAdvancedSettings.doubleWideCheckbox.isChecked(),
                            self.alchemicalAdvancedSettings.minBeforeSampleCheckbox.isChecked(),
                            self.alchemicalAdvancedSettings.memProCheckbox.isChecked(),
                            self.alchemicalAdvancedSettings.neutralizeLigOnlyCombobox.currentText(),
                            self.alchemicalAdvancedSettings.pinDownProCheckbox.isChecked(),
                            self.alchemicalAdvancedSettings.useOldCvCheckbox.isChecked(),
                            self.mainSettings.vmdLineEdit.text(),
                            self.alchemicalAdvancedSettings.OPLSMixingRuleCheckbox.isChecked(),
                        )
                    except PermissionError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Cannot read input files due to the permission reason!\n\
Restart the program or check the authority of the files!'
                        )
                        return
                    except fileParser.SelectionError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Selection corresponding to nothing!\n\
Check your selection again!'
                        )
                        if os.path.exists(f'{path}/BFEE'):
                            shutil.rmtree(f'{path}/BFEE')
                        return
                    except inputGenerator.DirectoryExistError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
./BFEE directory already exists!'
                        )
                        return
                    except inputGenerator.FileTypeUnknownError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Unkwn input file types! The following file types are supported:\n\
psf, parm, prm7, parm7, prmtop, pdb, coor, rst, rst7, inpcrd '
                        )
                        return
                    except Exception as e:
                        print(e)
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Unknown error! The error message is: \n\
{e}\n'
                        )
                        return

            # gromacs
            if self.preTreatmentMainTabs.currentIndex() == 1:

                QMessageBox.warning(
                    self, 'Warning', f'Any setting in "Advanced settings" is not supported \n\
when using Gromacs-formatted files as inputs!'
                )

                for item in [
                        self.topLineEdit.text(), 
                        self.gromacsPdbLineEdit.text(),
                        self.gromacsLigandOnlyPdbLineEdit.text(),
                        self.gromacsLigandOnlyTopLineEdit.text()
                    ]:
                    if not os.path.exists(item):
                        QMessageBox.warning(self, 'Error', f'file {item} does not exist!')
                        return

                if self.selectStrategyCombobox.currentText() == 'Geometric':
                    try:              
                        iGenerator.generateGromacsGeometricFiles(
                            path=path,
                            topFile=self.topLineEdit.text(),
                            pdbFile=self.gromacsPdbLineEdit.text(),
                            pdbFileFormat=self.gromacsStructureFileFormatCombobox.currentText(),
                            ligandOnlyTopFile=self.gromacsLigandOnlyTopLineEdit.text(),
                            ligandOnlyPdbFile=self.gromacsLigandOnlyPdbLineEdit.text(),
                            ligandOnlyPdbFileFormat=self.gromacsLigandOnlyStructureFileFormatCombobox.currentText(),
                            selectionPro=self.selectProteineLineEdit.text(),
                            selectionLig=self.selectLigandLineEdit.text(),
                            temperature=float(self.temperatureLineEdit.text())
                        )
                    except inputGenerator.DirectoryExistError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
./BFEE directory already exists!'
                        )
                        return
                    except BFEEGromacs.SelectionError:
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Selection corresponding to nothing!\n\
Check your selection again!'
                        )
                        if os.path.exists(f'{path}/BFEE'):
                            shutil.rmtree(f'{path}/BFEE')
                        return
                    except Exception as e:
                        print(e)
                        QMessageBox.warning(
                                self, 
                                'Error', 
                                f'\
Unknown error!'
                        )
                        return
                elif self.selectStrategyCombobox.currentText() == 'Alchemical':
                    QMessageBox.warning(self, 'Error', f'Alchemical route is not supported using Gromacs!')
                    return

            QMessageBox.information(self, 'Input generation', f'Input files have been generated successfully!')
            
            del iGenerator

        return f 

    def _mergePMFs(self):
        """merge a series of overlapped pmfs
        
        Returns:
            function obj: a slot function that merges a series of overlapped pmfs
        """
        
        def f():
            if self.mergePmfBox.count() == 0:
                QMessageBox.warning(self, 'Warning', f'Warning, no PMF selected!')
                return

            path, _ = QFileDialog.getSaveFileName(None, 'Set the name of merged PMF')

            pmfFiles = [self.mergePmfBox.item(i).text() for i in range(self.mergePmfBox.count())]

            if not ploter.isGaWTM(pmfFiles):
                pmfs = [ploter.readPMF(item) for item in pmfFiles]
            else:
                pmfFiles = [item for item in pmfFiles if not item.endswith('.reweightamd1.cumulant.pmf')]
                try:
                    pmfs = [ploter.correctGaWTM(item) for item in pmfFiles]
                except ploter.NoCorrectionFileError:
                    QMessageBox.warning(self, 'Error', f'A PMF file does not have corresponding correction!')
                    return

            mergedPMF = ploter.mergePMF(pmfs)
            ploter.writePMF(path, mergedPMF)
            QMessageBox.information(self, 'Merge PMFs', f'PMF merged successfully!')
        return f

    def _plotPMFs(self):
        """plot a series of overlapped pmfs
        
        Returns:
            function obj: a slot function that plots a series of overlapped pmfs
        """
        
        def f():
            if self.plotPmfBox.count() == 0:
                QMessageBox.warning(self, 'Warning', f'Warning, no PMF selected!')
                return

            pmfFiles = [self.plotPmfBox.item(i).text() for i in range(self.plotPmfBox.count())]

            if not ploter.isGaWTM(pmfFiles):
                pmfs = [ploter.readPMF(item) for item in pmfFiles]
            else:
                pmfFiles = [item for item in pmfFiles if not item.endswith('.reweightamd1.cumulant.pmf')]
                try:
                    pmfs = [ploter.correctGaWTM(item) for item in pmfFiles]
                except ploter.NoCorrectionFileError:
                    QMessageBox.warning(self, 'Error', f'A PMF file does not have corresponding correction!')
                    return
            
            mergedPMF = ploter.mergePMF(pmfs)
            ploter.plotPMF(mergedPMF)
        return f

    def _plotRMSDConvergence(self):
        """plot time evolution of PMF rmsd with respect to zero array
        
        Returns:
            function obj: a slot function that plots time evolution of PMF rmsd with respect to zero array
        """
        
        def f():
            path = self.plotPmfConvergenceBox.text()
            if not os.path.exists(path):
                QMessageBox.warning(self, 'Error', f'file {path} does not exist!')
                return
            
            rmsds = ploter.parseHistFile(path)
            ploter.plotConvergence(rmsds)
        return f
    
    def _plotHysteresis(self):
        """plot hysteresis between forward and backward alchemical transformations
        
        Returns:
            function obj: a slot function that plot hysteresis between forward and backward alchemical transformations
        """
        
        def f():
            forwardFilePath = self.plotHysteresisForwardLineEdit.text()
            backwardFilePath = self.plotHysteresisBackwardLineEdit.text() 
            if not os.path.exists(forwardFilePath):
                QMessageBox.warning(self, 'Error', f'file {forwardFilePath} does not exist!')
                return
            if not os.path.exists(backwardFilePath):
                QMessageBox.warning(self, 'Error', f'file {backwardFilePath} does not exist!')
                return

            forwardPostfix = os.path.splitext(forwardFilePath)[-1]
            backwardPostfix = os.path.splitext(backwardFilePath)[-1]

            if forwardPostfix != '.fepout' and forwardPostfix != '.log' \
                and backwardPostfix != '.fepout' and forwardPostfix != '.log':
                QMessageBox.warning(self, 'Error', f'File type not correct!')
                return

            if forwardPostfix != backwardPostfix:
                QMessageBox.warning(self, 'Error', f'File types of forward and backward simulations are not the same!')
                return
            
            pTreat = postTreatment.postTreatment(float(self.alchemicalPostTemperatureLineEdit.text()), 'namd', 'alchemical')
            if forwardPostfix == '.fepout':
                forwardProfile = np.transpose(pTreat._fepoutFile(forwardFilePath))
                backwardProfile = np.transpose(pTreat._fepoutFile(backwardFilePath))
                backwardProfile[:,1] *= -1
            elif forwardPostfix == '.log':
                forwardProfile = np.transpose(pTreat._tiLogFile(forwardFilePath))
                backwardProfile = np.transpose(pTreat._tiLogFile(backwardFilePath))
            ploter.plotHysteresis(forwardProfile, backwardProfile)
        return f


    def _initSingalsSlots(self):
        """initialize (connect) singals and slots
        """

        # pre-treatment tab
        self.selectStrategyAdvancedButton.clicked.connect(self._advancedSettings(self.selectStrategyCombobox))
        self.preTreatmentMainTabs.currentChanged.connect(self._changeStrategySettingStateForOldGromacs)

        # NAMD tab
        self.psfButton.clicked.connect(commonSlots.openFileDialog('psf/parm', self.psfLineEdit))
        self.coorButton.clicked.connect(commonSlots.openFileDialog('pdb/rst', self.coorLineEdit))

        # force field selection
        self.forceFieldCombobox.currentTextChanged.connect(self._changeFFButtonState)
        self.forceFieldAddButton.clicked.connect(commonSlots.openFilesDialog('prm', self.forceFieldFilesBox))
        self.forceFieldClearButton.clicked.connect(self.forceFieldFilesBox.clear)
        
        # MD engine
        self.selectMDEngineCombobox.currentTextChanged.connect(self._changeStrategySettingState)

        # gromacs tab
        self.gromacsPdbButton.clicked.connect(commonSlots.openFileDialog('pdb', self.gromacsPdbLineEdit))
        self.topButton.clicked.connect(commonSlots.openFileDialog('top', self.topLineEdit))
        self.gromacsLigandOnlyPdbButton.clicked.connect(commonSlots.openFileDialog('pdb', self.gromacsLigandOnlyPdbLineEdit))
        self.gromacsLigandOnlyTopButton.clicked.connect(commonSlots.openFileDialog('top', self.gromacsLigandOnlyTopLineEdit))

        # geometric tab
        self.rmsdBoundButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.rmsdBoundLineEdit))
        self.rmsdUnboundButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.rmsdUnboundLineEdit))
        self.ThetaButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.ThetaLineEdit))
        self.PhiButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.PhiLineEdit))
        self.PsiButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.PsiLineEdit))
        self.thetaButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.thetaLineEdit))
        self.phiButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.phiLineEdit))
        self.rButton.clicked.connect(commonSlots.openFileDialog('czar.pmf/UI.pmf', self.rLineEdit))

        # alchemical tab
        self.alchemicalForwardButton1.clicked.connect(commonSlots.openFileDialog('log', self.alchemicalForwardLineEdit1))
        self.alchemicalBackwardButton1.clicked.connect(commonSlots.openFileDialog('log', self.alchemicalBackwardLineEdit1))
        self.alchemicalForwardButton2.clicked.connect(commonSlots.openFileDialog('log', self.alchemicalForwardLineEdit2))
        self.alchemicalBackwardButton2.clicked.connect(commonSlots.openFileDialog('log', self.alchemicalBackwardLineEdit2))
        self.alchemicalForwardButton3.clicked.connect(commonSlots.openFileDialog('fepout', self.alchemicalForwardLineEdit3))
        self.alchemicalBackwardButton3.clicked.connect(commonSlots.openFileDialog('fepout', self.alchemicalBackwardLineEdit3))
        self.alchemicalForwardButton4.clicked.connect(commonSlots.openFileDialog('fepout', self.alchemicalForwardLineEdit4))
        self.alchemicalBackwardButton4.clicked.connect(commonSlots.openFileDialog('fepout', self.alchemicalBackwardLineEdit4))

        # generate input files
        self.generateInputButton.clicked.connect(self._generateInputFiles())

        # calculate binding free energy
        self.calculateButton.clicked.connect(self._showFinalResults())

        # quick-plot tab
        self.plotPmfAddButton.clicked.connect(commonSlots.openFilesDialog('pmf', self.plotPmfBox))
        self.plotPmfClearButton.clicked.connect(self.plotPmfBox.clear)
        self.plotPmfPlotButton.clicked.connect(self._plotPMFs())
        self.mergePmfAddButton.clicked.connect(commonSlots.openFilesDialog('pmf', self.mergePmfBox))
        self.mergePmfClearButton.clicked.connect(self.mergePmfBox.clear)
        self.mergePmfmergeButton.clicked.connect(self._mergePMFs())
        self.plotPmfConvergenceBrowseButton.clicked.connect(commonSlots.openFileDialog('pmf', self.plotPmfConvergenceBox))
        self.plotPmfConvergencePlotButton.clicked.connect(self._plotRMSDConvergence())
        self.plotHysteresisForwardButton.clicked.connect(commonSlots.openFileDialog('fepout/log', self.plotHysteresisForwardLineEdit))
        self.plotHysteresisBackwardButton.clicked.connect(commonSlots.openFileDialog('fepout/log', self.plotHysteresisBackwardLineEdit))
        self.plotHysteresisPlotButton.clicked.connect(self._plotHysteresis())
