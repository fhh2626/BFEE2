# generate NAMD/Gromacs/Colvars config files

class configTemplate:
    ''' generate NAMD/Gromacs/Colvars config files
        In the Colvars config file, ndx and xyz files are used to indicate the group of atoms
        The non-hydrogen atoms of protein are labeled as 'protein' in complex.ndx
        The non-hydrogen atoms of ligand are labeled as 'ligand' in complex.ndx
        The non-hydrogen atoms of user-defined reference are labeled as 'reference' in complex.ndx '''

    def __init__(self):
        pass

    def namdConfigTemplate(
                            self,
                            forceFieldType,
                            forceFieldFiles,
                            topFile,
                            coorFile,
                            NAMDRestartCoor,
                            NAMDRestartVel,
                            NAMDRestartXsc,
                            PBCCondition,
                            outputPrefix,
                            temperature,
                            numSteps,
                            cvFile = '',
                            cvDefinitionFile = '',
                            CVRestartFile = '',
                            fepFile = '',
                            fepWindowNum = 20,
                            fepForward = True,
                            fepDoubleWide = False,
                            fepMinBeforeSample = False,
                            membraneProtein = False
                            ):
        ''' the namd config file template
            Inputs:
                forceFieldType (string): 'charmm' or 'amber'
                forceFieldFiles (list of strings): name of charmm force field files
                topFile (string): name of the topology file (psf, parm) 
                coorFile (string): name of the coordinate file (pdb, rst)
                NAMDRestartCoor (string): name of namd binary restart coor file (if restart from a previous simulation)
                NAMDRestartVel (string): name of namd binary restart vel file
                NAMDRestartXsc (string): name of namd restart xsc file
                PBCCondition (np.array, 2*3): PBC vector, ((lengthX, lengthY, lengthZ),(centerX, centerY, centerZ))
                outputPrefix (string): prefix of output file
                temperature (float): temperature of the simulation
                numSteps (int): number of steps of the simulation
                cvFile (string): name of Colvars file
                cvDefinitionFile (string): name of a TCL file defining new CVs
                CVRestartFile (string): name of Colvars restart file
                fepFile (string): name of fep file, indicating which atoms will be generated/removed (if run alchemical simulation)
                fepWindowNum (int): number of fep windows
                fepForward (bool): whether this is a forward fep simulation
                fepDouble (bool): whether this is a double-wide fep simulation
                fepMinBeforeSample (bool): whether do minimization before sampling in each FEP window
                membraneProtein (bool): whether simulating a membrame protein
            Return:
                string: a NAMD config string if succeed, and empty string otherwise
            '''

        assert(forceFieldType == 'charmm' or forceFieldType == 'amber')

        configString = f'\
coordinates    {coorFile}                   \n'

        # force field files
        if forceFieldType == 'charmm':
            configString += f'\
structure      {topFile}                \n\
paraTypeCharmm    on                    \n'
            for ff in forceFieldFiles:
                configString += f'\
parameters    {ff}                      \n'
        elif forceFieldType == 'amber':
            configString += f'\
parmFile      {topFile}                 \n\
amber    on                             \n'
        else:
            # error
            return ''

        # structure
        if forceFieldType == 'charmm':
            configString += f'\
exclude    scaled1-4                    \n\
1-4scaling    1.0                       \n\
switching            on                 \n\
switchdist           10.0               \n\
cutoff               12.0               \n\
pairlistdist         14.0               \n'
        elif forceFieldType == 'amber':
            configString += f'\
exclude    scaled1-4                    \n\
1-4scaling    0.83333333                \n\
switching            on                 \n\
switchdist           8.0                \n\
cutoff               9.0                \n\
pairlistdist         11.0               \n'
        else:
            # error
            return ''

        if NAMDRestartCoor == '' and NAMDRestartVel == '' and NAMDRestartXsc == '' and PBCCondition != '':
            # set temperature
            configString += f'\
temperature    {temperature}                      \n\
cellBasisVector1    {PBCCondition[0][0]} 0 0         \n\
cellBasisVector2    0 {PBCCondition[0][1]} 0         \n\
cellBasisVector3    0 0 {PBCCondition[0][2]}         \n\
cellOrigin    {PBCCondition[1][0]} {PBCCondition[1][1]} {PBCCondition[1][2]}  \n'
        elif NAMDRestartCoor != '' and NAMDRestartVel != '' and NAMDRestartXsc != '' and PBCCondition == '':
            # restart from files
            configString += f'\
bincoordinates    {NAMDRestartCoor}                          \n\
binvelocities    {NAMDRestartVel}                            \n\
ExtendedSystem    {NAMDRestartXsc}                           \n'
        else:
            # error
            return ''

        # other parameters
        configString += f'\
binaryoutput         yes                        \n\
binaryrestart        yes                        \n\
outputname           {outputPrefix}             \n\
dcdUnitCell          yes                        \n\
outputenergies       5000                       \n\
outputtiming         5000                       \n\
outputpressure       5000                       \n\
restartfreq          5000                       \n\
XSTFreq              5000                       \n\
dcdFreq              5000                       \n\
hgroupcutoff         2.8                        \n\
wrapAll              off                        \n\
wrapWater            on                         \n\
langevin             on                         \n\
langevinDamping      1                          \n\
langevinTemp         300                        \n\
langevinHydrogen     no                         \n\
langevinpiston       on                         \n\
langevinpistontarget 1.01325                    \n\
langevinpistonperiod 200                        \n\
langevinpistondecay  100                        \n\
langevinpistontemp   {temperature}              \n\
usegrouppressure     yes                        \n\
PME                  yes                        \n\
PMETolerance         10e-6                      \n\
PMEInterpOrder       4                          \n\
PMEGridSpacing       1.0                        \n\
timestep             2.0                        \n\
fullelectfrequency   2                          \n\
nonbondedfreq        1                          \n\
rigidbonds           all                        \n\
rigidtolerance       0.00001                    \n\
rigiditerations      400                        \n\
stepspercycle        10                         \n\
splitpatch           hydrogen                   \n\
margin               2                          \n'

        # membrane protein
        if membraneProtein:
            configString += f'\
useflexiblecell      yes                        \n\
useConstantRatio     yes                        \n'
        else:
            configString += f'\
useflexiblecell      no                         \n\
useConstantRatio     no                         \n'

        # colvars definition
        if cvFile != '':
            configString += f'\
colvars    on                                   \n\
colvarsConfig    {cvFile}                       \n'

            if CVRestartFile != '':
                configString += f'\
colvarsInput     {CVRestartFile}                \n'


        if cvDefinitionFile != '':
            configString += f'\
source     {cvDefinitionFile}                   \n'

        # fep
        if fepFile == '':
            if NAMDRestartCoor == '' and NAMDRestartVel == '' and NAMDRestartXsc == '':
                configString += f'\
minimize    500                                 \n\
reinitvels    {temperature}                     \n'
            configString += f'\
run    {numSteps}                               \n'
        else:
            # currently the alchemical route is somewhat hard-coded
            # this will be improved in the future
            configString += f'\
source ../fep.tcl                                  \n\
alch on                                         \n\
alchType FEP                                    \n\
alchFile {fepFile}                              \n\
alchCol B                                       \n\
alchOutFile {outputPrefix}.fepout               \n\
alchOutFreq 50                                  \n\
alchVdwLambdaEnd 0.7                            \n\
alchElecLambdaStart 0.5                         \n\
alchEquilSteps 100000                           \n'

            if fepForward:
                if not fepDoubleWide:
                    if fepMinBeforeSample:
                        # minimize before sampling
                        configString += f'\
runFEPmin 0.0 1.0 {1.0/fepWindowNum} 500000 1000 {temperature}\n'
                    else:
                        configString += f'\
runFEP 0.0 1.0 {1.0/fepWindowNum} 500000\n'

                else:
                    # double wide simulation
                    configString += f'\
runFEP 0.0 1.0 {1.0/fepWindowNum} 500000 true\n'

            else:
                # backward
                if not fepDoubleWide:
                    if fepMinBeforeSample:
                        # minimize before sampling
                        configString += f'\
runFEPmin 1.0 0.0 {-1.0/fepWindowNum} 500000 1000 {temperature}\n'
                    else:
                        configString += f'\
runFEP 1.0 0.0 {-1.0/fepWindowNum} 500000\n'

                else:
                    # double wide simulation
                    configString += f'\
runFEP 1.0 0.0 {-1.0/fepWindowNum} 500000 true\n'

        return configString


    def cvRMSDTemplate(self, setBoundary, lowerBoundary, upperBoundary, refFile):
        ''' RMSD CV template
            Inputs:
                setBoundary (bool): whether set boundary (for free-energy calculation)
                lowerBoundary (float): lower boundary of free-energy calculaton
                upperboundary (float): upper boundary of free-energy calculation
                refFile (string): path to the reference file
            Return:
                string: string of RMSD definition '''

        string = f'\
colvar {{                                    \n\
    name RMSD                                \n'

        if setBoundary:
            string += f'\
    width 0.05                               \n\
    lowerboundary {lowerBoundary:.1f}            \n\
    upperboundary {upperBoundary:.1f}            \n\
    subtractAppliedForce on                  \n\
    expandboundaries  on                     \n\
    extendedLagrangian on                    \n\
    extendedFluctuation 0.05                 \n'

        string += f'\
    rmsd {{                                  \n\
        atoms {{                             \n\
            indexGroup  ligand               \n\
        }}                                   \n\
        refpositionsfile  {refFile}          \n\
    }}                                       \n\
}}                                           \n'
        return string

    def cvEulerAngleTemplate(self, setBoundary, lowerBoundary, upperBoundary, angle, refFile):
        ''' Euler angle template 
            Inputs:
                setBoundary (bool): whether set boundary (for free-energy calculation)
                lowerBoundary (float): lower boundary of free-energy calculaton
                upperboundary (float): upper boundary of free-energy calculation
                angle (string): 'eulerTheta', 'eulerPhi' of 'eulerPsi'
                refFile (string): path to the reference file
            Return:
                string: string of Euler angle definition '''

        assert(angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi')
        
        string = f'\
colvar {{                              \n\
    name {angle}                   \n'

        if angle == 'eulerTheta':
            string += f'\
    customFunction asin(2 * (q1*q3-q4*q2)) * 180 / 3.1415926\n'
        elif angle == 'eulerPhi':
            string += f'\
    customFunction atan2(2*(q1*q2+q3*q4), 1-2*(q2*q2+q3*q3)) * 180 / 3.1415926\n'
        elif angle == 'eulerPsi':
            string += f'\
    customFunction atan2(2*(q1*q4+q2*q3), 1-2*(q3*q3+q4*q4)) * 180 / 3.1415926\n'

        if setBoundary:
            string += f'\
    width 1                            \n\
    lowerboundary {lowerBoundary:.1f}      \n\
    upperboundary {upperBoundary:.1f}      \n\
    subtractAppliedForce on            \n\
    expandboundaries  on               \n\
    extendedLagrangian on              \n\
    extendedFluctuation 1              \n'

        string += f'\
    Orientation {{                             \n\
        name  q                                \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
	        enableFitGradients no              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
         refpositionsfile  {refFile}           \n\
    }}                                         \n\
}}                                             \n'

        return string

    def cvPolarAngleTemplate(self, setBoundary, lowerBoundary, upperBoundary, angle, refFile):
        ''' Polar angle template 
            Inputs:
                setBoundary (bool): whether set boundary (for free-energy calculation)
                lowerBoundary (float): lower boundary of free-energy calculaton
                upperboundary (float): upper boundary of free-energy calculation
                angle (string): 'polarTheta' or 'polarPhi'
                refFile (string): path to the reference file
            Return:
                string: string of polar angle definition '''

        assert(angle == 'polarTheta' or angle == 'polarPhi')
        
        string = f'\
colvar {{                                   \n\
    name {angle}                            \n'

        if angle == 'polarTheta':
            string += f'\
    customFunction acos(-i2) * 180 / 3.1415926\n'
        elif angle == 'polarPhi':
            string += f'\
    customFunction atan2(i3, i1) * 180 / 3.1415926\n'

        if setBoundary:
            string += f'\
    width 1                                 \n\
    lowerboundary {lowerBoundary:.1f}           \n\
    upperboundary {upperBoundary:.1f}           \n\
    subtractAppliedForce on                 \n\
    expandboundaries  on                    \n\
    extendedLagrangian on                   \n\
    extendedFluctuation 1                   \n'

        string += f'\
    distanceDir {{                          \n\
        name  i                             \n\
        group1 {{                           \n\
            indexGroup  reference           \n\
            centerReference    on           \n\
            rotateReference    on           \n\
            enableFitGradients no           \n\
            fittingGroup {{                 \n\
                indexGroup  protein         \n\
            }}                              \n\
            refpositionsfile  {refFile}     \n\
        }}                                  \n\
        group2 {{                           \n\
            indexGroup  ligand              \n\
            centerReference    on           \n\
            rotateReference    on           \n\
            enableFitGradients no           \n\
            fittingGroup {{                 \n\
                indexGroup  protein         \n\
            }}                              \n\
            refpositionsfile  {refFile}     \n\
        }}                                  \n\
    }}                                      \n\
}}                                          \n'
        return string

    def cvRTemplate(self, setBoundary, lowerBoundary, upperBoundary):
        ''' r distance template
            Inputs:
                setBoundary (bool): whether set boundary (for free-energy calculation)
                lowerBoundary (float): lower boundary of free-energy calculaton
                upperboundary (float): upper boundary of free-energy calculation
            Return:
                string: string of distance r definition '''
        
        string = f'\
colvar {{                            \n\
    name    r                        \n'
        if setBoundary:
            string += f'\
    width 0.1                        \n\
    lowerboundary {lowerBoundary:.1f}    \n\
    upperboundary {upperBoundary:.1f}    \n\
    subtractAppliedForce on          \n\
    expandboundaries  on             \n\
    extendedLagrangian on            \n\
    extendedFluctuation 0.1          \n'

        string += f'\
    distance {{                            \n\
        forceNoPBC       yes               \n\
        group1 {{                          \n\
            indexGroup  reference          \n\
	    }}                                 \n\
        group2 {{                          \n\
            indexGroup  ligand             \n\
        }}                                 \n\
    }}                                     \n\
}}                                         \n'

        return string

    def cvHeadTemplate(self, indexFile):
        ''' return the head of colvars file
            Inputs:
                indexFile (string): name of ndx file
            Return:
                string: head of colvars file '''
        return f'\
colvarsTrajFrequency      500              \n\
colvarsRestartFrequency   50000            \n\
indexFile                 {indexFile}      \n'

    def cvHarmonicWallsTemplate(self, cv, lowerWall, upperWall):
        ''' template of harmonic wall bias
            Inputs:
                cv (string): name of the colvars
                lowerWall (float): lower wall of the bias
                upperWall (float): upper wall of the bias
            Return:
                string: string of the harmonic wall bias definition '''
        string = f'\
harmonicWalls {{                           \n\
    colvars           {cv}                 \n\
    lowerWalls        {lowerWall:.1f}      \n\
    upperWalls        {upperWall:.1f}      \n\
    lowerWallConstant 0.2                  \n\
    upperWallConstant 0.2                  \n\
}}                                         \n'
        return string

    def cvHarmonicTemplate(self, cv, constant, center, tiWindows=0, tiForward=True, targetForceConstant = 0):
        ''' template for harmonic restraint bias
            Inputs:
                cv (string): name of the colvars
                constant (float): force constant of the restraint
                center (float): center of the restraint
                tiWindows (int): number of windows of the TI simulation (if runs a TI simulation)
                tiForward (bool): whether the TI simulation is forward (if runs a TI simulation)
                targetForceConstant (int): target force constant of the restraint in TI simulation (if runs a TI simulation)
            Return:
                string: string of the harmonic restraint definition '''
        string = f'\
harmonic {{                          \n\
    colvars         {cv}             \n\
    forceConstant   {constant:.1f}   \n\
    centers         {center:.1f}     \n'
        
        if tiWindows != 0:
            string += f'\
    targetNumSteps      500000                       \n\
    targetEquilSteps    100000                       \n\
    targetForceConstant {targetForceConstant}        \n\
    targetForceExponent 4                            \n'

            schedule = ''
            if tiForward:
                schedule += ' '.join([str(float(i) / float(tiWindows)) for i in range(tiWindows+1)])
            else:
                schedule += ' '.join([str(float(i) / float(tiWindows)) for i in range(tiWindows, -1, -1)])
            
            string += f'    lambdaSchedule {schedule}\n'

        string += '}\n'
        return string

    def cvABFTemplate(self, cv):
        ''' template for WTM-eABF bias
            Inputs:
                cv (string): name of the colvars
            Return:
                string: string of the WTM-eABF definition '''
        string = f'\
abf {{                            \n\
    colvars        {cv}           \n\
    FullSamples    10000          \n\
    historyfreq    50000          \n\
    writeCZARwindowFile           \n\
}}                                \n\
metadynamics {{                   \n\
    colvars           {cv}        \n\
    hillWidth         3.0         \n\
    hillWeight        0.05        \n\
    wellTempered      on          \n\
    biasTemperature   4000        \n\
}}                                \n'
        return string

    def cvProteinTemplate(self, centerCoor, refFile):
        ''' the template of restraining the protein
            Inputs:
                centerCoor (np.array, 3): (x,y,z), center of the protein 
                refFile (string): path of the reference file
            Return:
                string: string of the restraining the protein '''
        string = f'\
colvar {{                         \n\
  name translation                \n\
  distance {{                     \n\
    group1 {{                     \n\
      indexGroup  protein         \n\
    }}                            \n\
    group2 {{                     \n\
      dummyAtom ({centerCoor[0]}, {centerCoor[1]}, {centerCoor[2]})    \n\
    }}                            \n\
  }}                              \n\
}}                                \n\
harmonic {{                       \n\
  colvars       translation       \n\
  centers       0.0               \n\
  forceConstant 10.0              \n\
}}                                \n\
                                  \n\
colvar {{                         \n\
  name orientation                \n\
  orientation {{                  \n\
    atoms {{                      \n\
      indexGroup  protein         \n\
    }}                            \n\
    refPositionsFile   {refFile}  \n\
  }}                              \n\
}}                                \n\
harmonic {{                       \n\
  colvars       orientation       \n\
  centers       (1.0, 0.0, 0.0, 0.0)    \n\
  forceConstant 200.0             \n\
}}                                \n'
        return string
