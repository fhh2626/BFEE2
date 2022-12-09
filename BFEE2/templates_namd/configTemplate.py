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
                            membraneProtein = False,
                            OPLSMixingRule = False,
                            GaWTM = False
                            ):
        """the namd config file template

        Args:
            forceFieldType (str): 'charmm' or 'amber'
            forceFieldFiles (list of str): name of charmm force field files
            topFile (str): name of the topology file (psf, parm) 
            coorFile (str): name of the coordinate file (pdb, rst)
            NAMDRestartCoor (str): name of namd binary restart coor file (if restart from a previous simulation)
            NAMDRestartVel (str): name of namd binary restart vel file
            NAMDRestartXsc (str): name of namd restart xsc file
            PBCCondition (np.array, 2*3): PBC vector, ((lengthX, lengthY, lengthZ),(centerX, centerY, centerZ))
            outputPrefix (str): prefix of output file
            temperature (float): temperature of the simulation
            numSteps (int): number of steps of the simulation
            cvFile (str): name of Colvars file. Defaults to ''.
            cvDefinitionFile (str, optional): name of a TCL file defining new CVs. Defaults to ''.
            CVRestartFile (str, optional): name of Colvars restart file. Defaults to ''.
            fepFile (str, optional): name of fep file, indicating which atoms will be generated/removed 
                                     (if run alchemical simulation). Defaults to ''.
            fepWindowNum (int, optional): number of fep windows. Defaults to 20.
            fepForward (bool, optional): whether this is a forward fep simulation. Defaults to True.
            fepDoubleWide (bool, optional): whether this is a double-wide fep simulation. Defaults to False.
            fepMinBeforeSample (bool, optional): whether do minimization before sampling in each FEP window.
                                                 Defaults to False.
            membraneProtein (bool, optional): whether simulating a membrame protein. Defaults to False.
            OPLSMixingRule (bool, optional): whether use the OPLS mixing rules. Defaults to False.
            GaWTM (bool, optional): Whether this is an GaWTM-eABF simulation. Default to False
            
        Returns:
            str: a NAMD config string if succeed, and empty string otherwise
        """

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
        
        # OPLS mixing rule
        if OPLSMixingRule:
            configString += f'\
vdwGeometricSigma    yes               \n'

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
langevinTemp         {temperature}              \n\
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
colvars    on                                   \n'
            if not GaWTM:
                configString += f'\
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
            if not GaWTM:
                configString += f'\
run    {numSteps}                               \n'
            else:
                configString += f'\
accelMD                         on                     \n\
accelMDG                        on                     \n\
accelMDGcMDSteps                500000                 \n\
accelMDGcMDPrepSteps            100000                 \n\
accelMDGEquiPrepSteps           100000                 \n\
accelMDGEquiSteps               500000                 \n\
accelMDdihe                     on                     \n\
accelMDOutFreq                  1000                   \n\
accelMDGSigma0D                 6.0                    \n\
for {{set stage 0}} {{$stage < 2}} {{incr stage}} {{   \n\
    if {{$stage == 0}} {{                              \n\
        puts "Probing the GaMD parameters..."          \n\
        cv configfile {cvFile}                         \n\
        run norepeat 1000000                           \n\
    }} elseif {{$stage == 1}} {{                       \n\
        puts "Starting eABF + GaMD..."                 \n\
        cv reset                                       \n\
        cv configfile {cvFile + ".amd"}                \n\
        run norepeat   {numSteps}                      \n\
    }}                                                 \n\
}}                                                     \n'
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
    
    def gromacsMinimizeConfigTemplate(self):
        """the gromacs config file template for minimization
        
        Returns:
            str: a Gromacs config string if succeed, and empty string otherwise
        """
        configString = f'\
integrator      = steep     \n\
emtol           = 418.4     \n\
emstep          = 0.01      \n\
nsteps          = 20000     \n\
nstlist         = 1         \n\
cutoff-scheme   = Verlet    \n\
ns_type         =           \n\
coulombtype     = PME       \n\
rcoulomb        = 1.2       \n\
rvdw            = 1.2       \n\
pbc             = xyz       \n'
        return configString
    
    def gromacsConfigTemplate(
                            self,
                            forceFieldType,
                            temperature,
                            numSteps,
                            membraneProtein = False,
                            generateVelocities = True
                            ):
        """the gromacs config file template

        Args:
            forceFieldType (str): 'charmm' or 'amber'
            temperature (float): temperature of the simulation
            numSteps (int): number of steps of the simulation
            membraneProtein (bool, optional): whether simulating a membrame protein. Defaults to False.
            generateVelocities (bool, optional): whether generate velocities. Defaults to True.
            
        Returns:
            str: a Gromacs config string if succeed, and empty string otherwise
        """
        
        configString = f'\
integrator              = md                \n\
nsteps                  = {numSteps}        \n\
dt                      = 0.002             \n\
nstxout                 = 5000              \n\
nstvout                 = 5000              \n\
nstenergy               = 5000              \n\
nstlog                  = 5000              \n\
continuation            = no                \n\
constraint_algorithm    = lincs             \n\
constraints             = h-bonds           \n\
cutoff-scheme           = Verlet            \n\
nstlist                 = 10                \n\
coulombtype             = PME               \n'

        if forceFieldType == 'charmm':
            configString += f'\
rcoulomb                = 1.2               \n\
rvdw                    = 1.2               \n\
rvdw-switch             = 1.0               \n'

        elif forceFieldType == 'amber':
            configString += f'\
rcoulomb                = 0.9               \n\
rvdw                    = 0.9               \n\
rvdw-switch             = 0.8               \n'

        configString += f'\
vdwtype                 = Cut-off                \n\
vdw-modifier            = force-switch           \n\
DispCorr                = EnerPres               \n\
tcoupl                  = V-rescale              \n\
tc-grps                 = System                 \n\
tau_t                   = 1                      \n\
ref_t                   = {temperature}          \n\
pbc                     = xyz                    \n\
pcoupl                  = C-rescale              \n\
tau-p                   = 1                      \n'

        if not membraneProtein:
            configString += f'\
ref-p                   = 1.01325                \n\
compressibility         = 4.5e-5                 \n\
refcoord-scaling        = com                    \n\
pcoupltype              = isotropic              \n'

        elif membraneProtein:
            configString += f'\
ref-p                   = 1.01325 1.01325        \n\
compressibility         = 4.5e-5 4.5e-5          \n\
refcoord-scaling        = com                    \n\
pcoupltype              = semiisotropic          \n'
    
        if generateVelocities:
            configString += f'\
gen_vel                 = yes                    \n\
gen_temp                = {temperature}          \n\
gen_seed                = -1                     \n'

        return configString


    def cvRMSDTemplate(
            self, setBoundary, lowerBoundary, upperBoundary, refFile, 
            extendedLagrangian = True, reflectingBoundary = False,
            unit = 'namd'
        ):
        """RMSD CV template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            refFile (str): path to the reference file
            extendedLagrangian (bool, optional): Whether extended Lagrangian is added. Default to True
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
            unit (str, optional): unit, 'namd' or 'gromacs'. Default to namd.
        
        Returns:
            str: string of RMSD definition
        """
        
        if unit == 'namd':
            scaleFactor = 1
        elif unit == 'gromacs':
            scaleFactor = 0.1

        string = f'\
colvar {{                                    \n\
    name RMSD                                \n'

        if setBoundary:
            string += f'\
    width {0.05 * scaleFactor:.2f}                               \n\
    lowerboundary {lowerBoundary * scaleFactor:.2f}            \n\
    upperboundary {upperBoundary * scaleFactor:.2f}            \n'
    
        if setBoundary and reflectingBoundary:
            string += f'\
    reflectingLowerboundary  on              \n\
    reflectingUpperboundary  on              \n'
    
        if extendedLagrangian:
            string += f'\
    subtractAppliedForce on                  \n\
    expandboundaries  on                     \n\
    extendedLagrangian on                    \n\
    extendedFluctuation {0.05 * scaleFactor:.2f}                 \n'

        string += f'\
    rmsd {{                                  \n\
        atoms {{                             \n\
            indexGroup  ligand               \n\
        }}                                   \n\
        refpositionsfile  {refFile}          \n\
    }}                                       \n\
}}                                           \n'
        return string
    
    def cvAngleTemplate(
            self, setBoundary, lowerBoundary, upperBoundary, angle, 
            refFile, oldDefinition = True, extendedLagrangian = True, reflectingBoundary = False
        ):
        """Eulaer and polar angle template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperBoundary (float): upper boundary of free-energy calculation
            angle (str): 'eulerTheta', 'eulerPhi', 'eulerPsi', 'polarTheta' or 'polarPhi'
            refFile (str): path to the reference file
            oldDefinition (bool, optional): Whether use old definition of angles
                                            for compatibility. Defaults to True
            extendedLagrangian (bool, optional): Whether extended Lagrangian is added. Default to True
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
        """
        
        assert(
            angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi' or \
            angle == 'polarTheta' or angle == 'polarPhi'
        )
        
        if angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi':
            if oldDefinition:
                return self.cvEulerAngleTemplate(
                           setBoundary, lowerBoundary, upperBoundary, angle, refFile, extendedLagrangian, reflectingBoundary
                       )
            else:
                return self.newCvEulerAngleTemplate(
                           setBoundary, lowerBoundary, upperBoundary, angle, refFile, extendedLagrangian, reflectingBoundary
                       )
        elif angle == 'polarTheta' or angle == 'polarPhi':
            if oldDefinition:
                return self.cvPolarAngleTemplate(
                           setBoundary, lowerBoundary, upperBoundary, angle, refFile, extendedLagrangian, reflectingBoundary
                       )
            else:
                return self.newCvPolarAngleTemplate(
                           setBoundary, lowerBoundary, upperBoundary, angle, refFile, extendedLagrangian, reflectingBoundary
                       )

    def cvEulerAngleTemplate(
            self, setBoundary, lowerBoundary, upperBoundary, angle, 
            refFile, extendedLagrangian = True, reflectingBoundary = False
        ):
        """Euler angle template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'eulerTheta', 'eulerPhi' or 'eulerPsi'
            refFile (str): path to the reference file
            extendedLagrangian (bool, optional): Whether extended Lagrangian is added. Default to True
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
            
        Returns:
            string: string of Euler angle definition
        """

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
    upperboundary {upperBoundary:.1f}      \n'
    
        if setBoundary and reflectingBoundary:
            string += f'\
    reflectingLowerboundary  on              \n\
    reflectingUpperboundary  on              \n'
        
        if extendedLagrangian:
            string += f'\
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

    def cvPolarAngleTemplate(
            self, setBoundary, lowerBoundary, upperBoundary, angle, 
            refFile, extendedLagrangian = True, reflectingBoundary = False
        ):
        """Polar angle template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'polarTheta' or 'polarPhi'
            refFile (str): path to the reference file
            extendedLagrangian (bool, optional): Whether extended Lagrangian is added. Default to True
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
            
        Return:
            str: string of polar angle definition
        """

        assert(angle == 'polarTheta' or angle == 'polarPhi')
        
        string = f'\
colvar {{                                   \n\
    name {angle}                            \n'

        if angle == 'polarTheta':
            string += f'\
    customFunction acos(-i2) * 180 / 3.1415926\n'
        elif angle == 'polarPhi':
            string += f'\
    customFunction atan2(i3, i1) * 180 / 3.1415926\n\
    period  360                             \n\
    wrapAround 0.0                          \n'

        if setBoundary:
            string += f'\
    width 1                                 \n\
    lowerboundary {lowerBoundary:.1f}           \n\
    upperboundary {upperBoundary:.1f}           \n'
    
        if setBoundary and reflectingBoundary:
            string += f'\
    reflectingLowerboundary  on              \n\
    reflectingUpperboundary  on              \n'
    
        if extendedLagrangian:
            string += f'\
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
    
    def newCvEulerAngleTemplate(
            self, setBoundary, lowerBoundary, upperBoundary, angle, 
            refFile, extendedLagrangian = True, reflectingBoundary = False
        ):
        """new definition Euler angle template, probably the pinning down the protein is not required

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'eulerTheta', 'eulerPhi' of 'eulerPsi'
            refFile (str): path to the reference file
            extendedLagrangian (bool, optional): Whether extended Lagrangian is added. Default to True
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
            
        Returns:
            string: string of Euler angle definition
        """

        assert(angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi')
        
        string = f'\
colvar {{                              \n\
    name {angle}                   \n'

        if setBoundary:
            string += f'\
    width 1                            \n\
    lowerboundary {lowerBoundary:.1f}      \n\
    upperboundary {upperBoundary:.1f}      \n'
    
        if setBoundary and reflectingBoundary:
            string += f'\
    reflectingLowerboundary  on              \n\
    reflectingUpperboundary  on              \n'
        
        if extendedLagrangian:
            string += f'\
    subtractAppliedForce on            \n\
    expandboundaries  on               \n\
    extendedLagrangian on              \n\
    extendedFluctuation 1              \n'

        string += f'\
    {angle} {{                             \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
            centerToOrigin     on              \n\
	        enableFitGradients on              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
         refpositionsfile  {refFile}           \n\
    }}                                         \n\
}}                                             \n'

        return string
    
    def newCvPolarAngleTemplate(
            self, setBoundary, lowerBoundary, upperBoundary, angle, 
            refFile, extendedLagrangian = True, reflectingBoundary = False
        ):
        """new definition of Polar angle template, probably the pinning down the protein is not required

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'polarTheta' or 'polarPhi'
            refFile (str): path to the reference file
            extendedLagrangian (bool, optional): Whether extended Lagrangian is added. Default to True
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
            
        Return:
            str: string of polar angle definition
        """

        assert(angle == 'polarTheta' or angle == 'polarPhi')
        
        string = f'\
colvar {{                                   \n\
    name {angle}                            \n'
    
        if angle == 'polarTheta':
            string += f'\
    customFunction acos(-sin(t / 180 * 3.1415926) * sin(p / 180 * 3.1415926)) * 180 / 3.1415926\n'
        elif angle == 'polarPhi':
            string += f'\
    customFunction atan2(cos(t / 180 * 3.1415926), cos(p / 180 * 3.1415926) * sin(t / 180 * 3.1415926)) * 180 / 3.1415926\n\
    period  360                             \n\
    wrapAround 0.0                          \n'

        if setBoundary:
            string += f'\
    width 1                                 \n\
    lowerboundary {lowerBoundary:.1f}           \n\
    upperboundary {upperBoundary:.1f}           \n'
    
        if setBoundary and reflectingBoundary:
            string += f'\
    reflectingLowerboundary  on              \n\
    reflectingUpperboundary  on              \n'
        
        if extendedLagrangian:
            string += f'\
    subtractAppliedForce on                 \n\
    expandboundaries  on                    \n\
    extendedLagrangian on                   \n\
    extendedFluctuation 1                   \n'

        string += f'\
    polarTheta {{                             \n\
        name        t                          \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
            centerToOrigin     on              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
    }}                                         \n\
    polarPhi {{                             \n\
        name        p                      \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
            centerToOrigin     on              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
    }}                                         \n\
}}                                             \n'
        return string

    def cvRTemplate(
            self, setBoundary, lowerBoundary, upperBoundary, 
            extendedLagrangian = True, reflectingBoundary = False,
            unit = 'namd'
        ):
        """r distance template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy
            extendedLagrangian (bool, optional): Whether extended Lagrangian is added. Default to True
            reflectingBoundary (bool, optional): Whether use reflecting boundaries, requires setBoundary on. Default to False
            unit (str, optional): unit, 'namd' or 'gromacs'. Default to namd.
        
        Returns:
            str: string of distance r definition
        """
        
        if unit == 'namd':
            scaleFactor = 1
        elif unit == 'gromacs':
            scaleFactor = 0.1
        
        string = f'\
colvar {{                            \n\
    name    r                        \n'
        if setBoundary:
            string += f'\
    width {0.1 * scaleFactor:.2f}                        \n\
    lowerboundary {lowerBoundary * scaleFactor:.2f}    \n\
    upperboundary {upperBoundary * scaleFactor:.2f}    \n'
    
        if setBoundary and reflectingBoundary:
            string += f'\
    reflectingLowerboundary  on              \n\
    reflectingUpperboundary  on              \n'
        
        if extendedLagrangian:
            string += f'\
    subtractAppliedForce on          \n\
    expandboundaries  on             \n\
    extendedLagrangian on            \n\
    extendedFluctuation {0.1 * scaleFactor:.2f}          \n'

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

    def cvHeadTemplate(self, indexFile, reweightAMD=False):
        """return the head of colvars file

        Args:
            indexFile (str): name of ndx file
            reweightAMD (bool, optional): if reweightAMD, there should be "smp off" in the 
                                          head. Default to False.

        Returns:
            str: head of colvars file
        """
        
        string = f'\
colvarsTrajFrequency      5000             \n\
colvarsRestartFrequency   5000            \n\
indexFile                 {indexFile}      \n'
        if reweightAMD:
            string += f'\
smp                       off              \n'
        return string

    def cvHarmonicWallsTemplate(self, cv, lowerWall, upperWall, unit = 'namd'):
        ''' template of harmonic wall bias
        
        Args:
            cv (str): name of the colvars
            lowerWall (float): lower wall of the bias
            upperWall (float): upper wall of the bias
            unit (str, optional): unit, 'namd' or 'gromacs'. Default to namd.
                
        Returns:
            str: string of the harmonic wall bias definition '''
        
        if cv == 'RMSD' or cv == 'r':
            distanceCV = True
        else:
            distanceCV = False
            
        if unit == 'gromacs' and distanceCV:
            scaleFactor = 0.1
        else:
            scaleFactor = 1
            
        string = f'\
harmonicWalls {{                           \n\
    colvars           {cv}                 \n\
    lowerWalls        {lowerWall * scaleFactor:.1f}      \n\
    upperWalls        {upperWall * scaleFactor:.1f}      \n\
    lowerWallConstant 0.2                  \n\
    upperWallConstant 0.2                  \n\
}}                                         \n'
        return string

    def cvHarmonicTemplate(
            self, cv, constant, center, tiWindows=0, tiForward=True, targetForceConstant = 0,
            unit = 'namd'
        ):
        """template for a harmonic restraint

        Args:
            cv (str): name of the colvars
            constant (float): force constant of the restraint
            center (float): center of the restraint
            tiWindows (int): number of windows of the TI simulation (if runs a TI simulation). Defaults to 0.
            tiForward (bool, optional): whether the TI simulation is forward (if runs a TI simulation). Defaults to True.
            targetForceConstant (int, optional): targeted force constant of the restraint in TI simulation (if runs a TI simulation).
                                                 Defaults to 0.
            unit (str, optional): unit, 'namd' or 'gromacs'. Default to namd.
        
        Returns:
            str: string of the harmonic restraint definition
        """
        
        if cv == 'RMSD' or cv == 'r':
            distanceCV = True
        else:
            distanceCV = False
        
        if unit == 'gromacs' and distanceCV:
            forceConstantScaleFactor = 418.4
            centerScaleFactor = 0.1
        elif unit == 'gromacs' and not distanceCV:
            forceConstantScaleFactor = 4.164
            centerScaleFactor = 1
        else:
            forceConstantScaleFactor = 1
            centerScaleFactor = 1
            

        string = f'\
harmonic {{                          \n\
    colvars         {cv}             \n\
    forceConstant   {constant * forceConstantScaleFactor:.1f}   \n\
    centers         {center * centerScaleFactor:.1f}     \n'
        
        if tiWindows != 0:
            string += f'\
    targetNumSteps      500000                       \n\
    targetEquilSteps    100000                       \n\
    targetForceConstant {targetForceConstant * forceConstantScaleFactor:.1f}        \n\
    targetForceExponent 4                            \n'

            schedule = ''
            if tiForward:
                schedule += ' '.join([str(float(i) / float(tiWindows)) for i in range(tiWindows+1)])
            else:
                schedule += ' '.join([str(float(i) / float(tiWindows)) for i in range(tiWindows, -1, -1)])
            
            string += f'    lambdaSchedule {schedule}\n'

        string += '}\n'
        return string

    def cvABFTemplate(self, cv, unit = 'namd'):
        ''' template for WTM-eABF bias
        
        Args:
            cv (str): name of the colvars
            unit (str, optional): unit, 'namd' or 'gromacs'. Default to namd.
            
        Returns:
            str: string of the WTM-eABF definition '''
            
        if unit == 'namd':
            scaleFactor = 1
        elif unit == 'gromacs':
            scaleFactor = 4.184
            
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
    hillWeight        {0.05 * scaleFactor:.2f}        \n\
    wellTempered      on          \n\
    biasTemperature   4000        \n\
}}                                \n'
        return string
    
    def cvHistogramTemplate(self, cv):
        """ template for outputting the histogram of the CV

        Args:
            cv (str): name of the colvars
            
        Returns:
            str: string of the histogram definition
        """
        
        string = f'\
histogram {{                     \n\
  colvars    {cv}                \n\
  outputFileDX  none             \n\
  outputFreq 10000               \n\
}}                               \n'
        return string

    def cvReweightAMDTemplate(self, cv):
        """ template for reweight aMD results

        Args:
            cv (str): name of the colvars
            
        Returns:
            str: string of the reweightamd definition
        """
        
        string = f'\
reweightamd {{                  \n\
  colvars  {cv}                 \n\
}}                              \n'
        return string

    def cvProteinTemplate(self, centerCoor, refFile, unit = 'namd'):
        """the template of restraining the protein

        Args:
            centerCoor (np.array, 3): (x,y,z), center of the protein 
            refFile (str): path of the reference file
            unit (str, optional): unit, 'namd' or 'gromacs'. Default to namd.
        
        Returns:
            str: string of the restraining the protein
        """
        
        if unit == 'namd':
            translationCoorScaleFactor = 1
            translationForceScaleFactor = 1
            orientationForceScaleFactor = 1
        elif unit == 'gromacs':
            translationCoorScaleFactor = 0.1
            translationForceScaleFactor = 418.4
            orientationForceScaleFactor = 4.184
        
        string = f'\
colvar {{                         \n\
  name translation                \n\
  distance {{                     \n\
    group1 {{                     \n\
      indexGroup  protein         \n\
    }}                            \n\
    group2 {{                     \n\
      dummyAtom ({centerCoor[0] * translationCoorScaleFactor}, {centerCoor[1] * translationCoorScaleFactor}, {centerCoor[2] * translationCoorScaleFactor})    \n\
    }}                            \n\
  }}                              \n\
}}                                \n\
harmonic {{                       \n\
  colvars       translation       \n\
  centers       0.0               \n\
  forceConstant {100.0 * translationForceScaleFactor:.1f}             \n\
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
  forceConstant {2000.0 * orientationForceScaleFactor:.1f}            \n\
}}                                \n'
        return string
