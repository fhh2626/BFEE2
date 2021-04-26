# post-treatment of BFEE

import math

import numpy as np

# Boltazann constant for NAMD unit convention
BOLTZMANN = 0.0019872041
BOLTZMANN_GMX = 0.0019872041 * 4.184

# standard concentration
CSTAR = 1661
CSTAR_GMX = 1.661

class postTreatment:
    """the post-treatment of BFEE outputs
    """

    def __init__(self, temperature, unit, jobType='geometric'):
        """do post treatment, internally, all the unit should be converted into
           the NAMD/Colvars unit

        Args:
            temperature (float): temperature of the simulation
            unit (str): unit convention used by MD engine, 'namd' or 'gromacs'
            jobType (str): 'geometric' or 'alchemical'. Actually this arg is not used yet. Default to geometric.
        """

        if unit == 'namd':
            self.BOLTZMANN = BOLTZMANN
            self.CSTAR = CSTAR
        elif unit == 'gromacs':
            self.BOLTZMANN = BOLTZMANN_GMX
            self.CSTAR = CSTAR_GMX

        self.unit = unit
        self.beta = 1 / (self.BOLTZMANN * temperature)
        self.temperature = float(temperature)

    def _readPMF(self, filePath):
        """read a 1D PMF file

        Args:
            filePath (str): the path of the PMF file

        Returns:
            np.array (float, 2*N): ((x0,x1,x2, ...), (y0, y1, y2, ...))
        """
        
        data = np.loadtxt(filePath)
        x = data[:,0]
        y = data[:,1]

        return np.array((x, y))

    def _geometricRestraintContribution(self, pmf, forceConstant, rmsd=False, unbound=False):
        """calculate the contribution of RMSD and angle restraints

        Args:
            pmf (np.array, float, 2*N): ((x0,x1,x2, ...), (y0, y1, y2, ...))
            forceConstant (float): the force constant of the restraint
            rmsd (bool): whether the contribution of RMSD is being calculated. Defaults to False.
            unbound (bool, optional): whether unbound-state contribution is being calculated. Defaults to False.

        Returns:
            float: contribution of the geometric restraint
        """

        width = pmf[0][1] - pmf[0][0]

        if rmsd:
            # for RMSD, the restraintCenter is zero
            restraintCenter = 0
        else:
            # the minimum of pmf
            restraintCenter = pmf[0][np.argmin(pmf[1])]

        # integration
        numerator = 0
        denominator = 0
        for x, y in zip(pmf[0], pmf[1]):
            numerator += math.exp(-self.beta * y)
            denominator += math.exp((-self.beta) * (y + 0.5 * forceConstant * ((x - restraintCenter)**2)))
        
        contribution = math.log(numerator / denominator) / self.beta
        
        if unbound:
            return contribution
        else:
            return -contribution

    def _geometricRestraintContributionBulk(
        self, theta, forceConstantTheta, forceConstantPhi, forceConstantPsi
    ):
        """contribution of rotational restraints in the unbounded state

        Args:
            theta (float): restraining center of the theta angle
            forceConstantTheta (float): restraining force constant for Theta
            forceConstantPhi (float): restraining force constant for Phi
            forceConstantPsi (float): restraining force constant for Psi

        Returns:
            float: contribution of the geometric restraint in the unbound state
        """

        # all the units in radian
        theta0 = math.radians(theta)
        # periodic CV then the u(phi) and u(psi) should be the same in all cases
        phi0 = math.radians(180)
        psi0 = math.radians(180)

        forceConstantTheta *= (180 / math.pi)**2
        forceConstantPhi *= (180 / math.pi)**2
        forceConstantPsi *= (180 / math.pi)**2

        contributionTheta = 0
        contributionPhi = 0
        contributionPsi = 0
        # integration
        for i in range(1000):
            theta = i / 1000.0 * math.pi - math.pi / 2
            contributionTheta += 1.0/1000.0 * math.pi * math.sin(theta + math.pi / 2) * \
                 math.exp(-self.beta * 0.5 * forceConstantTheta * ((theta - theta0)**2))

            phi = i / 1000.0 * 2 * math.pi
            contributionPhi += 1.0/1000.0 * 2 * math.pi * \
                 math.exp(-self.beta * 0.5 * forceConstantPhi* ((phi - phi0)**2))

            psi = i / 1000.0 * 2 * math.pi
            contributionPsi += 1.0/1000.0 * 2 * math.pi * \
                 math.exp(-self.beta * 0.5 * forceConstantPsi* ((psi - psi0)**2))

        return -math.log((contributionTheta * contributionPhi * contributionPsi) / 8 / math.pi**2) / self.beta

    def _geometricJacobianCorrection(self, pmf):
        """correct the Jacobian contribution of separation pmf

        Args:
            pmf (np.array, float, 2*N): ((x0,x1,x2, ...), (y0, y1, y2, ...)), separation pmf
        """
        
        for i in range(len(pmf[0])):
            pmf[1][i] += 2 * self.BOLTZMANN * self.temperature * math.log(pmf[0][i])

        pmf[1] -= np.min(pmf[1])

    def _geometricCalculateSI(
        self, rStar, pmf, polarTheta, polarPhi, forceConstantPolarTheta, forceConstantPolarPhi):
        """calculation the contribution of S* and I* in the separation simulation

        Args:
            rStar (float): r* in integration
            pmf (np.array, float, 2*N): ((x0,x1,x2, ...), (y0, y1, y2, ...)), separation pmf
            polarTheta0 (float): restraining center of polarTheta
            polarPhi0 (float): restraining center of polarPhi
            forceConstantPolarTheta (float): restraining force constant for polarTheta
            forceConstantPolarPhi (float): restraining force constant for polarPhi

        Returns:
            float: contribution of S* and I* in the separation simulation
        """

        assert(rStar <= pmf[0][-1])

        polarTheta0 = math.radians(polarTheta)
        polarPhi0 = math.radians(polarPhi)

        forceConstantPolarTheta *= (180 / math.pi)**2
        forceConstantPolarPhi *= (180 / math.pi)**2

        contributionPolarTheta = 0
        contributionPolarPhi = 0
        # integration
        for i in range(1000):
            polarTheta = i / 1000.0 * math.pi
            contributionPolarTheta += 1.0 / 1000.0 * math.pi * math.sin(polarTheta) * \
                math.exp(-self.beta * 0.5 * forceConstantPolarTheta * (polarTheta - polarTheta0)**2)

            polarPhi = i / 1000.0 * 2 * math.pi - math.pi
            contributionPolarPhi += 1.0 / 1000.0 * 2 * math.pi * \
                math.exp(-self.beta * 0.5 * forceConstantPolarPhi * (polarPhi - polarPhi0)**2)

        S = rStar**2 * contributionPolarTheta * contributionPolarPhi

        # w(r*)
        wrStar = pmf[1][0]
        for x, y in zip(pmf[0], pmf[1]):
            if x >= rStar:
                wrStar = y
                break

        # integration
        width = pmf[0][1] - pmf[0][0]
        I = 0
        for x, y in zip(pmf[0], pmf[1]):
            I += width * math.exp(-self.beta * (y - wrStar))
            if x >= rStar:
                break
        
        return -1 / self.beta * math.log(S * I / self.CSTAR)

    def geometricBindingFreeEnergy(self, filePathes, parameters):
        """calculate binding free energy for geometric route

        Args:
            filePathes (list of strings, 8): pathes of PMF files for step1 - step8
            parameters (np.array, floats, 8): (forceConstant1, FC2, FC3, FC4, FC5, FC6, r*, FC8)

        Returns:
            np.array, float, 10: (contributions for step1, 2, 3, 4 ... 8, bulk restraining contribution, free energy)
        """

        assert len(parameters) == 8
        assert len(filePathes) == 8

        pmfs = []
        for path in filePathes:
            pmfs.append(self._readPMF(path))
        self._geometricJacobianCorrection(pmfs[6])

        contributions = np.zeros(10)
        contributions[0] = self._geometricRestraintContribution(pmfs[0], parameters[0], True, False)
        contributions[1] = self._geometricRestraintContribution(pmfs[1], parameters[1], False, False)
        contributions[2] = self._geometricRestraintContribution(pmfs[2], parameters[2], False, False)
        contributions[3] = self._geometricRestraintContribution(pmfs[3], parameters[3], False, False)
        contributions[4] = self._geometricRestraintContribution(pmfs[4], parameters[4], False, False)
        contributions[5] = self._geometricRestraintContribution(pmfs[5], parameters[5], False, False)
        contributions[6] = self._geometricCalculateSI(
            parameters[6], pmfs[6], pmfs[4][0][np.argmin(pmfs[4][1])], pmfs[5][0][np.argmin(pmfs[5][1])],
            parameters[4], parameters[5]
        )
        contributions[7] = self._geometricRestraintContribution(pmfs[7], parameters[7], True, True)
        contributions[8] = self._geometricRestraintContributionBulk(
            pmfs[1][0][np.argmin(pmfs[1][1])], parameters[1], parameters[2], parameters[3]
        )

        contributions[9] = np.sum(contributions[:9])

        if self.unit == 'namd':
            return contributions
        elif self.unit == 'gromacs':
            return contributions / 4.184

    def _alchemicalRestraintContributionBulk(
        self, eulerTheta, polarTheta, R, 
        forceConstantTheta=0.1, forceConstantPhi=0.1, forceConstantPsi=0.1,
        forceConstanttheta=0.1, forceConstantphi=0.1, forceConstantR=10
    ):
        """contribution of (standard concentration corrected) rotational 
           and orienetational restraints in the unbounded state

        Args:
            eulerTheta (float): restraining center of the Euler angle theta
            polarTheta (float): restraining center of the polar angle theta
            R (float): restraining center of anger R
            forceConstantTheta (float): restraining force constant for euler Theta. Defaults to 0.1.
            forceConstantPhi (float, optional): restraining force constant for euler Phi. Defaults to 0.1.
            forceConstantPsi (float, optional): restraining force constant for euler Psi. Defaults to 0.1.
            forceConstanttheta (float, optional): restraining force constant for polar theta. Defaults to 0.1.
            forceConstantphi (float, optional): restraining force constant for polar phi. Defaults to 0.1.
            forceConstantR (int, optional): restraining force constant for distance R. Defaults to 10.

        Returns:
            float: contribution of the geometric restraint in the unbound state
        """

        # degrees to rad
        eulerTheta = math.radians(eulerTheta + 90)
        polarTheta = math.radians(polarTheta)
        forceConstantTheta *= (180 / math.pi)**2
        forceConstantPhi *= (180 / math.pi)**2
        forceConstantPsi *= (180 / math.pi)**2
        forceConstanttheta *= (180 / math.pi)**2
        forceConstantphi *= (180 / math.pi)**2

        contribution = self.BOLTZMANN * self.temperature * math.log(
            8 * (math.pi**2) * self.CSTAR / ((R**2) * math.sin(eulerTheta) * math.sin(polarTheta)) * \
            math.sqrt(forceConstantTheta * forceConstantPhi * forceConstantPsi * forceConstanttheta * \
            forceConstantphi * forceConstantR ) / ((2 * math.pi * self.BOLTZMANN * self.temperature)**3)
        )
        return contribution

    def _alchemicalFepoutFile(self, filePath, fileType = 'fepout'):
        """parse a fepout/log file and return the total free energy change

        Args:
            filePath (str): path of the fepout file
            fileType (str): 'fepout' (decouping atoms) or 'log' (decoupling restraints). Defaults to 'fepout'.

        Returns:
            float: free-energy change
        """
        
        freeEnergy = 0

        # for TI simulation results of the .log file
        Lambda = []
        dA_dLambda = []

        with open(filePath, 'r', encoding='utf-8') as fepoutFile:
            for line in fepoutFile.readlines():
                if fileType == 'fepout':
                    if not line.startswith('#Free energy'):
                        continue
                if fileType == 'log':
                    if not ('dA/dLambda' in line):
                        continue
                splitedLine = line.strip().split()
                if fileType == 'fepout':
                    freeEnergy += float(splitedLine[11])
                if fileType == 'log':
                    Lambda.append(float(splitedLine[4]))
                    dA_dLambda.append(float(splitedLine[6]))

        # integraion of TI results
        if fileType == 'log':
            # backward simulation
            if Lambda[0] > Lambda[7]:
                Lambda.reverse()
                dA_dLambda.reverse()

            # step 2
            if Lambda[0] == Lambda[1]:
                for i in range(7, len(Lambda)):
                    freeEnergy += dA_dLambda[i] * (Lambda[i] - Lambda[i-7])
            else:
                # step 4
                for i in range(1, len(Lambda)):
                    freeEnergy += dA_dLambda[i] * (Lambda[i] - Lambda[i-1])

        return freeEnergy

    def alchemicalBindingFreeEnergy(self, filePathes, parameters):
        """calculate binding free energy for geometric route

        Args:
            filePathes (list of strings, 8): pathes of alchemical output files
                                             (step1-forward, step1-backward, step2-forward ...)
            parameters (np.array, floats, 9): (eulerTheta, polarTheta, r, forceConstant1, FC2, FC3, FC4, FC5, FC6)

        Returns:
            tuple:
                np.array, float, 6: (contributions for step1, 2, 3, 4, bulk restraining contribution, free energy)
                np.array, float, 6: errors corresponding each contribution
        """

        assert len(parameters) == 9
        assert len(filePathes) == 8

        # get free energies from fep outputs
        freeEnergies = []
        for i in range(len(filePathes)):
            if filePathes[i] != '':
                if (i // 2) % 2 == 0:
                    freeEnergies.append(self._alchemicalFepoutFile(filePathes[i], 'fepout'))
                else:
                    freeEnergies.append(self._alchemicalFepoutFile(filePathes[i], 'log'))
            else:
                # backward file can be empty
                freeEnergies.append(None)

        contributions = np.zeros(6)
        errors = np.zeros(6)

        if freeEnergies[1] is not None:
            contributions[0] = (freeEnergies[0] - freeEnergies[1]) / 2
            errors[0] = abs((freeEnergies[0] + freeEnergies[1]) / 1.414)
        else:
            contributions[0] = freeEnergies[0]
            errors[0] = 99999

        if freeEnergies[3] is not None:
            contributions[1] = -(freeEnergies[2] + freeEnergies[3]) / 2
            errors[1] = abs((freeEnergies[2] - freeEnergies[3]) / 1.414)
        else:
            contributions[1] = -freeEnergies[2]
            errors[1] = 99999

        if freeEnergies[5] is not None:
            contributions[2] = -(freeEnergies[4] - freeEnergies[5]) / 2
            errors[2] = abs((freeEnergies[4] + freeEnergies[5]) / 1.414)
        else:
            contributions[2] = -freeEnergies[4]
            errors[2] = 99999

        if freeEnergies[7] is not None:
            contributions[3] = (freeEnergies[6] + freeEnergies[7]) / 2
            errors[3] = abs((freeEnergies[6] - freeEnergies[7]) / 1.414)
        else:
            contributions[3] = freeEnergies[6]
            errors[3] = 99999

        contributions[4] = self._alchemicalRestraintContributionBulk(*parameters)
        errors[4] = 0

        contributions[5] = contributions[0] + contributions[1] + contributions[2] + contributions[3] + contributions[4]
        errors[5] = math.sqrt(errors[0]**2 + errors[1]**2 +errors[2]**2 + errors[3]**2 + errors[4]**2)

        return contributions, errors
