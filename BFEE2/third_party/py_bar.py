import sys, math
from typing import List, Literal, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

# VERSION
VERSION = 1.0

# boltzmann constant
BOLTZMANN = 0.0019872041
CSTAR = 1661

# type of FEP files
fep_type = Literal['forward', 'backward', 'double-wide']

class NAMDParser:
    """ parse NAMD .fepout files and get the necessary data
    """
    
    def __init__(self, forward_file: str, backward_file: Optional[str] = '') -> None:
        """ Read NAMD .fepout files. The end-user can either provide the bidirectional 
            .fepout files, or provide an .fepout file of the double-wide FEP simulation

        Args:
            forward_file (str): the fepout file for forward or double-wide simulation 
            backward_file (Optional[str], optional): the fepout file for backward simulation.
                                                     Defaults to ''.
        """
        
        self._double_wide = False
        if (backward_file == '') or (backward_file is None):
            self._double_wide = True
        
        # List[Tuple[float, float]], recording the boundaries of each window
        # List[Tuple[NDArray, NDArray]], recording the deltaU of forward of 
        # backward simulations of each window
        if not self._double_wide:
            self._windows, self._deltaU_data = self._pair_bidirectionalData(
                                                   self._read_fepout(forward_file),
                                                   self._read_fepout(backward_file)
                                               )
        else:
            self._windows, self._deltaU_data = self._read_double_wide_fepout(forward_file)
            
        if len(self._windows) != len(self._deltaU_data):
            raise RuntimeError('Internal numbers of windows and deltaU do not match! This is a bug!')
        
    def _read_fepout(self, fepout_file: str) -> Tuple[List[Tuple[float, float]], List[NDArray]]:
        """ Read an NAMD fepout file. Return the window and deltaU information

        Args:
            fepout_file (str): the path of the fepout file

        Returns:
            Tuple[List[Tuple[float, float]], List[NDArray]]: List[Tuple[float, float]], recording 
                                                             the boundaries of each window, and
                                                             List[NDArray], recording the deltaU
        """
        
        windows = []
        deltaU = []
        with open(fepout_file, 'r') as input_fepout:
            while True:
                line = input_fepout.readline()
                if not line:
                    break
                
                if line.startswith('#NEW FEP WINDOW:'):
                    splitedLine = line.strip().split()
                    windows.append((float(splitedLine[6]), float(splitedLine[8])))
                    continue
                
                if line.startswith('#STARTING COLLECTION'):
                    # collecting deltaU
                    deltaU_per_window = []
                    while True:
                        line = input_fepout.readline()
                        if line.startswith('FepEnergy:'):
                            splitedLine = line.strip().split()
                            deltaU_per_window.append(float(splitedLine[6]))
                        else:
                            deltaU.append(np.array(deltaU_per_window))
                            break
                        
        return windows, deltaU
    
    def _read_double_wide_fepout(self, fepout_file: str) -> Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]:
        """ Read an NAMD double-wide fepout file. Return the window and bidirectional deltaU information 

        Args:
            fepout_file (str): the path of the fepout file

        Returns:
            Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]: recording the boundary and 
                                                                             deltaU of each window of
                                                                             bidirectional simulations
        """
        
        windows = []
        deltaU_forward = []
        deltaU_backward = []
        with open(fepout_file, 'r') as input_fepout:
            # The first window samples forward only, the last window backward only
            first_window = True
            last_window = False
            while True:
                line = input_fepout.readline()
                if not line:
                    break
                
                if line.startswith('#NEW FEP WINDOW: LAMBDA SET TO 1') and windows != []:
                    last_window = True
                if line.startswith('#NEW FEP WINDOW: LAMBDA SET TO 0 ') and windows != []:
                    last_window = True
                
                if not last_window:
                    if line.startswith('#NEW FEP WINDOW:'):
                        splitedLine = line.strip().split()
                        windows.append((float(splitedLine[6]), float(splitedLine[8])))
                        continue
                
                if line.startswith('#STARTING COLLECTION'):
                    # collecting deltaU
                    deltaU_forward_per_window = []
                    deltaU_backward_per_window = []
                    while True:
                        line = input_fepout.readline()
                        if line.startswith('FepEnergy:'):
                            splitedLine = line.strip().split()
                            if not last_window:
                                deltaU_forward_per_window.append(float(splitedLine[6]))
                            else:
                                deltaU_backward_per_window.append(float(splitedLine[6]))
                        elif line.startswith('FepE_back:'):
                            splitedLine = line.strip().split()
                            deltaU_backward_per_window.append(float(splitedLine[6]))
                        else:
                            if deltaU_forward_per_window != []:
                                deltaU_forward.append(np.array(deltaU_forward_per_window))
                            if deltaU_backward_per_window != []:
                                deltaU_backward.append(np.array(deltaU_backward_per_window))
                            break
                
                if line.startswith('#Free energy change for lambda window'):
                    first_window = False
                    last_window = False
        
        if len(deltaU_forward) != len(deltaU_backward):
            raise RuntimeError('Forward and backward data do not match!')
        
        deltaU = []
        for i in range(len(deltaU_forward)):
            deltaU.append((np.array(deltaU_forward[i]), np.array(deltaU_backward[i])))
        
        return windows, deltaU
    
    def _pair_bidirectionalData(
        self, 
        forward_data: Tuple[List[Tuple[float, float]], List[NDArray]], 
        backward_data: Tuple[List[Tuple[float, float]], List[NDArray]]
    ) -> Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]:
        """ pair data from bidirectional simulations

        Args:
            forward_data (Tuple[List[Tuple[float, float]], List[NDArray]]): data from forward simulation
            backward_data (Tuple[List[Tuple[float, float]], List[NDArray]]): data from backward simulation
        
        Returns:
            Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]: recording the boundary and 
                                                                             deltaU of each window of
                                                                             bidirectional simulations 
        """
                
        merged_data = []
        
        for i in range(len(forward_data[0])):
            for j in range(len(backward_data[0])):
                if forward_data[0][i][0] == backward_data[0][j][1] and \
                    forward_data[0][i][1] == backward_data[0][j][0]:
                    merged_data.append((forward_data[1][i], backward_data[1][j]))
                    break
            else:
                raise RuntimeError('Error! the forward and backward files do not match!')
        
        return forward_data[0], merged_data
    
    def get_data(self) -> Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]:
        """ return the boundary and deltaU of each window

        Returns:
            Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]: recording the boundary and 
                                                                             deltaU of each window of
                                                                             bidirectional simulations 
        """
        
        return self._windows, self._deltaU_data

class ColvarsParser:
    """ parse Colvars cvtrj files and get the necessary data
    """
    def __init__(self, 
            cvtrj_file: str, 
            step_per_window: int, 
            equilibration_per_window: int,
            force_constants: List[int],
            centers: List[int],
            lambda_list: List[float]) -> None:
        
        self._windows, self._deltaU_data = self._read_double_wide_cvtrj(
                                                    cvtrj_file,
                                                    step_per_window, 
                                                    equilibration_per_window,
                                                    force_constants,
                                                    centers,
                                                    lambda_list
                                                )
            
        if len(self._windows) != len(self._deltaU_data):
            print(len(self._windows))
            print(len(self._deltaU_data))
            raise RuntimeError('Internal numbers of windows and deltaU do not match! This is a bug!')
        
        self._restaint_contribution = \
            self._alchemicalRestraintContributionBulk(centers[0], centers[3], centers[5], *force_constants)
        
    def _read_double_wide_cvtrj(
            self, 
            cvtrj_file: str, 
            step_per_window: int, 
            equilibration_per_window: int,
            force_constants: List[int],
            centers: List[int],
            lambda_list: List[float]
        ) -> Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]:
        """ Read an Colvars cvtrj file and regard it as double-wide free-energy calculation. 
            Return the window and bidirectional deltaU information.

        Args:
            cvtrj_file (str): the path of the fepout file
            step_per_window (int): total steps of a window 
            equilibration_per_window (int): steps of equilibration of a window
            force_constants (List[int]) force constants of each CV.
            centers (List[int]): the center (lambda=1) of each CV. The first N CVs are considered.
            lambda_list (List[float]): Lambda schedule [0, ..., 1] or [1, ..., 0]
            
        Returns:
            Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]: recording the boundary and 
                                                                             deltaU of each window of
                                                                             bidirectional simulations
        """

        assert len(force_constants) == len(centers), "Error, The lengths of force_constants of centers are not equal!"

        force_constants = np.array(force_constants)
        centers = np.array(centers)
        lambda_list = np.array(lambda_list)
        
        windows = []
        for i in range(len(lambda_list) - 1):
            windows.append((lambda_list[i], lambda_list[i + 1]))

        deltaU_forward = []
        deltaU_backward = []

        num_CVs = len(force_constants)
        num_windows = len(lambda_list)

        with open(cvtrj_file, 'r') as input_fepout:
            # The first window samples forward only, the last window backward only
            window_index = 0

            while True:
                line = input_fepout.readline()
                if not line:
                    break
                if line.startswith("#"):
                    continue
                splitedLine = line.strip().split()
                step = int(splitedLine[0])
                
                if step % step_per_window == 0:
                    # collecting deltaU
                    deltaU_forward_per_window = []
                    deltaU_backward_per_window = []

                    while True:
                        line = input_fepout.readline()
                        if not line:
                            break
                        if line.startswith("#"):
                            continue
                        splitedLine = line.strip().split()
                        step = int(splitedLine[0])

                        if step % step_per_window == 0:
                            if window_index == 0:
                                deltaU_forward.append(deltaU_forward_per_window)
                            elif window_index == num_windows - 1:
                                deltaU_backward.append(deltaU_backward_per_window)
                            else:
                                deltaU_forward.append(deltaU_forward_per_window)
                                deltaU_backward.append(deltaU_backward_per_window)
                            window_index += 1
                            break

                        # equilibration
                        if step < window_index * step_per_window + equilibration_per_window:
                            continue
                        else:
                            dis = abs(np.array(splitedLine[1:1+num_CVs]).astype(float) - centers)
                            dis[dis>180] -= 360
                            if window_index == 0:
                                deltaU_forward_per_window.append(0.5 * np.sum((lambda_list[window_index + 1] - lambda_list[window_index]) * force_constants * (dis**2)))
                            elif window_index == num_windows - 1:
                                deltaU_backward_per_window.append(0.5 * np.sum((lambda_list[window_index - 1] - lambda_list[window_index]) * force_constants * (dis**2)))
                            else:
                                deltaU_forward_per_window.append(0.5 * np.sum((lambda_list[window_index + 1] - lambda_list[window_index]) * force_constants * (dis**2)))
                                deltaU_backward_per_window.append(0.5 * np.sum((lambda_list[window_index - 1] - lambda_list[window_index]) * force_constants * (dis**2)))
        
        if len(deltaU_forward) != len(deltaU_backward):
            raise RuntimeError('Forward and backward data do not match!')
        
        deltaU = []
        for i in range(len(deltaU_forward)):
            deltaU.append((np.array(deltaU_forward[i]), np.array(deltaU_backward[i])))
        
        return windows, deltaU
    
    def get_data(self) -> Tuple[List[Tuple[float, float]], List[Tuple[NDArray, NDArray]]]:
        return self._windows, self._deltaU_data
    
    def get_restraint_contribution(self) -> float:
        return self._restaint_contribution
    
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

        contribution = BOLTZMANN * 300 * math.log(
            8 * (math.pi**2) * CSTAR / ((R**2) * math.sin(eulerTheta) * math.sin(polarTheta)) * \
            math.sqrt(forceConstantTheta * forceConstantPhi * forceConstantPsi * forceConstanttheta * \
            forceConstantphi * forceConstantR ) / ((2 * math.pi * BOLTZMANN * 300)**3)
        )
        return contribution

class FEPAnalyzer:
    """ Analyze FEP simulations
    """
    
    def __init__(
        self, 
        window_boundaries: List[Tuple[float, float]], 
        deltaU_data: List[Tuple[NDArray, NDArray]],
        temperature: float,
    ) -> None:
        """_summary_

        Args:
            window_boundaries (List[Tuple[float, float]]): boundaries of each window
            deltaU_data (List[Tuple[NDArray, NDArray]]): deltaU for forward and backward simulations
                                                         of each window. The two deltaU should be opposite
                                                         numbers.
            temperature (float): temperature of the simulation
        """
        self._windows, self._deltaU_data = window_boundaries, deltaU_data
        self._temperature = temperature
        if len(self._windows) != len(self._deltaU_data):
            raise RuntimeError('Internal numbers of windows and deltaU do not match! This is a bug!')
        
    def MergeData(
            self, 
            window_boundaries: List[Tuple[float, float]], 
            deltaU_data: List[Tuple[NDArray, NDArray]]
    ) -> bool:
        """ Merge the exist dU with another one. Used in FEP + Colvars_FEP tasks.
            window boundaries and the number of deltaU in a window must be the same

        Args:
            window_boundaries (List[Tuple[float, float]]): boundaries of each window
            deltaU_data (List[Tuple[NDArray, NDArray]]): deltaU for forward and backward simulations
                                                         of each window. The two deltaU should be opposite
                                                         numbers.

        Returns:
            bool: Whether merge is successful
        """
        if len(window_boundaries) != len(self._windows):
            print(1)
            return False
        
        #for i in range(len(window_boundaries)):
        #    print(f'FEP: {self._windows[i]}   Colvars: {window_boundaries[i]}')
        #    print(f'FEP: {len(deltaU_data[i][0])}   Colvars: {len(self._deltaU_data[i][0])}')

        for i in range(len(window_boundaries)):
            if (len(deltaU_data[i][0]) - 1) != len(self._deltaU_data[i][0]) \
                and (len(deltaU_data[i][0]) - 1) != len(self._deltaU_data[i][0]) / 2 \
                and len(deltaU_data[i][0]) != len(self._deltaU_data[i][0]) \
                and len(deltaU_data[i][0]) != len(self._deltaU_data[i][0]) / 2:
                return False
            if (len(deltaU_data[i][1]) - 1) != len(self._deltaU_data[i][1]) \
                and (len(deltaU_data[i][1]) - 1) != len(self._deltaU_data[i][1]) / 2 \
                and len(deltaU_data[i][1]) != len(self._deltaU_data[i][1]) \
                and len(deltaU_data[i][1]) != len(self._deltaU_data[i][1]) / 2:
                return False
            if (len(deltaU_data[i][0]) - 1) == len(self._deltaU_data[i][0]):
                temp_forward = self._deltaU_data[i][0] + deltaU_data[i][0][:-1]
            elif (len(deltaU_data[i][0]) - 1) == len(self._deltaU_data[i][0]) / 2:
                temp_forward = self._deltaU_data[i][0][1::2] 
                temp_forward += deltaU_data[i][0][:-1]
            elif len(deltaU_data[i][0]) == len(self._deltaU_data[i][0]):
                temp_forward = self._deltaU_data[i][0] + deltaU_data[i][0]
            elif len(deltaU_data[i][0]) == len(self._deltaU_data[i][0]) / 2:
                temp_forward = self._deltaU_data[i][0][1::2] 
                temp_forward += deltaU_data[i][0]

            if (len(deltaU_data[i][1]) - 1) == len(self._deltaU_data[i][1]):
                temp_backward = self._deltaU_data[i][1] + deltaU_data[i][1][:-1]
            elif (len(deltaU_data[i][1]) - 1) == len(self._deltaU_data[i][1]) / 2:
                temp_backward = self._deltaU_data[i][1][::2] 
                temp_backward += deltaU_data[i][1][:-1]
            elif len(deltaU_data[i][1]) == len(self._deltaU_data[i][1]):
                temp_backward = self._deltaU_data[i][1] + deltaU_data[i][1]
            elif len(deltaU_data[i][1]) == len(self._deltaU_data[i][1]) / 2:
                temp_backward = self._deltaU_data[i][1][::2] 
                temp_backward += deltaU_data[i][1]
            self._deltaU_data[i] = (temp_forward, temp_backward)

        return True
    
    def FEP_free_energy(self) -> Tuple[List[Tuple[float, float]], List[NDArray], List[NDArray]]:
        """ Calculate and return the free-energy change using the FEP equation
            
        Returns:
            Tuple[List[Tuple[float, float]], List[NDArray], List[NDArray]]: window boundaries, free energies and errors
        """
        
        free_energies = []
        errors = []
        for i in range(len(self._windows)):
            forward_free_energy = -BOLTZMANN * self._temperature * \
                np.log(np.mean(np.exp(-self._deltaU_data[i][0] / (BOLTZMANN * self._temperature))))
            backward_free_energy = -BOLTZMANN * self._temperature * \
                np.log(np.mean(np.exp(-self._deltaU_data[i][1] / (BOLTZMANN * self._temperature))))
            free_energies.append((forward_free_energy - backward_free_energy) / 2)
            errors.append(np.abs(forward_free_energy + backward_free_energy) / np.sqrt(2))
        return self._windows, free_energies, errors
    
    def BAR_free_energy(
        self, 
        tolerance: float = 1e-6,
        block_size: int = 20,
        n_bootstrap: int = 20,
    ) -> Tuple[List[Tuple[float, float]], List[NDArray], List[NDArray]]:
        """ Calculate and return the free-energy change using the BAR estimator
        
        Args:
            tolerance (float): tolerance of the SCF. Default to 1e-6.
            block_size (int): the size of the block in block bootstrap. Default to 10.
            n_bootstrap (int): number of bootstrap samples. Default to 20.
            
        Returns:
            Tuple[List[Tuple[float, float]], List[NDArray], List[NDArray]]: window boundaries, free energies and errors
        """
        
        free_energies = []
        errors = []
        for i in range(len(self._windows)):
            dA = self._BAR_estimator(self._deltaU_data[i], tolerance)
            err = self._BAR_error_estimator(self._deltaU_data[i], tolerance, block_size, n_bootstrap)
            
            free_energies.append(dA)
            errors.append(err)
        return self._windows, free_energies, errors
    
    def Window_boundaries(self) -> List[Tuple[float, float]]:
        """ Get the boundaries of windows

        Returns:
            List[Tuple[float, float]]: windows boundaries
        """
        return self._windows
    
    def _BAR_estimator(self, deltaU: Tuple[NDArray, NDArray], tolerance: float = 1e-6) -> float:
        """ Estimate the free energy of a window using the BAR estimator

        Args:
            Tuple[NDArray, NDArray]: deltaU data of forward and backward simulations
            tolerance (float): tolerance of the SCF. Default to 1e-6.

        Returns:
            float: free energy change
        """
        
        def fermi(x):
            return 1 / (1 + np.exp(x))
        
        beta = 1 / (BOLTZMANN * self._temperature)
        c = 0
        # BAR estimator
        exp_beta_dA = np.mean(fermi(beta * (deltaU[0] - c))) / np.mean(fermi(beta * (deltaU[1] + c)))
        dA = np.log(exp_beta_dA) / (-beta) + c
        
        while np.abs(c - dA) > tolerance:
            c = dA
            # BAR estimator
            exp_beta_dA = np.mean(fermi(beta * (deltaU[0] - c))) / np.mean(fermi(beta * (deltaU[1] + c)))
            dA = np.log(exp_beta_dA) / (-beta) + c
        
        return dA
    
    def _BAR_error_estimator(
        self, 
        deltaU: Tuple[NDArray, NDArray], 
        tolerance: float = 1e-6,
        block_size: int = 20,
        n_bootstrap: int = 20
    ) -> float:
        """ Estimate the error of the free energy estimate of a window 
            using the BAR estimator and block bootstrap method

        Args:
            Tuple[NDArray, NDArray]: deltaU data of forward and backward simulations
            tolerance (float): tolerance of the SCF. Default to 1e-6.
            block_size (int): the size of the block in block bootstrap. Default to 10.
            n_bootstrap (int): number of bootstrap samples. Default to 20.

        Returns:
            float: error of the free-energy estimator
        """
        
        forward_size = len(deltaU[0])
        backward_size = len(deltaU[1])
        bootstrap_samples = int(np.max((forward_size, backward_size)) / block_size)
        
        if bootstrap_samples < 1:
            raise RuntimeError('Error! block_size larger than sample size!')
        
        # block bootstrap
        estimates = np.zeros(n_bootstrap)

        for i in range(n_bootstrap):
            forward_bootstrap = np.zeros(bootstrap_samples * block_size, dtype=int)
            for idx, j in enumerate(np.random.randint(0, forward_size - block_size - 1, bootstrap_samples)):
                forward_bootstrap[idx*block_size:idx*block_size+block_size] = j + np.arange(block_size)
                    
            backward_bootstrap = np.zeros(bootstrap_samples * block_size, dtype=int)
            for idx, j in enumerate(np.random.randint(0, backward_size - block_size - 1, bootstrap_samples)):
                backward_bootstrap[idx*block_size:idx*block_size+block_size] = j + np.arange(block_size)
            
            estimates[i] = self._BAR_estimator(
                            (deltaU[0][forward_bootstrap], deltaU[1][backward_bootstrap]),
                            tolerance
                        )

        return np.std(estimates)
