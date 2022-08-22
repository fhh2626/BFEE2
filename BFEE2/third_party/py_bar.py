from typing import List, Literal, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

# boltzmann constant
BOLTZMANN = 0.0019872041

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
                
                if line.startswith('#NEW FEP WINDOW: LAMBDA SET TO 1'):
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
            errors.append((forward_free_energy + backward_free_energy) / np.sqrt(2))
        return self._windows, free_energies, errors
    
    def BAR_free_energy(
        self, 
        tolerance: float = 1e-6,
        block_size: int = 20,
        n_bootstrap: int = 20
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
        
        # bootstrap
        estimates = []
        for i in range(n_bootstrap):
            forward_bootstrap = []
            for j in np.random.randint(0, forward_size - block_size - 1, bootstrap_samples):
                for k in range(block_size):
                    forward_bootstrap.append(j + k)
                    
            backward_bootstrap = []
            for j in np.random.randint(0, backward_size - block_size - 1, bootstrap_samples):
                for k in range(block_size):
                    backward_bootstrap.append(j + k)
            
            estimates.append(
                self._BAR_estimator(
                    (deltaU[0][forward_bootstrap], deltaU[1][backward_bootstrap]),
                    tolerance
                )
            )
        return np.std(estimates)
