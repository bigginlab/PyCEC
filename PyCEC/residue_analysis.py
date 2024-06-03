import MDAnalysis as mda

class ResidueAnalysis:
    """
    Class to analyse the residues of the CEC system.
    """

    def __init__(self, universe, resid):
        # Attributes
        self.u = universe
        self.resid = resid
        self.times = self.get_traj_times()
        self.n_frames = len(self.u.trajectory)



    def get_traj_times(self):
        """
        Get the time of each frame in the trajectory.
        
        """
        times = []
        for ts in self.u.trajectory:
            times.append(ts.time)

        return times


    def get_torsions(self, dihedral=[], verbose=False):
        """
        Get the torsion angles of the selection.
        
        """
        # Selection
        selection = self.u.select_atoms(f'resid {self.resid} and name {' '.join([f"{atomtype}" for atomtype in dihedral])}')

        torsions = []
        for ts in self.u.trajectory:
            torsions.append(selection.dihedral.value())
        
        if verbose:
            # Combine the index numbers, times, and torsions
            torsions_h = []
            for i, t, d in zip(range(self.n_frames), self.times, torsions):
                torsions_h.append([i, t, d])

            # Print the output out line by line
            for i, t, d in torsions_h:
                print(f"Frame {i} at time {t} ps is {d:.1f} deg.")
        else:
            return torsions
        
    def get_two_torsions(self, dihedral1=[], dihedral2=[], verbose=False):
        """
        Get the torsion angles of the selection.
        
        """
        # Selection
        selection1 = self.u.select_atoms(f'resid {self.resid} and name {' '.join([f"{atomtype}" for atomtype in dihedral1])}')
        selection2 = self.u.select_atoms(f'resid {self.resid} and name {' '.join([f"{atomtype}" for atomtype in dihedral2])}')

        torsions1 = []
        torsions2 = []
        for ts in self.u.trajectory:
            torsions1.append(selection1.dihedral.value())
            torsions2.append(selection2.dihedral.value())
        
        if verbose:
            # Combine the index numbers, times, and torsions
            torsions_h = []
            for i, t, d1, d2 in zip(range(self.n_frames), self.times, torsions1, torsions2):
                torsions_h.append([i, t, d1, d2])

            # Print the output out line by line
            for i, t, d1, d2 in torsions_h:
                print(f"Frame {i} at time {t} ps is {d1:.1f} deg and {d2:.1f} deg.")
        else:
            return torsions1, torsions2

    def get_rmsd(self, reference):
        """
        Get the RMSD of the selection.
        
        """
        # Selection
        selection = self.u.select_atoms(f'resid {self.resid}')

        rmsd = []
        for ts in self.u.trajectory:
            rmsd.append(selection.rmsd(reference))
        
        return rmsd


if __name__ == "__main__":

    # Directory and title
    dir1 = '/biggin/b222/catz0163/pept/dynamics/pept_holo/pept_AF_H87P_D342P_v2/qmmm'
    title1 = 'PepT2 with AF H87P D342P'

    # Load the universe
    u1 = mda.Universe(f'{dir1}/prod-s100.pdb', f'{dir1}/prod-s100.xtc')


    # Residue analysis
    d342 = ResidueAnalysis(u1, 342)
    #d342_torsions = d342.get_torsions(['N', 'CA', 'CB', 'CG'], verbose=True)
    #d342_torsions_c = d342.get_torsions(['CA', 'CB', 'CG', 'OD1'], verbose=True)
    d342_torsions = d342.get_two_torsions(['N', 'CA', 'CB', 'CG'], ['CA', 'CB', 'CG', 'OD1'], verbose=True)