# PyCEC
from PyCEC.cec_system import CECSystem

# MDAnalysis
import MDAnalysis as mda  # TODO: This needs to be fixed tbh, not efficient
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

from tqdm import tqdm

from matplotlib import pyplot as plt
import numpy as np


class CVAnalysis:
    """
    Class to analyse the collective variable of the CEC system.
    """
    def __init__(self, universe, initial_resid, target_resid,
                 other_resids=None, ligand=None, cyzone_dim=[8, 8, -8],
                 frame_n=0):
        # Attributes
        self.u = universe
        self.initial_resid = initial_resid
        self.target_resid = target_resid
        self.other_resids = other_resids
        self.ligand = ligand
        self.cyzone_dim = cyzone_dim
        self.frame_n = frame_n

        # Trajectory information
        self.n_frames = len(self.u.trajectory)
        self.timestep = self.u.trajectory.dt

        self.times = self.get_traj_times()
        self.cv = CECSystem(self.u, self.initial_resid, self.target_resid,
                            self.other_resids, self.ligand, self.cyzone_dim,
                            self.frame_n)

    def get_traj_times(self):
        """
        Get the time of each frame in the trajectory.

        """
        times = []
        for ts in self.u.trajectory:
            times.append(ts.time)

        return times

    def get_water_counts(self):
        """
        Iterate through the frames and get the water counts.

        """
        water_counts = []
        for ts in tqdm(self.u.trajectory, desc="Getting water counts..."):
            self.cv.set_frame(ts.frame)
            water_counts.append(len(self.cv.waters.atoms))

        return water_counts

    def get_water_counts_h(self):
        """
        Iterate through the frames and get water counts in a readable format.

        """
        water_counts = self.get_water_counts()

        # Combine the index numbers, times, and water counts efficiently
        water_c = []
        for i, t, w in zip(range(self.n_frames), self.times, water_counts):
            water_c.append([i, t, w])

        # Print the output out line by line
        for i, t, w in water_c:
            print(f"Frame {i} at time {t} ps has {w} water molecules.")

        # return water_c


class WaterAnalysis:
    """
    Class to analyse the water molecules in the CEC system.
    """
    def __init__(self, universe, initial_resid, target_resid, other_resids,
                 ligand, cyzone_dim, frame_n):
        # Attributes
        self.u = universe
        self.initial_resid = initial_resid
        self.target_resid = target_resid
        self.other_resids = other_resids
        self.ligand = ligand
        self.cyzone_dim = cyzone_dim
        self.frame_n = frame_n  # TODO: is frame relevant here?

        # Trajectory information
        self.n_frames = len(self.u.trajectory)
        self.timestep = self.u.trajectory.dt

        self.times = self.get_traj_times()

        self.cv = CECSystem(self.u, self.initial_resid, self.target_resid,
                            self.other_resids, self.ligand, self.cyzone_dim,
                            self.frame_n)

        self.waters = self.cv.waters
        self.waters_sele_str = self.cv.waters_sele_str

    def get_traj_times(self):
        """
        Get the time of each frame in the trajectory.

        """
        times = []
        for ts in self.u.trajectory:
            times.append(ts.time)

        return times

    def get_hydrogen_bonds(self, donors_sel=None, hydrogens_sel=None,
                           acceptors_sel=None, distance_cutoff=3.5,
                           angle_cutoff=120.0):
        """
        Get the hydrogen bonds between the a pair of molecules.

        """
        # Create the HBA object
        hba = HBA(universe=self.u, donors_sel=donors_sel,
                  hydrogens_sel=hydrogens_sel, acceptors_sel=acceptors_sel,
                  d_h_cutoff=distance_cutoff, d_h_a_angle_cutoff=angle_cutoff)
        hba.run()  # run the analysis

        # Get counts
        counts = hba.count_by_time()

    def get_water_h_bonds(self, distance_cutoff=3.5, angle_cutoff=120.0):
        """
        Get the hydrogen bonds between the water molecules.

        """
        # Get the hydrogen bonds
        counts = self.get_hydrogen_bonds(donors_sel=self.waters_sele_str,
                                         hydrogens_sel=self.waters_sele_str,
                                         acceptors_sel=self.waters_sele_str,
                                         distance_cutoff=distance_cutoff,
                                         angle_cutoff=angle_cutoff)
        return counts

    def plot_hbond_counts(self, counts, save=False):
        """
        Plot the hydrogen bond counts.

        """
        fig, ax = plt.subplots()
        ax.plot(self.times, counts)
        ax.set_ylabel("Number of hydrogen bonds")
        ax.set_xlabel("Frame number")
        ax.set_title("Hydrogen bond counts over time")

        if save:
            fig.savefig("hbond_counts.png", dpi=300)

        plt.show()

    def get_max_hbond_counts(self, counts):
        """
        Get the frames and times of the 10 frames with the most hydrogen bonds.
        """
        # Get the frames with the most hydrogen bonds
        max_counts = sorted(counts, reverse=True)[:10]  # get 10 highest counts
        print(f"\nMax counts: {max_counts}")
        max_frames = []
        for c in max_counts:
            max_frames.append(np.where(counts == c)[0][0])
        max_times = [self.times[f] for f in max_frames]

        print(f"\nFrames with the most hydrogen bonds: {max_frames}")
        print(f"Times of the frames with the most hydrogen bonds: {max_times}")

        return max_frames, max_times

    def find_consistent_waters(self):
        """
        Find waters consistently included in the waters list.
        """
        pass

    def get_residence_times(self):
        """
        Get the residence time of the water molecules.
        """
        pass


if __name__ == "__main__":

    # Directory and title
    dir1 = '/biggin/b222/catz0163/pept/dynamics/pept_holo/pept_AF_H87P_D342P_v2/qmmm'
    dir1 = '/Users/nfo24278/Documents/dphil/proton_transfer/PyCEC/simulations/pept_AF_H87P_D342P'
    title1 = 'PepT2 with AF H87P D342P'

    # Load the universe
    u1 = mda.Universe(f'{dir1}/prod-s200.pdb', f'{dir1}/prod-s200.xtc')

    # # Initialise class
    # cv = CVAnalysis(u1, initial_resid=342, target_resid=56,
    #                 other_resids=[53, 57, 61, 622, 324, 345, 161, 60, 87],
    #                 ligand=[1, 2], cyzone_dim=[8, 8, -8],
    #                 frame_n=160)
    # water_c = cv.get_water_counts_h()
    # print(f"\n{water_c}")

    # Test the water analysis
    wa = WaterAnalysis(u1, initial_resid=342, target_resid=56,
                       other_resids=[53, 57, 61, 622, 324, 345, 161, 60, 87],
                       ligand=[1, 2], cyzone_dim=[8, 8, -8],
                       frame_n=160)

    counts = wa.get_water_h_bonds()

    wa.plot_hbond_counts(counts, save=True)
