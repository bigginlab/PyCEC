# PyCEC
from PyCEC.cec_system import CECSystem

# MDAnalysis
import MDAnalysis as mda  # TODO: This needs to be fixed tbh, not efficient
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

from tqdm import tqdm

from matplotlib import pyplot as plt
import seaborn as sns
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
        self.all_resids = initial_resid + target_resid + other_resids
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
        self.all_resids = [initial_resid] + [target_resid] + other_resids
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
        # hba = HBA(universe=self.u, donors_sel=donors_sel,
        #           hydrogens_sel=hydrogens_sel, acceptors_sel=acceptors_sel,
        #           d_h_cutoff=distance_cutoff, d_h_a_angle_cutoff=angle_cutoff)
        hba = HBA(universe=self.u, donors_sel=donors_sel,
                  hydrogens_sel=hydrogens_sel, acceptors_sel=acceptors_sel,
                  d_h_cutoff=distance_cutoff, d_h_a_angle_cutoff=angle_cutoff)
        hba.run()  # run the analysis

        # Get counts
        counts = hba.count_by_time()

        return counts


    def get_water_h_bonds(self, distance_cutoff=3.5, angle_cutoff=120.0, waters_only=True):
        """
        Get the hydrogen bonds between the water molecules.

        """
        # Get the hydrogen bonds
        if waters_only:
            counts = self.get_hydrogen_bonds(donors_sel=self.cv.waters_sele_str,
                                             hydrogens_sel=self.cv.waters_sele_str,
                                             acceptors_sel=self.cv.waters_sele_str,
                                             distance_cutoff=distance_cutoff,
                                             angle_cutoff=angle_cutoff)
        else:
            counts = self.get_hydrogen_bonds(donors_sel=self.cv.qm_all_sele_br_str,
                                             hydrogens_sel=self.cv.qm_all_sele_br_str,
                                             acceptors_sel=self.cv.qm_all_sele_br_str,
                                             distance_cutoff=distance_cutoff,
                                             angle_cutoff=angle_cutoff)

        print(len(counts))
        return counts

    def plot_hbond_counts(self, counts, save=False, filename=None, title=None):
        """
        Plot the hydrogen bond counts.
        """
        fig, ax = plt.subplots()
        ax.plot(self.times, counts, color='black', lw=1)
        ax.set_ylabel("Number of hydrogen bonds")
        ax.set_xlabel("Frame number")
        if title:
            ax.set_title(title)
        else:
            ax.set_title("Hydrogen bond counts over time")
        ax.ticklabel_format(style='plain')

        if save:
            if filename:
                fig.savefig(filename, dpi=300)
            else:
                fig.savefig('hbond_counts.png', dpi=300)

        plt.show()

    def plot_multiple_hbond_counts_grid(
            self,
            parameter='cyzone_dim',
            ncols=2,
            values=[[9, 8, -8], [9, 8, -8]],
            save=False,
            filename=None,
            title=None,
            distance_cutoff=3.5,
            angle_cutoff=120.0,
            waters_only=True
            ):
        """
        Plot multiple hydrogen bond counts in a grid.
        """

        import math
        num_plots = len(values)
        nrows = math.ceil(num_plots / ncols)  # Calculate the required number of rows
        
        # Create subplots with nrows and ncols
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(6 * ncols, 4 * nrows),
            sharey=True
            )
        
        # Flatten axes if it's a 2D array for easier indexing
        axes = axes.flatten()
        
        for i, value in enumerate(values):
            if parameter == 'cyzone_dim':
                print(f"Setting cyzone_dim to {value}")
                self.cv.set_cyzone_dim(value)
                print(f"{len(self.cv.waters.atoms)} waters selected.")
                print(f"{len(self.cv.qm_all.atoms)} atoms selected.\n")
                counts = self.get_water_h_bonds(
                    distance_cutoff=distance_cutoff,
                    angle_cutoff=angle_cutoff,
                    waters_only=waters_only
                    )
                _, _ = self.get_max_hbond_counts(counts)
                
            else:
                pass
            ax = axes[i]
            ax.plot(self.times, counts, color='black', lw=1)
            ax.set_ylabel("Number of hydrogen bonds")
            ax.set_xlabel("Frame number")
            
            if waters_only:
                ax.set_title(f"Water H-bond counts: {parameter} = {value}, {len(self.cv.waters.atoms)} waters")
            else:
                ax.set_title(f"H-bond counts: {parameter} = {value}, {len(self.cv.qm_all.atoms)} atoms")
            ax.ticklabel_format(style='plain')
        
        # Hide any unused subplots (if num_plots < nrows * ncols)
        for j in range(num_plots, len(axes)):
            axes[j].set_visible(False)
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the figure if needed
        if save:
            if filename:
                fig.savefig(filename, dpi=300)
            else:
                fig.savefig('multiple_hbond_counts_grid.png', dpi=300)

        plt.show()

    def plot_multiple_hbond_counts_overlay(
            self,
            parameter='cyzone_dim',
            values=[[9, 8, -8], [9, 8, -8]],
            save=False,
            filename=None,
            title=None,
            distance_cutoff=3.5,
            angle_cutoff=120.0,
            waters_only=True
            ):
        """
        Plot and overlay multiple hydrogen bond counts.
        """        
        # Create subplots with nrows and ncols
        fig, ax = plt.subplots()
        all_counts = []
        for i, value in enumerate(values):
            if parameter == 'cyzone_dim':
                print(f"Setting cyzone_dim to {value}")
                self.cv.set_cyzone_dim(value)
                print(f"{len(self.cv.waters.atoms)} waters selected.")
                print(f"{len(self.cv.qm_all.atoms)} atoms selected.\n")
                counts = self.get_water_h_bonds(
                    distance_cutoff=distance_cutoff,
                    angle_cutoff=angle_cutoff,
                    waters_only=waters_only
                    )
                _, _ = self.get_max_hbond_counts(counts)
                all_counts.append(counts)
                
            else:
                pass
            if waters_only:
                ax.plot(self.times, counts, lw=1, label=f"{parameter} = {value}, {len(self.cv.waters.atoms)} waters")
            else:
                ax.plot(self.times, counts, lw=1, label=f"{parameter} = {value}, {len(self.cv.qm_all.atoms)} atoms")

        ax.legend(frameon=False)
        ax.set_ylabel("Number of hydrogen bonds")
        ax.set_xlabel("Time [ns]")
        
        if waters_only:
            ax.set_title(f"Water H-bond counts: Variable {parameter}")
        else:
            ax.set_title(f"H-bond counts: Variable {parameter}")
        ax.ticklabel_format(style='plain')
        
        
        # Adjust layout
        plt.tight_layout()
        
        # Save the figure if needed
        if save:
            if filename:
                fig.savefig(filename, dpi=300)
            else:
                fig.savefig('multiple_hbond_counts_grid.png', dpi=300)

        plt.show()

        return np.array(all_counts)
    
    def plot_standard_deviation(self, counts):
        """
        Plot the standard deviation of the hydrogen bond counts.
        """
        fig, ax = plt.subplots()
        ax.plot(self.times, np.std(counts, axis=0), color='black', lw=1)
        ax.set_ylabel("Standard deviation of hydrogen bond counts")
        ax.set_xlabel("Frame number")
        ax.set_title("Standard deviation of hydrogen bond counts over time")
        ax.ticklabel_format(style='plain')
        plt.tight_layout()
        plt.show()

    def plot_autocorrelation(self, counts):
        """
        Plot the autocorrelation of the hydrogen bond counts.
        """
        from statsmodels.graphics.tsaplots import plot_acf
        fig, ax = plt.subplots()
        plot_acf(counts, ax=ax, lags=10, zero=False, bartlett_confint=False)
        ax.set_title("Autocorrelation of hydrogen bond counts")
        plt.tight_layout()
        plt.show()

    def water_bridge_analysis(self):
        """
        Analyse the water bridges between the water molecules.
        MDAnalysis.analysis.hydrogenbonds.WaterBridgeAnalysis()
        Very slow and not super useful.
        """
        pass

    def hydrogen_bond_lifetimes(self):
        """
        Get the hydrogen bond lifetimes.
        """
        palette = sns.color_palette("flare", n_colors=len(self.all_resids))
        fig, ax = plt.subplots()
        for i, resid in enumerate(self.all_resids):

            hba = HBA(
                universe=self.u,
                donors_sel=self.cv.waters_sele_str + ' or ' + f'resid {resid}',
                hydrogens_sel=self.cv.waters_sele_str + ' or ' + f'resid {resid}',
                acceptors_sel=self.cv.waters_sele_str + ' or ' + f'resid {resid}'
                )
            hba.run()  # run the analysis

            # Get counts
            # counts = hba.count_by_time()
            tau_timeseries, timeseries = hba.lifetime(intermittency=5)
            print(f"\nHydrogen bond tau: {tau_timeseries}")
            print(f"Timeseries: {timeseries}")

            ax.plot(tau_timeseries, timeseries, color=palette[i], lw=1, label=f"Resid {resid}")

        ax.set_ylabel("Autocorrelation")
        ax.set_xlabel("Tau")
        ax.set_title("Hydrogen bond lifetimes")
        ax.ticklabel_format(style='plain')
        ax.legend(frameon=False)
        plt.tight_layout()
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

    # counts = wa.get_water_h_bonds()
    # wa.plot_autocorrelation(counts)

    # # wa.plot_hbond_counts(counts, save=True)
    # # wa.plot_multiple_hbond_counts_grid(parameter='cyzone_dim', values=[[7, 8, -8], [10, 8, -8]], save=False, waters_only=False)
    # all_counts = wa.plot_multiple_hbond_counts_overlay(parameter='cyzone_dim', values=[[7, 8, -8], [8, 8, -8], [9, 8, -8], [10, 8, -8]], save=False, waters_only=False)
    # wa.plot_standard_deviation(all_counts)

    wa.hydrogen_bond_lifetimes()