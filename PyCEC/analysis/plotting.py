from PyCEC.cec_system import CECSystem
from PyCEC.analysis.analysis import CVAnalysis

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from tqdm import tqdm

import MDAnalysis as mda


class CVPlotting:
    """
    Class to plot the collective variable of the CEC system.
    """

    def __init__(self, cv_analysis):
        # Attributes
        self.cv_analysis = cv_analysis

    def plot_water_counts(self, save=False):
        """Plot the water counts over time.

        Parameters
        ----------
        save : bool, optional
            _description_, by default False
        """        
        water_counts = self.cv_analysis.get_water_counts()

        fig, ax = plt.subplots()
        ax.plot(self.cv_analysis.times, water_counts)
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("Number of water molecules")
        ax.set_title("Number of water molecules over time")

        if save:
            fig.savefig("water_counts.png")

        plt.show()

    def plot_water_counts_line(self, save=False, filename='water_counts_line_hist', test=False, show=False):
        """
        Plot the water counts over time, with a histogram of frequency to the right of the main plot.
        
        """

        if test:
            # Generate some random data
            water_counts = np.random.randint(0, 100, self.cv_analysis.n_frames)
        else:
            # Get the water counts
            water_counts = self.cv_analysis.get_water_counts()

        # Create a figure and a grid of subplots
        fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [4, 1], 'wspace': 0.02}, figsize=(10, 5))

        # Line plot on the left
        axs[0].plot(self.cv_analysis.times, water_counts, c='b', dash_capstyle='round')
        axs[0].set_title('Collective Variable: Water Counts over Time')
        axs[0].set_xlabel('Time (ps)')
        axs[0].set_ylabel('Water Counts per Frame')
        axs[0].ticklabel_format(style='plain')

        # Histogram on the right
        axs[1].hist(water_counts, bins=self.cv_analysis.n_frames//10, orientation='horizontal', 
                    color='b', edgecolor='white', lw=0.4)
        axs[1].axis('off')

        # Adjust layout to prevent overlap
        plt.tight_layout()

        # Save the plot
        if save:
            fig.savefig("water_counts2.png", dpi=300)

        # Show the plot
        if show:
            plt.show()


    def plot_water_counts_boxplot(self, save=False):
        """
        Plot the boxplot of water counts.
        
        """
        water_counts = self.cv_analysis.get_water_counts()

        fig, ax = plt.subplots()
        ax.boxplot(water_counts)
        ax.set_ylabel("Number of water molecules")
        ax.set_title("Boxplot of water counts")

        if save:
            fig.savefig("water_counts_boxplot.png")

        plt.show()

    def plot_water_counts_violinplot(self, save=False):
        """
        Plot the violin plot of water counts.
        
        """
        water_counts = self.cv_analysis.get_water_counts()

        fig, ax = plt.subplots()
        ax.violinplot(water_counts)
        ax.set_ylabel("Number of water molecules")
        ax.set_title("Violin plot of water counts")

        if save:
            fig.savefig("water_counts_violinplot.png")

        plt.show()


if __name__ == "__main__":
    # Directory and title
    dir1 = '/biggin/b222/catz0163/pept/dynamics/pept_holo/pept_AF_H87P_D342P_v2/qmmm'
    title1 = 'PepT2 with AF H87P D342P'

    # Load the universe
    u1 = mda.Universe(f'{dir1}/prod-s200.pdb', f'{dir1}/prod-s200.xtc')

    # Initialise class
    cv_analysis = CVAnalysis(u1, initial_resid=342, target_resid=56,
                    other_resids=[53, 57, 61, 622],
                    ligand=[1, 2], cyzone_dim=[8, 8, -8],
                    frame_n=160)

    # Create the CVPlotting object
    cv_plotting = CVPlotting(cv_analysis)

    # Plot the water counts over time
    cv_plotting.plot_water_counts_line(save=True, test=True)