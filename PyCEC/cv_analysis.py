# PyCEC
from collective_variable import CECCollectiveVariable

# MDAnalysis
import MDAnalysis as mda

from tqdm import tqdm
import time




class CVAnalysis:
    """
    Class to analyse the collective variable of the CEC system.
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
        self.frame_n = frame_n

        # Trajectory information
        self.n_frames = len(self.u.trajectory)
        self.timestep = self.u.trajectory.dt

        self.times = self.get_traj_times()
        self.cv = CECCollectiveVariable(self.u, self.initial_resid, self.target_resid,
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
        for ts in tqdm(self.u.trajectory, desc="Getting water counts..."): # tqdm progress bar
            self.cv.set_frame(ts.frame)
            water_counts.append(len(self.cv.waters.atoms))

        return water_counts
    
    def get_water_counts_h(self):
        """
        Iterate through the frames and get the water counts in a readable format.

        """
        water_counts = self.get_water_counts()

        # Combine the index numbers, times, and water counts efficiently
        water_c = []
        for i, t, w in zip(range(self.n_frames), self.times, water_counts):
            water_c.append([i, t, w])

        # Print the output out line by line
        for i, t, w in water_c:
            print(f"Frame {i} at time {t} ps has {w} water molecules.")

        #return water_c




# #print(dir(CECCollectiveVariable))

# frame_test = 160

# # Directory and title
# dir1 = '/biggin/b222/catz0163/pept/dynamics/pept_holo/pept_AF_H87P_D342P_v2/qmmm'
# title1 = 'PepT2 with AF H87P D342P'

# # Load the universe
# u1 = mda.Universe(f'{dir1}/prod-s100.pdb', f'{dir1}/prod-s100.xtc')

# print("\n")
# # Initialise class
# cv = CECCollectiveVariable(u1, initial_resid=342, target_resid=56,
#                             other_resids=[53, 57, 61, 622],
#                             ligand=[1, 2], cyzone_dim=[8, 8, -8],
#                             frame_n=frame_test)

# # The atom groups
# # print(cv.cv_atom_group)
# # print(cv.qm_all)

# cv.set_frame(100)
# cv.get_sele_info()
# print("\n")
# cv.set_frame(140)
# cv.get_sele_info()
# print("\n")
# cv.set_frame(140)


if __name__ == "__main__":

    # Directory and title
    dir1 = '/biggin/b222/catz0163/pept/dynamics/pept_holo/pept_AF_H87P_D342P_v2/qmmm'
    title1 = 'PepT2 with AF H87P D342P'

    # Load the universe
    u1 = mda.Universe(f'{dir1}/prod-s100.pdb', f'{dir1}/prod-s100.xtc')

    # Initialise class
    cv = CVAnalysis(u1, initial_resid=342, target_resid=56,
                    other_resids=[53, 57, 61, 622, 324, 345, 161, 60, 87],
                    ligand=[1, 2], cyzone_dim=[8, 8, -8],
                    frame_n=160)
    water_c = cv.get_water_counts_h()
    print(f"\n{water_c}")