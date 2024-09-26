"""
PyCEC.collective_variable is a module that contains the CECCollectiveVariable\
class, which is used to generate and analyse a collective variable for the\
CEC system. The class is used to select a set of residues from a\
MDAnalysis Universe object that define a collective variable. The class\
can also be used to write the selection to a PDB file, generate a GROMACS\
compatible index file, and get information about the selection.
"""
# LOCAL
if __name__ == "__main__":
    import utils
else:
    from PyCEC import utils


# CORE
import os

# PLOTTING
import matplotlib.pyplot as plt

# MDANALYISIS
import warnings  # MDAnalysis warnings are annoying
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")

    import MDAnalysis as mda

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"


# CLASS: CECCollectiveVariable
class CECSystem:
    """
    Class to generate and analyse a collective variable.
    """
    # ATTRIBUTES
    # Atom and residue types
    water_types_t = ('TIP3', 'SOL')
    heavy_atom_types_t = ('O*', 'N*')
    light_atom_types_t = ('HD*', 'HE*')

    def __init__(self, universe, initial_resid, target_resid,
                 other_resids=None, backbone=False, ligand=None,
                 ligand_bb=False, cyzone_dim=[8, 8, -8], frame_n=0,
                 warnings_on=True):
        """
        Initialise the class.

        Parameters
        ----------
        universe : MDAnalysis Universe
            Universe object to select residues from.
        initial_resid : int
            Residue id of the protonated residue (CEC start state).
        target_resid : int
            Residue id of the deprotonated residue (CEC end state).
        other_resids : list of int, optional
            List of residue ids to include in selection. (default is None)
        backbone : bool, optional
            Include backbone atoms in the selection. (default is False)
        ligand : list of str, optional
            List of ligand resids to include in selection. (default is None)
        ligand_bb : bool, optional
            Include backbone atoms of ligand in selection. (default is False)
        cyzone_dim : list of int, optional
            Dimensions of the cyzone around the protein (default is [8, 8, -8])
        frame_n : int, optional
            Frame number to write to the PDB file (default is 0)
        """
        self.universe = universe
        self.initial = initial_resid
        self.target = target_resid
        self.other_resids = other_resids
        self.backbone = backbone
        self.ligand = ligand
        self.ligand_bb = ligand_bb
        self.cyzone_dim = cyzone_dim

        self.frame_n = frame_n

        # Set the frame
        self.universe.trajectory[frame_n]

        # Generate water selection at initialisation
        self.water_types = ' or '.join([f'resname {wt}' for wt in self.water_types_t])

        # Warnings
        if warnings_on is False:
            warnings.filterwarnings('ignore')
            # warnings.filterwarnings("ignore", category=UserWarning)

        # Self-attributes
        # CV Selection
        self.cv_selection_str_all, resids_sele, lig_str = self.generate_cv_selection()
        self.cv_atom_group = self.universe.select_atoms(self.cv_selection_str_all)  # all resids
        self.resids_atom_group = self.universe.select_atoms(resids_sele)  # initial, target and other resids
        self.ligand_atom_group = self.universe.select_atoms(lig_str)  # ligands
        self.protein = self.universe.select_atoms('protein')  # protein

        # QM atoms
        self.qm_all, self.qm_heavy, self.qm_light, self.qm_all_sele_str, \
            self.qm_all_sele_br_str, self.waters, self.waters_sele_str = self.select_QM_atoms()

    def generate_cv_selection(self):
        """
        Function to select a set of residues from resIDs that define a CV.
        """

        # CV residues
        cv_resids_sele = f"(resid {self.initial} {self.target}"
        if self.other_resids is not None:
            for resid in self.other_resids:
                cv_resids_sele += f' {resid}'
            cv_resids_sele += f') and not ({self.water_types})'
        if self.backbone is False:
            # Exclude backbone atoms and ALL protons
            cv_resids_sele += ' and not backbone and not name H*'
        # print(f"cv_resids_sele: {cv_resids_sele}")

        # Ligands
        if self.ligand is not None:
            lig_str = '(resid'
            for lig in self.ligand:
                lig_str += f' {lig}'
            lig_str += f') and not ({self.water_types})'
            if self.ligand_bb is False:
                lig_str += '  and not name H*'  # remove backbone and protons
            # print(f"ligand_sele: {lig_str}")
            # Combine selections
            sel_str = f'({cv_resids_sele}) or ({lig_str})'

        # print(f"\nSelection string: {sel_str}\n")

        # Return the selection string
        return sel_str, cv_resids_sele, lig_str

    def select_QM_atoms(self):
        """
        Select the atoms and waters to be treated with QM from CV selection.
        """
        # CV residues - initial and target
        cv_resids_sele = f"resid {self.initial} {self.target}"

        # Heavy atoms selection
        ha_types_sele = ' or '.join([f'name {ha}' for ha in self.heavy_atom_types_t])
        ha_sele_str = (
            f'({cv_resids_sele}) and ({ha_types_sele}) and protein'
            f' and not backbone and not name H*'
            )  # exclude backbone and protons
        ha_sele = self.universe.select_atoms(ha_sele_str)  # ha atomgroup

        # Light atoms selection
        la_types_sele = ' or '.join([f'name {la}' for la in self.light_atom_types_t])
        la_sele_str = (
            f'({cv_resids_sele}) and ({la_types_sele}) '
            f'and protein and not backbone'
            )  # exclude backbone
        la_sele = self.universe.select_atoms(la_sele_str)  # la atomgroup

        # WATERS SELECTION
        # Heavy water atoms
        ha_wat_sele_str = (
            f'name O* and ({self.water_types}) and cyzone '
            f'{' '.join([f"{dim}" for dim in self.cyzone_dim])}'
            f' (({cv_resids_sele}) and protein)'
            )
        # print(f"ha_wat_sele_str: {ha_wat_sele_str}")
        ha_water_sele = self.universe.select_atoms(ha_wat_sele_str)  # ha water

        # Light water atoms
        la_wat_sele_str = (
            f'(byres (name O* and ({self.water_types}) and cyzone '
            f'{' '.join([f"{dim}" for dim in self.cyzone_dim])} '
            f'(({cv_resids_sele}) and protein))) and name H*'
            )
        # print(f"la_wat_sele_str: {la_wat_sele_str}")
        la_water_sele = self.universe.select_atoms(la_wat_sele_str)  # la water

        # All water atoms - atomgroup
        waters = ha_water_sele + la_water_sele
        waters_sele_str = f'({ha_wat_sele_str}) or ({la_wat_sele_str})'

        # Combine the selections
        qm_heavy = ha_sele + ha_water_sele
        qm_light = la_sele + la_water_sele
        qm_all = ha_sele + ha_water_sele + la_sele + la_water_sele
        qm_all_sele_str = f'({ha_sele_str}) or ({ha_wat_sele_str}) or ({la_sele_str}) or ({la_wat_sele_str})'
        qm_all_sele_br_str = f'byres (({ha_sele_str}) or ({ha_wat_sele_str}) or ({la_sele_str}) or ({la_wat_sele_str}))'


        # Return the selections
        return qm_all, qm_heavy, qm_light, qm_all_sele_str, qm_all_sele_br_str, waters, waters_sele_str

    def set_frame(self, frame_n=None):
        """
        Function to set the frame number.

        Parameters
        ----------
        frame_n : int, optional
            Frame number to set the universe to (default is None)
        """
        # set frame to class frame if not provided or not the same
        if frame_n is None or frame_n == self.frame_n:
            frame_n = self.frame_n
        # set the frame to the provided frame number if different
        elif frame_n != self.frame_n:
            self.universe.trajectory[frame_n]  # set the frame
            self.frame_n = frame_n  # update the frame instance attribute
            # Update qm selections  waters) every time the frame is updated
            self.qm_all, self.qm_heavy, self.qm_light, self.waters = self.select_QM_atoms()

    def set_cyzone_dim(self, cyzone_dim):
        """
        Function to set the cyzone dimensions.

        Parameters
        ----------
        cyzone_dim : list of int
            Dimensions of the cyzone around the protein.
        """
        self.cyzone_dim = cyzone_dim
        # QM atoms
        self.qm_all, self.qm_heavy, self.qm_light, self.qm_all_sele_str, \
            self.qm_all_sele_br_str, self.waters, self.waters_sele_str = self.select_QM_atoms()

    def create_dir(self, dir_name):
        """
        Function to create a directory.

        Parameters
        ----------
        dir_name : str
            Name of the directory to create.
        """
        # Create the directory
        try:
            os.mkdir(dir_name)
        except FileExistsError:
            pass

    @utils.timeit
    def write_pdb(self, filename='cv-selection', dir_name='structures',
                  atom_group=None, frame_n=None, idx=True,
                  idx_grp_nm='CV-atoms'):
        """
        Function to write a selection to a PDB file.

        Parameters
        ----------

        """
        # Set the frame
        self.set_frame(frame_n=frame_n)

        # Create the directory
        utils.create_dir(dir_name) # TODO: follow up this utils section

        # Write the selection to a PDB file
        with mda.Writer(f"./{dir_name}/{filename}-f{self.frame_n}.pdb", atom_group) as writer:
            writer.write(atom_group)

        # Write an index file for the selection
        if idx:
            self.write_index(filename=filename, dir_name=dir_name,
                             idx_grp_nm=idx_grp_nm, atom_group=atom_group,
                             frame_n=self.frame_n)

    # Just make this generalisable to any selection
    def write_index(self, filename='cv-selection', dir_name='structures',
                    idx_grp_nm='CV-atoms', atom_group=None, frame_n=None,
                    verbose=False):
        """
        Function to generate a GROMACS compatible index (.ndx) file from a \
            selection's atomsIDs.

        Parameters
        ----------

        """
        # Set the frame
        self.set_frame(frame_n=frame_n)

        # Create the directory
        self.create_dir(dir_name)

        # Get the atom IDs
        atom_ids = [f"{atom.id:4d}" for atom in atom_group]  # assume min 4 digits per element
        rows = [atom_ids[i:i + 15] for i in range(0, len(atom_ids), 15)]  # 15 elements per row

        # Prepare index group (append everything together)
        index_grp = f"[ {idx_grp_nm} ]\n\n" + "\n".join(" ".join(row) for row in rows) + "\n"

        # Write to the file
        with open(f"./{dir_name}/{filename}-f{frame_n}.ndx", 'w') as file:
            file.write(index_grp)

        # Print the same content to the console if verbose
        if verbose:
            print(index_grp)

    def write_cv_pdb(self, filename='cv-selection', dir_name='structures',
                     frame_n=None):
        """
        Function to write the the full CV selection to a PDB file.

        Parameters
        ----------

        """
        self.write_pdb(filename=filename, dir_name=dir_name,
                       atom_group=self.cv_atom_group, frame_n=frame_n,
                       idx=True, idx_grp_nm='CV-atoms')

    def write_qm_pdb(self, filename='qm-selection', dir_name='structures',
                     frame_n=None):
        """
        Function to write the the QM atom selection to a PDB file.

        Parameters
        ----------

        """
        self.write_pdb(filename=filename, dir_name=dir_name,
                       atom_group=self.qm_all, frame_n=frame_n, idx=True,
                       idx_grp_nm='QM-atoms')

    def write_protein_pdb(self, filename='protein', dir_name='structures',
                          frame_n=None):
        """
        Function to write just the protein to a PDB file.

        Parameters
        ----------

        """
        self.write_pdb(filename=filename, dir_name=dir_name,
                       atom_group=self.protein, frame_n=frame_n, idx=True,
                       idx_grp_nm='Protein')

    def write_cv_and_qm(self, filename_cv='cv-selection',
                        filename_qm='qm-selection', filename_pt='protein',
                        dir_name='structures', frame_n=None):
        """
        Function to write the CV and QM selections to PDB files.

        Parameters
        ----------

        """
        self.write_cv_pdb(filename=filename_cv, dir_name=dir_name, frame_n=frame_n)
        self.write_qm_pdb(filename=filename_qm, dir_name=dir_name, frame_n=frame_n)
        self.write_protein_pdb(filename=filename_pt, dir_name=dir_name, frame_n=frame_n)

    def get_sele_info(self, verbose=False):
        """
        Function to get information about the selection.
        """
        sele = self.cv_atom_group
        print("\n--- RESIDS ---\n")
        #print(f"Selection string: {self.cv_selection_str_all}")
        print(f"Number of atoms: {len(sele.atoms)}")
        print(f"Number of residues: {len(sele.residues)}")
        print(f"Residue names: {', '.join(set(sele.resnames))}")
        print(f"Atom names: {', '.join(set(sele.atoms.names))}")
        if verbose:
            print(f"Residue ids: {', '.join(map(str, set(sele.resids)))}")
            print(f"Atom ids: {', '.join(map(str, set(sele.indices)))}")
        print("\n--- QM ATOMS ---\n")
        print(f"Number of QM atoms: {len(self.qm_all.atoms)}")
        print(f"Number of QM residues: {len(self.qm_all.residues)}")
        print(f"Residue names: {', '.join(set(self.qm_all.resnames))}")
        print(f"Atom names: {', '.join(set(self.qm_all.atoms.names))}")
        if verbose:
            print(f"Residue ids: {', '.join(map(str, set(self.qm_all.resids)))}")
            print(f"Atom ids: {', '.join(map(str, set(self.qm_all.indices)))}")
        print("\n--- WATERS ---\n")
        print(f"Number of water atoms: {len(self.waters.atoms)}")
        print(f"Number of water residues: {len(self.waters.residues)}")
        print(f"Residue names: {', '.join(set(self.waters.resnames))}")
        print(f"Atom names: {', '.join(set(self.waters.atoms.names))}")
        if verbose:
            print(f"Residue ids: {', '.join(map(str, set(self.waters.resids)))}")
            print(f"Atom ids: {', '.join(map(str, set(self.waters.indices)))}")
        print("\n")

###############################################################################


# Testing the class
if __name__ == "__main__":

    # Provisional - Due to excessive warnings from MDAnalysis
    warnings.filterwarnings('ignore')
    frame_test = 91

    # Directory and title
    # dir1 = '/biggin/b222/catz0163/pept/dynamics/pept_holo/pept_AF_H87P_D342P_v2'
    dir1 = './simulations/pept_AF_H87P_D342P'
    title1 = 'PepT2 with AF H87P D342P'

    # Load the universe
    u1 = mda.Universe(f'{dir1}/prod-s200.pdb', f'{dir1}/prod-s200.xtc')

    print("\n")
    # Initialise class
    cv = CECSystem(u1,
                   initial_resid=342,
                   target_resid=56,
                   other_resids=[53, 57, 61, 622, 324, 345, 161, 60, 87],
                   ligand=[1, 2],
                   cyzone_dim=[9, 8, -8],
                   frame_n=frame_test)

    # cv.cv_selection_string
    cv.get_sele_info()
    cv.write_cv_and_qm(frame_n=frame_test, dir_name='../structures_342_56')

    # Initialise class
    cv2 = CECSystem(u1,
                    initial_resid=87,
                    target_resid=56,
                    other_resids=[53, 57, 61, 622, 324, 345, 161, 60, 342],
                    ligand=[1, 2],
                    cyzone_dim=[9, 8, -8],
                    frame_n=frame_test)

    # cv.cv_selection_string
    cv2.get_sele_info()
    cv2.write_cv_and_qm(frame_n=frame_test, dir_name='../structures_87_56')
