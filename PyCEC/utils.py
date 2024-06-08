"""
Various utility functions for PyCEC.
"""

import time
import os
import numpy as np
from MDAnalysis import AtomGroup

def dict_to_str(dictionary):
    """
    Convert a dictionary to a string.

    Parameters
    ----------
    dictionary : dict
        Dictionary to convert to a string.

    Returns
    -------
    str
        String representation of the dictionary.
    """
    # If the value is a list type, print the type
    unwrapped = ', '.join(('<array>' if isinstance(val, (list, np.ndarray, AtomGroup)) else str(val)) for val in dictionary.values())
    return unwrapped


# Decorator to time functions
def timeit(func):
    """
    Decorator to time functions.
    """
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f"function {func.__name__} ({dict_to_str(kwargs)}) took {end - start:.2f}s to run.")
        return result
    return wrapper


def create_dir(dir_name):
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