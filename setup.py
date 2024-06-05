from setuptools import setup, find_packages

setup(
    name='PyCEC',
    version='0.1.0',
    author='Ronald Cvek',
    author_email='ronaldcvek@gmail.com',
    description='Python implementation of the Center of Excess Charge (CEC) collective variable.',
    url='https://github.com/roncv/PyCEC',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'MDAnalysis',
        'tqdm',
        'matplotlib'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
    ],
    python_requires='>=3.6',
    package_data={  # Optional: include data files
        'your_package_name': ['data/*']
    },
    entry_points={  # Optional: define command-line scripts
        'console_scripts': [
            'pycec=PyCEC.scripts.pycec:func_name',
        ],
    },
)