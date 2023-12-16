from setuptools import setup

setup(
    name='pynfam',

    version='2.0',

    url='https://bitbucket.org/evney/pynfam/src',

    author='Evan Ney',

    author_email='evan.ney@unc.edu',

    description='An MPI enabled Python workflow for large scale nuclear beta decay calculations.',

    long_description=open('README.md').read(),

    packages=[
        'pynfam',
        'pynfam.utilities',
        'pynfam.outputs',
        'pynfam.fortran',
        'pynfam.strength'
        ],

    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX',
        'License :: OSI Approved :: MIT License'
        ],

    install_requires=[
        'numpy',
        'pandas',
        'mpi4py',
        'f90nml',
        'scipy',
        'future'
        ],

    python_requires='>=2.7'
)
