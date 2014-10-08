Follow these steps to install LSPs in your system.


Substitute /path/ with the folder where you want to install the program.

1. cd /path/
2. git clone https://github.com/danielemichilli/LSPs.git
3. export PYTHONPATH=/path/LSPs/src:$PYTHONPATH


Dependencies (order matters!):
- python >1.7
- cython
- numpy
- pandas
  - dateutil
  - pytz
  - bottleneck (recommended)
  - tables
    - HDF5
    - numexpr
- matplotlib
  - libpng (optional)
    - zlib
  - pyparsing (optional)
  - six (optional)

Dependencies for Cartesius users (full instructions at https://surfsara.nl/systems/cartesius/software/python) - RECOMMENDED:
 1. cd $HOME
 2. module load python
 3. mkdir pythonpackages (if it doesn't exsist!)
 4. cd pythonpackages
 5. wget https://pypi.python.org/packages/source/B/Bottleneck/Bottleneck-0.8.0.tar.gz
 6. tar xzvf Bottleneck-0.8.0.tar.gz
 7. cd Bottleneck-0.8.0
 8. python setup.py install --home=$HOME/pythonpackages
 9. cd $HOME
10. rm pythonpackages/Bottleneck-0.8.0.tar.gz
11. rm -r pythonpackages/Bottleneck-0.8.0
12. export PYTHONPATH=$HOME/pythonpackages/lib/python:$PYTHONPATH


To launch the routines through command (OPTIONAL):
1. cd path/LSPs/bin
2. chmod +x *
3. export PATH=path/LSPs/bin:$PATH


To make the export commands permanent, write those in your .bash file