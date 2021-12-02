import os, sys
import argparse

date_object = datetime.date.today()
# define the name of the directory to be created
path = "data_" + str(date_object)

try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed. Might already exist!" % path)
else:
    print ("Successfully created the directory %s " % path)

# Functions to save and load python type objects to file using pickle.
def save_obj(data, filename ):
    with open(filename + '.pkl', 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)

def load_obj(filename ):
    with open(filename + '.pkl', 'rb') as f:
        return pickle.load(f)

def parse_args(string=None):
    """
    parse arguments
    """
    parser = argparse.ArgumentParser(description='Quantum Circuit Simulation')
    parser.add_argument('--filename', dest='filename',
                        help='the filename of the wavefunction'
                        'Default: None',
                        default=None, type=str)
    parser.add_argument('--L', dest='L',
                        help='system size. Default: 10',
                        default=10, type=int)
    parser.add_argument('--option', dest='option',
                        help='option to choose from parameter text file',
                        default=0, type=int)
    parser.add_argument('--g', dest='g',
                        help='system size. Default: 10',
                        default=1., type=float)
    parser.add_argument('--h', dest='h',
                        help='system size. Default: 10',
                        default=0., type=float)

    if len(sys.argv) == 1:
        pass
        # parser.print_help()
        # sys.exit(1)

    if string is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(string.split())
    return args


args = parse_args()

option = args.option
filename = args.filename

with open(filename) as f:
    line = f.readlines()[option-1]

args = parse_args(line)
L = args.L
g = args.g
h = args.h

print(g)
print(h)

#print(L)
#print(filename)