import struct
import numpy as np
import matplotlib.pyplot as plt
import sys


def read_benchmark_data(filename):
    bytes = None

    with open(filename, 'rb') as f:
        bytes = f.read()

    return bytes


def parse_benchmark_data(bytes):
    bpc = 1
    bpi = 4
    bpd = 8
    
    keys = ['ntrials', 'ndims', 'size_min_index', 'size_max_index',
            'step', 'offset', 'base', 'nstats', 'typelen', 'type',
            'runtimes', 'runstats']

    data_map = {}
    
    start, end = 0*bpi, 9*bpi
    indices = struct.unpack('@{0}i'.format(9), bytes[start:end])
    
    for k, v in zip(keys, indices):
        data_map[k] = v
    
    # start, end = 0*bpi, 1*bpi
    # data_map['ntrials'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 1*bpi, 2*bpi
    # data_map['ndims'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 2*bpi, 3*bpi
    # data_map['size_min_index'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 3*bpi, 4*bpi
    # data_map['size_max_index'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 4*bpi, 5*bpi
    # data_map['step'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 5*bpi, 6*bpi
    # data_map['offset'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 6*bpi, 7*bpi
    # data_map['base'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 7*bpi, 8*bpi
    # data_map['nstats'] = struct.unpack('@i', bytes[start:end])
    
    # start, end = 8*bpi, 9*bpi
    # data_map['typelen'] = struct.unpack('@i', bytes[start:end])
    
    start = end
    end = start + data_map['typelen'] * bpc
    fmt = '@{0}s'.format(data_map['typelen'])
    data_map['type'] = struct.unpack(fmt, bytes[start:end])[0].decode()

    for k, v in data_map.items():
        print('{0:20s}: {1}'.format(k, v))

    start = end
    end = start + data_map['ndims'] * data_map['ntrials'] * bpd
    fmt = '@{0}d'.format(data_map['ndims'] * data_map['ntrials'])
    data_map['runtimes'] = struct.unpack(fmt, bytes[start:end])
    data_map['runtimes'] = np.reshape(data_map['runtimes'], (data_map['ntrials'], data_map['ndims'])) # Check if this indexing fudges numbers.
    print(data_map['runtimes'].shape)
    
    start = end
    end = start + data_map['ndims'] * data_map['nstats'] * bpd
    fmt = '@{0}d'.format(data_map['ndims'] * data_map['nstats'])
    data_map['runstats'] = struct.unpack(fmt, bytes[start:end])
    data_map['runstats'] = np.reshape(data_map['runstats'], (data_map['nstats'], data_map['ndims']))
    print(data_map['runstats'])
    
    return data_map


def plot_runstats(data_map):
    offset = data_map['offset']
    size_min_index = data_map['size_min_index']
    size_max_index = data_map['size_max_index']
    step = data_map['step']
    base = data_map['base']
    dimtype = data_map['type']
    
    if dimtype == 'exponential':
        dimensions = np.array([offset + base**(si) for si in range(size_min_index, size_max_index + 1, step)])
    else:
        dimensions = np.array([offset + si * step for si in range(size_min_index, size_max_index + 1)])
    print('DEBUG: dimensions =', dimensions)

    runstats = data_map['runstats']
    minstats = runstats[0, :]
    maxstats = runstats[1, :]
    meanstats = runstats[2, :]
    stdstats = runstats[3, :]

    fig, ax = plt.subplots(1, figsize=(9, 6))
    if dimtype == 'exponential':
        ax.loglog(dimensions, maxstats, color='red', marker='x', label='Maximum')
        ax.loglog(dimensions, meanstats, color='blue', marker='x', label='Mean ($\\mu$)')
        ax.loglog(dimensions, meanstats + stdstats, color='green', marker='x', label='$\\mu + \\sigma$')
    else:
        ax.semilogy(dimensions, maxstats, color='red', marker='x', label='Maximum')
        ax.semilogy(dimensions, meanstats, color='blue', marker='x', label='Mean ($\\mu$)')
        ax.semilogy(dimensions, meanstats + stdstats, color='green', marker='x', label='$\\mu + \\sigma$')

    ax.legend(loc='best')
    ax.set_xlabel('Largest Dimension')
    ax.set_ylabel('Runtime (Seconds)')
    ax.set_title('Integral Loop Runtimes with {0} runs'.format(data_map['ntrials']))
    
    plt.show()
    fig.savefig('runtimes.png')


if __name__ == '__main__':
    filename = 'benchmark_data.dat'

    if len(sys.argv) > 1:
        filename = sys.argv[1]

    bytes = read_benchmark_data(filename)
    data_map = parse_benchmark_data(bytes)
    plot_runstats(data_map)