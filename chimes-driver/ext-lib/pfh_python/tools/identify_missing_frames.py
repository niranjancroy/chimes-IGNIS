#!/usr/bin/env python3

def missing_elements(L):
    """ stolen from stack overflow at some point"""
    import numpy as np
    L = np.sort(L)
    start, end = L[0], L[-1]
    return sorted(set(range(start, end + 1)).difference(L))

def contiguous_integer_blocks(data):
    """
    stolen from 
    https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
    """
    from itertools import groupby
    from operator import itemgetter
    import numpy as np
    data = np.array(data, dtype=int)
    blocks = []
    for k, g in groupby(enumerate(data), lambda i_x: i_x[0]-i_x[1]):
        blocks.append(np.array(list(map(itemgetter(1), g))))
    return blocks

import sys
usage = "usage:  python {} <input base name> [max num = 1365] [format = pdf]".format(sys.argv[0].split('/')[-1])

try:
    inbase = sys.argv[1]
except Exception:
    print(usage)
    sys.exit(1)

max_num = 1366
if len(sys.argv) > 2:
    max_num = int(sys.argv[2])

fmt = 'pdf'
if len(sys.argv) > 3:
    fmt = sys.argv[3]

import glob

files = glob.glob(inbase+'*'+fmt)
frame_numbers = [int(fn.split('/')[-1].split('_')[1][1:]) for fn in files]
if max(frame_numbers) < max_num:
    frame_numbers += [max_num]
frame_numbers.sort()
me = missing_elements(frame_numbers)
if len(me):
    print("Missing {} frames between {} and {}:".format(len(me), frame_numbers[0], frame_numbers[-1]))
    print(me)

    #now find contiguous blocks in the frames that we have:
    blocks = contiguous_integer_blocks(frame_numbers)
    ## these are the ones that we can skip
    print("\nCopy the following as frame_skip_list in the submit script:")
    string = 'frame_skip_list = ['
    for b in blocks:
        string += 'np.arange({}, {}), '.format(b[0], b[-1]+1)
    string = string.rstrip(', ') + ']'
    print(string)


    ## alternatively, find the blocks in the missing elements and use those as the break points
    #blocks = contiguous_integer_blocks(me)
    # for b in blocks:
    #     print("{} -- {},".format(b[0], b[-1]), end=' ')
    # print('\nCopy the following to split into {} frame blocks:'.format(len(blocks)))
    # for b in blocks:
    #     print("['--frame_min={}', '--frame_max={}'], ".format(b[0], b[-1]), end='')
    print('\n')


else:
    print("Looks like you have all frames ranging from {} to {}!".format(frame_numbers[0], frame_numbers[-1]))
