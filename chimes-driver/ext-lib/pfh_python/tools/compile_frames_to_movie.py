#!/usr/bin/env python3
from __future__ import print_function,division

import sys
from subprocess import call
from optparse import OptionParser 
from tempfile import mkdtemp
from os.path import abspath
from glob import glob
from shutil import rmtree
import tqdm
import os

parser = OptionParser(usage="%prog <input name for glob> <outname> [options]")

parser.add_option('--PngDir', dest='PngDir', help="directory to store the PNGs in, if converting from PDFs.  Otherwise, puts them in a temp directory, so they don't get saved.", default=None)
parser.add_option('-r','--rate',dest='r',help="Frame rate.  Default = %default",default=16)
parser.add_option('-s','--size',dest='s',help="Output size.  Default = %default to use input sizes",default=None)
parser.add_option('-q','--quality',dest='q',help="Quality.  Lower = better.  Default = %default",default=None)
parser.add_option('--crop',dest='crop', help="Pixels to crop off the outside.  Default = %default",default=None, type=int)
parser.add_option('--downsat',dest='downsat',action='store_true',help="Lower the saturation via the hardcoded filter y=y*0.75",default=False)
parser.add_option('--downcolor',dest='downcolor',action='store_true',help="Lower the colors via the hardcoded filer r*=0.7:g*=0.8:b*=0.9",default=False)
parser.add_option('--vf',dest='vf',help="Video filters.  Separate different filters by spaces.",default=None)  #  Probably most commonly something like:  crop='out_w=in_w-2:out_h=in_h-2:x=1:y=1' lutrgb='r=val*0.75:g=val*.75:b=val*.8'",default=None)
parser.add_option('--crf', dest='crf', help="constant rate factor; varies the bitrate.  lower = higher bitrate", default=None)

ops,args = parser.parse_args()
try:
    inname,outname = args
except ValueError:
    parser.print_help()
    sys.exit(1337)
    
tempdir = mkdtemp()
files = glob(inname+'*.png')
if len(files) == 0:
    files = glob(inname+'*.jpg')
if len(files) == 0:
    files = glob(inname+'*.jpeg')
if len(files) == 0:
    files = glob(inname+'*.pdf')
if len(files) == 0:
    print("No png, jpg, jpeg, or pdf files found!")
    sys.exit(1)
files.sort()
ext = files[0].split('.')[-1]

if ext == 'pdf':
    print("..converting the PDFs to PNGs...")
    if ops.s is None:
        density = '2048'
    else:
        density = ops.s.split('x')[0]
  
    ext = 'png'
    if ops.PngDir is not None:
        os.makedirs(ops.PngDir, exist_ok=True)
        for ii, fname in enumerate(tqdm.tqdm(files)):
            ## convert to a permanent directory, if it doesn't already exist
            png_out = ops.PngDir + '/' + fname.split('/')[-1].split('.')[0] + '.' + ext
            if not os.path.isfile(png_out):
                call(['convert', '-density', density, '-quality', '100', fname, png_out])
            else:
                print("-- cowardly failing to overwrite {}".format(png_out))
                
            ## now link the file:
            call(['ln', '-s', abspath(png_out), tempdir+'/frame_{:05}.{}'.format(ii, ext)]) 
    else:
        for ii, fname in enumerate(tqdm.tqdm(files)):
            ## just convert the files to a temp directory
            call(['convert', '-density', density, '-quality', '100', fname, ops.PngDir+'/frame_{:05}.{}'.format(ii, ext)])
else:
    # just have to link the files
    for ii in range(len(files)):
        call(['ln','-s',abspath(files[ii]),tempdir+'/frame_{:05}.'.format(ii)+ext])

tocall = ['avconv','-f','image2','-r',str(ops.r),'-i',str(tempdir+'/frame_%05d.'+ext), '-threads','auto']   #,'-q',str(ops.q)]
if ops.s is not None:
    tocall.append('-s')
    tocall.append(str(ops.s))
if ops.q is not None:
   tocall.append('-qscale')
   tocall.append(str(ops.q))
if ops.crop is not None:
   if 'crop' in ops.vf:
      raise ValueError("Asked to crop in both the filters and via crop.  Use one or the other")
   if ops.vf is None:      ops.vf = ''
   ops.vf += ' crop=out_w=in_w-{0}:out_h=in_h-{0}:x={0}:y={0}'.format(ops.crop)
if ops.downsat:
   if ops.vf is None:
      ops.vf = ''
   ops.vf += ' lutyuv=y=0.75*val'
if ops.downcolor:
   if ops.vf is None:
      ops.vf = ''
   ops.vf += ' lutrgb=r=0.75*val:g=0.7*val:b=0.9*val'

if ops.vf is not None:
   for filter in ops.vf.split():
      tocall.append('-vf')
      tocall.append(filter)

if ops.crf is not None:
  tocall.append('-crf')
  tocall.append(ops.crf)

if not outname.endswith('.mp4'):
    outname += '.mp4'
tocall.append(outname)

hashes = '########################################################'
print(hashes)
print("-- running:")
print("\t", ' '.join(tocall))
print(hashes)

e = call(tocall)

rmtree(tempdir)

sys.exit(e)
