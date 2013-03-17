#!/usr/bin/env python

import argparse as arg
import matplotlib.mlab as mlab

parser = arg.ArgumentParser('Calculate MRS spectra from P files')

parser.add_argument('in_file', action='store', metavar='File', 
                    help='MRS P file')

parser.add_argument('out_file', action='store', metavar='File', 
                    help='Output file (.csv)')
 
parser.add_argument('--sampling_rate', action='store', metavar='Float',
                    help='Sampling rate (default : 5000.0 Hz)',
                    default=5000.0)

parser.add_argument('--min_ppm', action='store', metavar='Float',
                    help='The minimal frequency (default: -0.7 ppm)',
                    default=-0.7)

parser.add_argument('--max_ppm', action='store', metavar='Float',
                    help='The minimal frequency (default: 4.3 ppm)',
                    default=4.3)

in_args = parser.parse_args()

if __name__ == "__main__":
    # This reads from file and does spectral analysis everything:
    f_ppm, e1, e2 = quant.reconstruct(in_args.in_file,
                                      in_args.max_ppm,
                                      in_args.min_ppm,
                                      in_args.sampling_rate)
    
    # Pack it into a recarray:
    names = ('ppm', 'echo1', 'echo2', 'diff', 'sum')
    formats = (float, float, float, float, float)
    dt = zip(names, formats)
    m_e1 = np.mean(e1, 0)
    m_e2 = np.mean(e2, 0)
    diff = m_e2 - m_e1
    sum = m_e2 + m_e1
    prep_arr = [(f_ppm[i], m_e1[i], m_e2[i], diff[i], sum[i]) for i in
                range(len(f_ppm))]

    out_array = np.array(prep_arr, dtype=dt)

    # And save to output:
    mlab.rec2csv(out_array, in_args.out_file)
    
    
