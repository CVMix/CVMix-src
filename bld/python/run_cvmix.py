#!/usr/bin/python
import sys, os, glob, shutil, numpy
import cvmix

# BL Parameters --
# Low  latitude: 6.50e-5, 1.15e-4, 4.5e-3, 2500
# High latitude: 7.50e-5, 0.95e-4, 4.5e-3, 2500
# Orig BL paper: 8.00e-5, 1.05e-4, 4.5e-3, 2500

# --- INPUT PARAMETERS --- #
nlev = 30
ncol = 2
ocn_depth = 5250.0

# --- BRYAN-LEWIS MIXING PARAMETERS FOR COLUMN 1 --- #
col1_vdc1 = 6.50E-5
col1_vdc2 = 1.15E-4
col1_linv = 4.50E-3
col1_dpth = 2500.0

# --- BRYAN-LEWIS MIXING PARAMETERS FOR COLUMN 2 --- #
col2_vdc1 = 7.50E-5
col2_vdc2 = 0.95E-4
col2_linv = 4.5E-3
col2_dpth = 2500.0

iface_depth = numpy.zeros([nlev+1])
visc = numpy.zeros([(nlev+1) * 2])
diff = numpy.zeros([1 * (nlev+1) * 2])

iface_depth.shape = (nlev+1)
visc.shape = (ncol, nlev+1)
diff.shape = (ncol, nlev+1, 1)

print diff.shape

for i in numpy.arange(1,nlev+1):
	iface_depth[i] = iface_depth[i-1] + (ocn_depth/nlev)

[visc, diff] = cvmix.cvmix_test(ncol, ocn_depth, col1_vdc1, col1_vdc2, col1_linv, col1_dpth, col2_vdc1, col2_vdc2, col2_linv, col2_dpth, iface_depth)
#cvmix.cvmix_test(ocn_depth, col1_vdc1, col1_vdc2, col1_linv, col1_dpth, col2_vdc1, col2_vdc2, col2_linv, col2_dpth, iface_depth, visc, diff)

print diff.shape

diff_out = open('diff.dat', 'w')
diff_out.write(repr(diff))
diff_out.close()

visc_out = open('visc.dat', 'w')
visc_out.write(repr(visc))
visc_out.close()

