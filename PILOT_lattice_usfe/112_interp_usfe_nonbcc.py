import numpy as np
import sys
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

args=sys.argv
lattice_file=args[1] #filename for atomic lattice
N1 = int(args[2])
N2 = int(args[3])
N3 = int(args[4])

#READ IN LATTICE
print('reading in lattice file')
coords = np.genfromtxt(lattice_file,usecols=[1,2,3],skip_header=2)
elems = np.genfromtxt(lattice_file,usecols=[0],dtype=str,skip_header=2)
lattice = np.zeros(elems.shape,dtype=int)


#print(f"coords = {coords}")
print(f"elms = {elems}")
#print(f"lattice = {lattice}")

#READ IN USFE VALUES
print('reading in USFE values')
filename = args[5] if len(args) > 5 else 'MoNbTi_112_1323K.csv'  # optional 5th arg
print('using SFE data file: %s' % filename)
# Read the CSV.  The header names the composition columns after SFE, e.g.
#   SFE,Ta,W          (2 components)   or   sfe,Ti,Nb,Mo   (3 components)
# Composition columns are stored in REVERSED header order internally -- an
# arbitrary but consistent choice that keeps existing 3-component MoNbTi runs
# byte-for-byte identical to the previous version.  utf-8-sig tolerates a BOM.
f = open(filename, 'r', encoding='utf-8-sig')
header = f.readline()
comp_elems = [c.strip() for c in header.strip().split(',')][1:][::-1]
ncomp = len(comp_elems)
if ncomp not in (2, 3):
    raise ValueError("CSV header must be SFE plus 2 or 3 element columns "
                     "(found %d: %s)" % (ncomp, comp_elems))

sfe_list, comp_list = [], []
for line in f:
    line = line.strip()
    if not line or line.startswith('#'):
        continue
    p = line.split(',')
    sfe_list.append(float(p[0]))
    frac = [float(p[1 + e]) if (1 + e) < len(p) and p[1 + e] != "" else 0.0
            for e in range(ncomp)]
    comp_list.append(frac[::-1])          # reverse to match comp_elems order
f.close()

sfe_vals = np.array(sfe_list, dtype=float)
comp_data = np.array(comp_list, dtype=float)   # columns follow comp_elems
if len(sfe_vals) < ncomp:
    raise ValueError("CSV needs at least the pure-element rows (one per corner)")

#print("USFE Data:")
#print(usfe_data)

# Map lattice element names -> composition index (same order as comp_data).
name_to_idx = {name: e for e, name in enumerate(comp_elems)}
missing = sorted(set(np.unique(elems).tolist()) - set(name_to_idx))
if missing:
    raise ValueError("Lattice has elements not in the CSV header: %s" % missing)
lattice = np.array([name_to_idx[e] for e in elems], dtype=int)

# Build the interpolant from the CSV points:
#   3 components -> embed the composition triangle in 2D, then grid-interpolate
#   2 components -> 1D interpolation of SFE vs the fraction of element index 1
if ncomp == 3:
    grid_x, grid_y = np.mgrid[0:1:2001j, 0:np.sqrt(3)/2:1001j]
    xy = np.zeros((len(comp_data), 2))
    xy[:, 0] = 0.5 * (2 * comp_data[:, 1] + comp_data[:, 2])
    xy[:, 1] = np.sqrt(3) / 2 * comp_data[:, 2]
    usfe_interp = griddata(xy, sfe_vals, (grid_x, grid_y), method='linear')
else:
    order = np.argsort(comp_data[:, 1])
    t_known = comp_data[order, 1]
    sfe_known = sfe_vals[order]

#print (lattice)

print('calculating local USFEs')
#DETERMINE LOCAL COMPOSITIONS
neigh_comp = np.zeros((3,len(lattice),16),dtype=int)
for i in range(N1):
    for j in range(N2):
        for k in range(N3):
            idx = N1*N2*k + N1*j + i

            #### (211) slip plane ####
            #1st NN
            neigh_comp[0,idx,0] = N1*N2*k + N1*((j+1)%N2) + i
            neigh_comp[0,idx,1] = N1*N2*k + N1*((j-1)%N2) + i
            #3rd NN
            neigh_comp[0,idx,2] = N1*N2*((k+1)%N3) + N1*j + i
            neigh_comp[0,idx,3] = N1*N2*((k-1)%N3) + N1*j + i
            neigh_comp[0,idx,4] = N1*N2*k + N1*j + (i+1)%N1
            neigh_comp[0,idx,5] = N1*N2*k + N1*j + (i-1)%N1
            #4th NN
            neigh_comp[0,idx,6] = N1*N2*((k+1)%N3) + N1*((j+1)%N2) + i
            neigh_comp[0,idx,7] = N1*N2*((k+1)%N3) + N1*((j-1)%N2) + i
            neigh_comp[0,idx,8] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + i
            neigh_comp[0,idx,9] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + i
            neigh_comp[0,idx,10] = N1*N2*k + N1*((j+1)%N2) + (i+1)%N1
            neigh_comp[0,idx,11] = N1*N2*k + N1*((j-1)%N2) + (i+1)%N1
            neigh_comp[0,idx,12] = N1*N2*k + N1*((j+1)%N2) + (i-1)%N1
            neigh_comp[0,idx,13] = N1*N2*k + N1*((j-1)%N2) + (i-1)%N1
            #5th NN
            neigh_comp[0,idx,14] = N1*N2*k + N1*((j+2)%N2) + i
            neigh_comp[0,idx,15] = N1*N2*k + N1*((j-2)%N2) + i

            #### (-1-21) slip plane ####
            #1st NN
            neigh_comp[1,idx,0] = N1*N2*k + N1*((j+1)%N2) + i
            neigh_comp[1,idx,1] = N1*N2*k + N1*((j-1)%N2) + i
            #3rd NN
            neigh_comp[1,idx,2] = N1*N2*k + N1*j + (i-1)%N1
            neigh_comp[1,idx,3] = N1*N2*((k-1)%N3) + N1*j + i
            neigh_comp[1,idx,4] = N1*N2*((k-1)%N3) + N1*j + (i+1)%N1
            neigh_comp[1,idx,5] = N1*N2*((k+1)%N3) + N1*j + (i-1)%N1
            #4th NN
            neigh_comp[1,idx,6] = N1*N2*k + N1*((j+1)%N2) + (i-1)%N1
            neigh_comp[1,idx,7] = N1*N2*k + N1*((j-1)%N2) + (i-1)%N1
            neigh_comp[1,idx,8] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + i
            neigh_comp[1,idx,9] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + i
            neigh_comp[1,idx,10] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + (i+1)%N1
            neigh_comp[1,idx,11] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + (i+1)%N1
            neigh_comp[1,idx,12] = N1*N2*((k+1)%N3) + N1*((j+1)%N2) + (i-1)%N1
            neigh_comp[1,idx,13] = N1*N2*((k+1)%N3) + N1*((j-1)%N2) + (i-1)%N1
            #5th NN
            neigh_comp[1,idx,14] = N1*N2*k + N1*((j+2)%N2) + i
            neigh_comp[1,idx,15] = N1*N2*k + N1*((j-2)%N2) + i

            #### (-11-2) ####
            #### slip plane ####
            #1st NN
            neigh_comp[2,idx,0] = N1*N2*k + N1*((j+1)%N2) + i
            neigh_comp[2,idx,1] = N1*N2*k + N1*((j-1)%N2) + i
            #3rd NN
            neigh_comp[2,idx,2] = N1*N2*((k+1)%N3) + N1*j + i
            neigh_comp[2,idx,3] = N1*N2*((k-1)%N3) + N1*j + i
            neigh_comp[2,idx,4] = N1*N2*k + N1*j + (i+1)%N1
            neigh_comp[2,idx,5] = N1*N2*((k-1)%N3) + N1*j + (i+1)%N1
            #4th NN
            neigh_comp[2,idx,6] = N1*N2*((k+1)%N3) + N1*((j+1)%N2) + i
            neigh_comp[2,idx,7] = N1*N2*((k+1)%N3) + N1*((j-1)%N2) + i
            neigh_comp[2,idx,8] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + i
            neigh_comp[2,idx,9] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + i
            neigh_comp[2,idx,10] = N1*N2*k + N1*((j+1)%N2) + (i+1)%N1
            neigh_comp[2,idx,11] = N1*N2*k + N1*((j-1)%N2) + (i+1)%N1
            neigh_comp[2,idx,12] = N1*N2*((k-1)%N3) + N1*((j+1)%N2) + (i+1)%N1
            neigh_comp[2,idx,13] = N1*N2*((k-1)%N3) + N1*((j-1)%N2) + (i+1)%N1
            #5th NN
            neigh_comp[2,idx,14] = N1*N2*k + N1*((j+2)%N2) + i
            neigh_comp[2,idx,15] = N1*N2*k + N1*((j-2)%N2) + i
            
#translate atomic coords to compositions
comp = np.zeros((3,N1*N2*N3,ncomp))
for i in range(N1):
    for j in range(N2):
        for k in range(N3):
            idx = N1*N2*k + N1*j + i
            for l in range(3):
                for n in neigh_comp[l,idx]:
                    comp[l,idx,int(lattice[n])] += 1
                comp[l,idx,int(lattice[idx])] += 1
                comp[l,idx] /= len(neigh_comp[l,idx])+1
print("Calculated local compositions")
#print(f"Dimensions of 'comp': {comp.shape}")

#CONVERT COMPOSITION TO USFE
usfe = np.zeros((3,N1*N2*N3))
if ncomp == 3:
    x = np.round(0.5*(2*comp[:,:,1]+comp[:,:,2])*2000).astype(int)
    y = np.round(comp[:,:,2]*1000).astype(int)
    usfe_near = None                 # nearest-neighbour grid, built only if needed
    n_extrap = 0
    for i in range(N1*N2*N3):
        for j in range(3):
            v = usfe_interp[x[j,i], y[j,i]]
            if v != v:               # NaN: local comp is outside the CSV convex hull
                if usfe_near is None:
                    usfe_near = griddata(xy, sfe_vals, (grid_x, grid_y),
                                         method='nearest')
                v = usfe_near[x[j,i], y[j,i]]
                n_extrap += 1
            usfe[j,i] = v
    if n_extrap:
        print("WARNING: %d/%d sites fell outside the CSV composition range and "
              "used nearest-neighbour extrapolation. Add pure-element corner rows "
              "to the CSV (e.g. sfe,0,0,1 / sfe,0,1,0 / sfe,1,0,0) for exact "
              "interpolation near the corners." % (n_extrap, 3*N1*N2*N3))
else:
    tq = comp[:,:,1]
    usfe = np.interp(tq, t_known, sfe_known)   # SFE vs element[1] fraction (clamps)
    n_clamp = int(np.count_nonzero((tq < t_known[0]) | (tq > t_known[-1])))
    if n_clamp:
        print("WARNING: %d/%d sites were outside the CSV composition range and "
              "were clamped to the nearest endpoint. Add rows nearer the pure "
              "elements for exact values." % (n_clamp, tq.size))
print("Calculated local USFE values")
#print(f"Dimensions of 'usfe': {usfe.shape}")

#WRITE PFDD MPEA USFE INPUT FILE
filename="%s_usfe" % (lattice_file)
f = open(filename,'w')
f.write("%d %d %d %d\n"%(N1,N2,N3,3))
test = 0
for i in range(N1):
    for j in range(N2):
        for k in range(N3):
            for l in range(3):
                f.write("%d %d %d %d %f\n"%(i,j,k,l,usfe[l,N1*N2*k + N1*j + i]))

f.close()


print("Wrote USFE file %s"%(filename))

print(f"Dimensions of 'comp': {comp.shape}")

#Write comp values to a file
output_file = 'monbti_112_1323K_comp_values.txt'
with open(output_file, 'w') as f:
    for i in range(N1 * N2 * N3):
        for j in range(1):
            f.write(f"({i}, {j}): {comp[j, i]}\n")

print(f"Comp values have been written to {output_file}")
