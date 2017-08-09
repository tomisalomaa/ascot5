import sys
sys.path.append('../ui')

from pylab import *
import ui_B_2D
import ui_B_3D

def read_magn_bkg(fn):
    str = dict()
    fh = open(fn)
    tmp = map(float,fh.readline().split())
    str['phi0']       = tmp[0]
    str['nSector']    = int(tmp[1])
    str['nPhi']       = int(tmp[2])
    str['nCoil']      = int(tmp[3])
    str['zeroAtCoil'] = int(tmp[4])    

    r1,r2,nr = map(float,fh.readline().split())
    nr = int(nr)
    z1,z2,nz = map(float,fh.readline().split())
    nz = int(nz)

    str['r'] = linspace(r1,r2,nr)
    str['z'] = linspace(z1,z2,nz)
    if(str['nPhi'] > 1):
        dphi = 360.0 / ( str['nSector'] * str['nPhi'] )
        str['phi'] = ( linspace( str['phi0'] + dphi * 0.5,
                                 str['phi0'] + 360.0 / str['nSector'] - dphi * 0.5,
                                 str['nPhi'] ) )

    str['phimap_tor'] = array(map(int,fh.readline().split()))
    str['phimap_pol'] = array(map(int,fh.readline().split()))

    #data = loadtxt(fh).flatten()
    data = array(fh.read().split(), dtype=float).flatten()

    sz2d = nr*nz
    sz3d = nr*nz*str['nPhi']
    str['psi']  = data[:sz2d].reshape(nz,nr)
    print shape(data),nz,str['nPhi'],nr,sz2d,sz3d
    str['br']   = data[sz2d+0*sz3d:sz2d+1*sz3d].reshape(nz,str['nPhi'],nr).squeeze()
    str['bphi'] = data[sz2d+1*sz3d:sz2d+2*sz3d].reshape(nz,str['nPhi'],nr).squeeze()
    str['bz']   = data[sz2d+2*sz3d:sz2d+3*sz3d].reshape(nz,str['nPhi'],nr).squeeze()

    fh.close()
    
    read_magn_header('input.magn_header',str)
    return str

def read_magn_header(fn,str):
    fh = open(fn)

    # first four lines contain nothing interesting
    fh.readline()
    fh.readline()
    fh.readline()
    fh.readline()

    # Next three lines contain axis psi, R, and z values
    tmp = map(float,fh.readline().split())
    str['psi0']       = tmp[0]
    str['psi1']       = tmp[1]

    tmp = map(float,fh.readline().split())
    str['axis_r']       = tmp[0]

    tmp = map(float,fh.readline().split())
    str['axis_z']       = tmp[0]

    fh.close()
    return str

def read_bfield_stellarator(f):
    str = dict()
    str['bphi'] = f['bfield/stellarator/']
    return str
    
def write_magn_bkg(f, m):
    if(m['nPhi'] > 1):
        ui_B_3D.write_hdf5('ascot.h5',
                           np.array([m['r'][0], m['r'][-1]]),
                           np.array([m['z'][0], m['z'][-1]]), m['nPhi'],
                           m['psi']/(2*np.pi), m['br'], m['bphi'], m['bz'],
                           np.array([m['axis_r'], m['axis_z']]), np.array([m['psi0'], m['psi1']]))
    else:
        ui_B_2D.write_hdf5('ascot.h5',
                           np.array([m['r'][0], m['r'][-1]]),
                           np.array([m['z'][0], m['z'][-1]]),
                           m['psi']/(2*np.pi), m['br']*0, m['bphi'], m['bz']*0,
                           np.array([m['axis_r'], m['axis_z']]), np.array([m['psi0'], m['psi1']]))
    return
