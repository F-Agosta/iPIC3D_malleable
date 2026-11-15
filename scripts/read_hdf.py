import scipy
import numpy as np
import tables
from matplotlib import pyplot as plt

# Directory of the simulation data
directory = 'data/'
# Which quantities to process
data_read = {'Fields': True, 'Currents': True, 'Rho': True, 'Pressure': True}
# Which cycles to process
first_cycle = 0
last_cycle = 6001
step = 60

#ezxp = []

def openProcFile(proc_filename, read_proc0, bounds):
    """Open single proc*.hdf file and return its handler

       Parameters
         proc_filename - name of the input proc*hdf file
         read_proc0 == list_read_proc[0], the bottom processor layer index. See below.
         bounds - oordinates of the edges of the domain saved in the proc file
    """
    try :
        h5proc_file = tables.open_file(proc_filename, mode = "r", title = "Proc_file")
        p_coordinate = h5proc_file.root.topology.cartesian_coord.read()
    except :
        print('Error opening proc file: ', proc_filename)

    p_coordinate[2] = p_coordinate[2] - read_proc0
    bounds[0] = p_coordinate[0]*(shape[0]-1) # x_start -1 to account for the 1 node overlapping
    bounds[1] = bounds[0] + shape[0]           # x_end
    bounds[2] = p_coordinate[1]*(shape[1]-1) # y_start
    bounds[3] = bounds[2] + shape[1]           # y_end
    bounds[4] = p_coordinate[2]*(shape[2]-1) # z_start
    bounds[5] = bounds[4] + shape[2]           # z_end
    return h5proc_file
#end openProcFile

""" Constants """
str_fields = 'fields'
str_rho = 'rho'
str_currents = 'currents'
str_pressure = 'pressure'
str_moments = 'moments'
str_species = 'species_'
current_comp_names = ['Jx', 'Jy', 'Jz']
field_comp_names = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez']
field_comp_count = len(field_comp_names)
field_vec_names = ['B', 'E']
current_vec_name = 'J'
press_names = ['pXX', 'pXY', 'pXZ', 'pYY', 'pYZ', 'pZZ']
press_comp_count = len(press_names)
list_cycle = np.arange(first_cycle, last_cycle+1, step).tolist()
bounds = np.zeros(6) # local domain bounds: x_start, x_end, y_start, ...

setting_filename = directory + "settings.hdf" # Setting file name
h5setting_file = tables.open_file(setting_filename, mode = "r", title = "Setting_file")
Nprocs = h5setting_file.root.topology.Nprocs.read()[0]
XLEN = h5setting_file.root.topology.XLEN.read()[0]
YLEN = h5setting_file.root.topology.YLEN.read()[0]
ZLEN = h5setting_file.root.topology.ZLEN.read()[0]
dx = h5setting_file.root.collective.Dx.read()[0]
dy = h5setting_file.root.collective.Dy.read()[0]
dz = h5setting_file.root.collective.Dz.read()[0]
Nxc = h5setting_file.root.collective.Nxc.read()[0]
Nyc = h5setting_file.root.collective.Nyc.read()[0]
Nzc = h5setting_file.root.collective.Nzc.read()[0]
ns = h5setting_file.root.collective.Ns.read()[0]
b0x = h5setting_file.root.collective.Bx0.read()[0]
h5setting_file.close()
# Shape of each individual subdomain. We assume that each processor has domain of the same size!
shape = (Nxc/XLEN+1, Nyc/YLEN+1, Nzc/ZLEN+1)
species = range(ns)

numproc = 1
rankproc = 0
Zbyproc = ZLEN/numproc # Minimum Nomber of Z you have to read
reste = ZLEN - Zbyproc*numproc
if rankproc < reste:
    Zrange = np.arange(rankproc*(Zbyproc+1),(rankproc+1)*(Zbyproc+1))
else:
    Zrange = np.arange(reste*(Zbyproc+1)+(rankproc-reste)*Zbyproc, reste*(Zbyproc+1)+(rankproc-reste+1)*Zbyproc)

list_read_proc = []
print(Zrange)
for Z in Zrange:
    print(type(Z), type(Nprocs), type(ZLEN))
    list_read_proc.extend(range(np.int32(Z),Nprocs,ZLEN))

list_read_proc = np.array(list_read_proc)
zread_proc = len(list_read_proc)/(YLEN*XLEN) #Number of proc file you want to read in the z direction

Nzlocal = zread_proc*(shape[2] - 1)
print(type(field_comp_count), type(Nxc), type(Nyc), type(Nzlocal))
Fields_local = np.zeros((field_comp_count,Nxc+1,Nyc+1,np.int32(Nzlocal)+1), dtype="float32")
numproc = min(numproc, ZLEN) # To account for specific case where more processors are used here than in the Z direction of the simulation
if len(list_read_proc) > 0:
    numpoints = (Nxc+1)*(Nyc+1)*(Nzlocal+1)
    RangeZlocal = rankproc*Nzc/numproc, (rankproc+1)*Nzc/numproc
    Xcoordinates = np.arange(Nxc+1).flatten().astype('float32')*dx
    Ycoordinates = np.arange(Nyc+1).flatten().astype('float32')*dy
    Zcoordinates = (np.arange(Nzlocal+1).flatten().astype('float32') + list_read_proc[0]*(shape[2]-1))*dz
    '''vtkXcoordinates = vtk.vtkFloatArray()
    vtkXcoordinates.SetNumberOfComponents(1)
    vtkXcoordinates.SetNumberOfTuples(Nxc+1)
    vtkXcoordinates.SetVoidArray(Xcoordinates,Nxc+1,1)
    vtkYcoordinates = vtk.vtkFloatArray()
    vtkYcoordinates.SetNumberOfComponents(1)
    vtkYcoordinates.SetNumberOfTuples(Nyc+1)
    vtkYcoordinates.SetVoidArray(Ycoordinates,Nyc+1,1)
    vtkZcoordinates = vtk.vtkFloatArray()
    vtkZcoordinates.SetNumberOfComponents(1)
    vtkZcoordinates.SetNumberOfTuples(Nzlocal+1)
    vtkZcoordinates.SetVoidArray(Zcoordinates,Nzlocal+1,1)'''
    list_proc_files = []
    for read_proc in list_read_proc:
        list_proc_files.append(directory + "proc" + repr(read_proc) + ".hdf")

    """ 1. Fields: magnetic and electric """
    if data_read['Fields']:
        for cycle in list_cycle:
            print(cycle)
            for proc_filename in list_proc_files:
                h5proc_file = openProcFile(proc_filename, list_read_proc[0], bounds)
                ## Cycle over field components
                for c in range(field_comp_count):
                    for node in h5proc_file.walk_nodes('/' + str_fields + '/'+ field_comp_names[c]):
                        if type(node) == tables.array.Array:
                            num = int(node.name[6:])
                            if num == cycle:
                                ###print(type(c), type(bounds[0]))
                                Fields_local[c, np.int32(bounds[0]):np.int32(bounds[1]), np.int32(bounds[2]):np.int32(bounds[3]), np.int32(bounds[4]):np.int32(bounds[5])] = node.read()
                        #end if
                    #end for
                #end for
                h5proc_file.close()
            #end for
            vector_fields = {}
            for c in range(len(field_vec_names)):
                if field_vec_names[c] == 'E':
                   print(field_vec_names[c])
                   vector_fields[field_vec_names[c]] = Fields_local[3*c:3*(c+1), :, :, :]
                   print(type(vector_fields[field_vec_names[c]]))
                   print(vector_fields[field_vec_names[c]].shape)
                   field = np.transpose(vector_fields[field_vec_names[c]][2,:,:,0]) 
                   ixp = Nxc/4
                   iyp = Nyc/4
                   #ezxp.append(field[3*int(iyp),3*int(ixp)])
                   fig = plt.figure()
                   ax = fig.add_subplot(111)

                   cc = ax.imshow(
                     field,
                     cmap='seismic',
                     ###extent=[xmin, xmax, ymin, ymax],
                     vmin=-0.006,
                     vmax=0.006,
                     interpolation='antialiased',
                     origin='lower'
                     )
                   ax.set_xlabel("x", fontsize=20)
                   ax.set_ylabel("y", fontsize=20)
                   ax.tick_params(labelsize=20)
                   cbar = fig.colorbar(cc)
                   cbar.set_label(field_vec_names[c], fontsize=20)
                   cbar.ax.tick_params(labelsize=20)

                   #fig.tight_layout()

                   plt.savefig(field_vec_names[c] + "_" + (str(cycle)).zfill(6) + ".png", format='png')

#f = open("ezxp.txt", "w")
#for item in ezxp:
#   f.write(str(item) +"\n")
#f.close()

#fig = plt.figure()
#plt.plot(ezxp/b0x**2)
#plt.savefig("Ezxp.png", format='png')
