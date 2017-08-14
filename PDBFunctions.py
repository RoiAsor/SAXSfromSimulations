import os, math
import numpy as np
from string import digits

def SimToDplusFormat(name):
    outputPath = os.path.dirname(name) + '\\D+Format\\'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    PDBRead = open(name, 'r+')
    PDB = PDBRead.readlines()
    PDBRead.close()

    outputfilename = outputPath + 'dplus_format_' + os.path.basename(name)
    WriteF = open(outputfilename, 'w+')

    for line in PDB:
        if (line[13:15] != 'MW'):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                type = Atomtype(line[12:16])
                res = Atomtype(line[17:20])
                if type == 'NA' and res == 'NA':
                    line  = line[:76] + 'NA1+' + '\n'
                elif type == 'CL' and res == 'CL':
                    line = line[:76] + 'CL1-' + '\n'

                else:
                    if type.startswith("H"):
                        if line.startswith("ATOM"):
                            type = " H"
                    elif type.startswith("C"):
                        if line.startswith("ATOM"):
                            type = " C"
                    elif type.startswith("O"):
                        if line.startswith("ATOM"):
                            type = " O"
                    elif type.startswith("N"):
                        if line.startswith("ATOM") and res != "SOL":
                            type = " N"
                    elif type.startswith("S"):
                        if line.startswith("ATOM") and res == "CYS":
                            type = " S"
                    if len(type) == 1:
                        line = line[:76] + ' ' + type + line[78:]
                    else:
                        line = line[:76] + type + line[78:]
            WriteF.write(line)
    WriteF.close()

    return outputfilename

def Atomtype(string):
    new_string = string
    remove_digits = str.maketrans('', '', digits)
    new_string = string.translate(remove_digits)
    new_string = new_string.strip()
    #print("old string = ", string , "new string = ", new_string)
    return new_string

def PDBtoCOM(PDBname):

    path = os.path.dirname(PDBname)
    basename = os.path.basename(PDBname)
    outputPath = path + '\\D+Format_GC\\'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    PDBRead = open(PDBname, 'r+')
    PDB = PDBRead.readlines()
    PDBRead.close()

    xcenter = 0.0
    zcenter = 0.0
    ycenter = 0.0
    totMass = 0.0
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if line[77] == 'C':
                Mass = 12.0
            elif line[77] == 'N':
                Mass = 14.0
            elif line[77] == 'O':
                Mass = 16.0
            elif line[77] == 'H':
                Mass = 1.0
            elif line[77] == 'S':
                Mass = 32.0
            else:
                Mass = 1.0
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            xcenter = xcenter + x * Mass
            ycenter = ycenter + y * Mass
            zcenter = zcenter + z * Mass
            totMass = totMass + Mass


    xcenter = xcenter / totMass
    ycenter = ycenter / totMass
    zcenter = zcenter / totMass

    outputfile = outputPath + 'dplus_format_COM_' + basename
    WriteF = open(outputfile, 'w+')
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            newX = x - xcenter
            newY = y - ycenter
            newZ = z - zcenter
            newline = coord_replace(line, line[30:38], line[38:46], line[46:54], newX, newY, newZ)
            WriteF.write(newline)
        else:
            WriteF.write(line)
    WriteF.close()
    return outputfile


def ProteintoCOM(PDBname):
    xyz = []
    PDBRead = open(PDBname, 'r+')
    PDB = PDBRead.readlines()
    PDBRead.close()
    i = 0  ## line index
    k = 0  ## atoms index
    xcenter = 0.0
    zcenter = 0.0
    ycenter = 0.0
    totMass = 0.0
    count = 0
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if line[17:20] != 'SOL' and line[18:20] != 'NA' and line[18:20] != 'CL':
                if line[77] == 'C':
                    Mass = 12.0
                elif line[77] == 'N':
                    Mass = 14.0
                elif line[77] == 'O':
                    Mass = 16.0
                elif line[77] == 'H':
                    Mass = 1.0
                elif line[77] == 'S':
                    Mass = 32.0
                else:
                    Mass = 1.0
                    count = count + 1

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                xcenter = xcenter + x*Mass
                ycenter = ycenter + y*Mass
                zcenter = zcenter + z*Mass
                totMass = totMass + Mass
                if line[13:15] == 'CA':
                    xyz.append([x,y,z])
                k = k + 1
        i = i + 1

    xcenter = xcenter / totMass
    ycenter = ycenter / totMass
    zcenter = zcenter / totMass
    out_center = []
    for i in range(0, len(xyz)):
        xyz[i][0] = xyz[i][0] - xcenter
        xyz[i][1] = xyz[i][1] - ycenter
        xyz[i][2] = xyz[i][2] - zcenter
    out_center.append([xcenter, ycenter, zcenter])

    return xyz, out_center


def SolventtoCOM(PDBname):
    xyz = []
    PDBRead = open(PDBname, 'r+')
    PDB = PDBRead.readlines()
    PDBRead.close()
    xcenter = 0.0
    zcenter = 0.0
    ycenter = 0.0
    totMass = 0.0
    count = 0
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):

            if line[77] == 'C':
                Mass = 12.0
            elif line[77] == 'N':
                Mass = 14.0
            elif line[77] == 'O':
                Mass = 16.0
            elif line[77] == 'H':
                Mass = 1.0
            elif line[77] == 'S':
                Mass = 32.0
            else:
                Mass = 1.0
                count = count + 1

            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            xcenter = xcenter + x*Mass
            ycenter = ycenter + y*Mass
            zcenter = zcenter + z*Mass
            totMass = totMass + Mass

            xyz.append([x,y,z])

    xcenter = xcenter / totMass
    ycenter = ycenter / totMass
    zcenter = zcenter / totMass
    out_center = []
    for i in range(0, len(xyz)):
        xyz[i][0] = xyz[i][0] - xcenter
        xyz[i][1] = xyz[i][1] - ycenter
        xyz[i][2] = xyz[i][2] - zcenter
    out_center.append([xcenter, ycenter, zcenter])

    return xyz, out_center



def CutSphericalEnvelope(name,COM, Rmax):

    PDBRead = open(name, 'r+')
    PDB = PDBRead.readlines()
    PDBRead.close()


    path = os.path.dirname(name)
    basename = os.path.basename(name)
    outputPath = path + '\\SphericalEnv\\'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)

    outputfilename = outputPath + 'Env_' +basename
    NewPDB = open(outputfilename, 'w+')

    error = False
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                RfromCOM = np.sqrt((x - COM[0][0])**2 + (y - COM[0][1])**2 + (z - COM[0][2])**2)
                if Rmax >= RfromCOM:
                    newline = coord_replace(line, line[30:38], line[38:46], line[46:54], x-COM[0][0], y - COM[0][1], z - COM[0][2])
                    NewPDB.write(newline)
                if is_protein(line) and RfromCOM > Rmax-12:
                    error = True

        else:
            NewPDB.write(line)
    NewPDB.close()
    return outputfilename, error

def is_protein(line):
    Protein = False
    if line[17:20] != 'SOL' and line[18:20] != 'NA' and line[18:20] != 'CL':
        Protein = True
    return Protein


def Euler_angles(xyz, ref_xyz):
    ref_coord = np.copy(ref_xyz)
    coord = np.copy(xyz)
    Matrix, res, rank, s = np.linalg.lstsq(ref_coord, coord)
    Engles = Euler_angles_abg(np.transpose(Matrix))
    return Engles[0], Engles[1], Engles[2]


def rotation_axisN(ux, uy, uz, theta):
    matrix = np.matrix(((math.cos(theta) + ux ** 2 * (1 - math.cos(theta)),
                         ux * uy * (1 - math.cos(theta)) - uz * math.sin(theta),
                         ux * uz * (1 - math.cos(theta)) + uy * math.sin(theta)),
                        (uy * ux * (1 - math.cos(theta)) + uz * math.sin(theta),
                         math.cos(theta) + uy ** 2 * (1 - math.cos(theta)),
                         uy * uz * (1 - math.cos(theta)) - ux * math.sin(theta)),
                        (uz * ux * (1 - math.cos(theta)) - uy * math.sin(theta),
                         uz * uy * (1 - math.cos(theta)) + ux * math.sin(theta),
                         math.cos(theta) + uz ** 2 * (1 - math.cos(theta)))))
    return matrix


def Euler_angles_abg(Matrix):
    epsilon = 10 ** (-5)
    Rmatrix = np.copy(Matrix)

    if (1.0 - np.abs(Rmatrix[0][2])) > epsilon:

        # finding theta:

        thet1 = np.arcsin(Rmatrix[0][2])
        thet2 = np.pi - thet1

        # finding psi:

        psi1 = np.arctan2(-Rmatrix[1][2] / np.cos(thet1), Rmatrix[2][2] / np.cos(thet1))
        psi2 = np.arctan2(-Rmatrix[1][2] / np.cos(thet2), Rmatrix[2][2] / np.cos(thet2))

        # finding phi:

        phi1 = np.arctan2(-Rmatrix[0][1] / np.cos(thet1), Rmatrix[0][0] / np.cos(thet1))
        phi2 = np.arctan2(-Rmatrix[0][1] / np.cos(thet2), Rmatrix[0][0] / np.cos(thet2))

        Eangles1 = np.asarray([psi1, thet1, phi1])
        Eangles2 = np.asarray([psi2, thet2, phi2])

        Eangles = Eangles1
    else:
        phi = 0.0
        if (1.0 - Rmatrix[0][2]) < epsilon:
            thet = np.pi / 2
            psi = phi + np.arctan2(Rmatrix[1][0], -Rmatrix[2][0])
        else:
            thet = -1 * np.pi / 2
            psi = - phi + np.arctan2(-Rmatrix[1][0], Rmatrix[2][0])
        Eangles2 = np.asarray([psi, thet, phi])
        Eangles = Eangles2

    for l in range(0, 3):
        if Eangles[l] < 0:
            Eangles[l] = Eangles[l] + 2 * np.pi

    return Eangles


def Euler_Matrix(alpha, beta, gamma):
    c1 = np.cos(alpha)
    c2 = np.cos(beta)
    c3 = np.cos(gamma)
    s1 = np.sin(alpha)
    s2 = np.sin(beta)
    s3 = np.sin(gamma)
    matrix = np.matrix(((c2 * c3, -c2 * s3, s2),
                        (c1 * s3 + c3 * s1 * s2, c1 * c3 - s1 * s2 * s3, -c2 * s1),
                        (s1 * s3 - c1 * c3 * s2, c3 * s1 + c1 * s2 * s3, c1 * c2)))
    return matrix

def read_pdb_Protein_xyz(pdb_name):
    """
read xyz from pdb
return:
[[x1 y1 z1]
[x2 y2 z2]
[.. .. ..]
[xn yn zn]]
   """
    xyz = []
    pdb_file = open(pdb_name, 'r+')
    count = 0
    xcenter = 0
    ycenter = 0
    zcenter = 0
    totMass = 0
    FullFile = []
    out_center = []

    for line in pdb_file:
        FullFile.append(line)
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if line[17:20] != 'SOL' and line[18:20] != 'NA' and line[18:20] != 'CL':
                if line[77] == 'C':
                    Mass = 12.0
                elif line[77] == 'N':
                    Mass = 14.0
                elif line[77] == 'O':
                    Mass = 16.0
                elif line[77] == 'H':
                    Mass = 1.0
                elif line[77] == 'S':
                    Mass = 32.0
                else:
                    Mass = 1.0
                    count = count + 1

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                xcenter = xcenter + x * Mass
                ycenter = ycenter + y * Mass
                zcenter = zcenter + z * Mass
                totMass = totMass + Mass
                if line[13:15] == 'CA':
                    xyz.append([x, y, z])

    pdb_file.close()
    xcenter = xcenter / totMass
    ycenter = ycenter / totMass
    zcenter = zcenter / totMass
    for i in range(0, len(xyz)):
        xyz[i][0] = xyz[i][0] - xcenter
        xyz[i][1] = xyz[i][1] - ycenter
        xyz[i][2] = xyz[i][2] - zcenter
    out_center.append([xcenter, ycenter, zcenter])
    return FullFile, xyz, out_center


def ExportRotatedPDB(filename, ref_xyz):
    path = os.path.dirname(filename)
    basename = os.path.basename(filename)
    outputPath = path + '\\Rotated\\'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    outputfile = outputPath + 'Rotated_' + basename
    PDBfile, xyz, COM = read_pdb_Protein_xyz(filename)
    alpha, beta, gamma = Euler_angles(xyz, ref_xyz)
    RotatePDB(PDBfile, outputfile, alpha, beta, gamma, COM)
    return outputfile, alpha, beta, gamma

def RotatePDB(PDB, filename, alpha, beta, gamma,center):
    R = Euler_Matrix(alpha, beta, gamma)
    PDB = np.copy(PDB)
    WriteF = open(filename, 'w+')

    vec = np.matrix(np.zeros((3,1)))
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            vec[0,0] = float(line[30:38]) - center[0][0]
            vec[1, 0] = float(line[38:46]) - center[0][1]
            vec[2, 0] = float(line[46:54]) - center[0][2]
            new_vec = np.transpose(R) * vec
            newX = new_vec[0,0] +  (center[0][0])
            newY = new_vec[1,0] + (center[0][1])
            newZ = new_vec[2,0] + (center[0][2])
            newline = coord_replace(line, line[30:38], line[38:46], line[46:54], newX, newY, newZ)
            WriteF.write(newline)
        else:
            WriteF.write(line)
    WriteF.close()


def RotateSolvent(filename, alpha, beta, gamma):

    path = os.path.dirname(filename)
    basename = os.path.basename(filename)
    outputPath = path + '\\Rotated\\'
    if not os.path.exists(outputPath):
        os.makedirs(outputPath)
    outputfile = outputPath + 'Rotated_' + basename

    R = Euler_Matrix(alpha, beta, gamma)

    PDB  = open(filename, 'r+')
    WriteF = open(outputfile, 'w+')

    vec = np.matrix(np.zeros((3,1)))
    for line in PDB:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            vec[0,0] = float(line[30:38])
            vec[1, 0] = float(line[38:46])
            vec[2, 0] = float(line[46:54])
            new_vec = np.transpose(R) * vec
            newX = new_vec[0,0]
            newY = new_vec[1,0]
            newZ = new_vec[2,0]
            newline = coord_replace(line, line[30:38], line[38:46], line[46:54], newX, newY, newZ)
            WriteF.write(newline)
        else:
            WriteF.write(line)
    PDB.close()
    WriteF.close()
    return outputfile


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False



def stringcut(s,size):
    s = "%f" % s
    s = s[0:size+1]
    index_s = s.index('.')
    num_s = float(s)
    num_s = round(num_s, int(size-index_s-1))
    num_s = '%f' % num_s
    num_s = num_s[0:size]
    if num_s.index('.'):
        while len(num_s) < size:
            num_s = num_s + '0'
    else:
        while len(num_s) < size:
            num_s = num_s + ' '
    return num_s


def coord_replace(a, old_x, old_y, old_z,x_, y_, z_):
    x = stringcut(x_,8)
    y = stringcut(y_,8)
    z = stringcut(z_, 8)
    string = x + y + z
    d = a.replace(a[30:54], string)
    return d