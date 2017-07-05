import os
import numpy as np
import subprocess
from dplus.CalculationRunner import LocalRunner
from dplus.CalculationInput import GenerateInput
import shutil


def CalcPDBcmd(exe_path,JsonPath, outputpath, name, outfile, errfile):

    name = name[0:name.index('.')]
    process = subprocess.Popen([exe_path, JsonPath], stdout=outfile, stderr= errfile)
    process.wait()
    process.terminate()

    if not os.path.exists(outputpath + '\\Amps'):
        os.makedirs(outputpath + '\\Amps')

    Aname = outputpath + '\\Amps\\' + name + '.amp'

    All_Amps = os.listdir(JsonPath + "\\cache")
    for i in All_Amps:
        shutil.move(JsonPath + "\\cache\\" + i, Aname)

    return Aname

def ReadStateFile(Path):
    statefile = open(Path + "\\args.json",'r+')
    state = statefile.readlines()
    statefile.close()
    return state


def CahngeParams(state, name,q,GridSize, path):
    name = repr(name).replace("'", "")
    pos1 = state[146].index(':')
    pos2 = state[146].index(',')

    state[146] = state[146].replace(state[146][(pos1+1):pos2], '"' + name + '"')

    pos1 = state[183].index(':')
    pos2 = state[183].index(',')
    state[183] = state[183].replace(state[183][(pos1+1):pos2], " " + str(GridSize))

    pos1 = state[185].index(':')
    pos2 = state[185].index(',')
    state[185] = state[185].replace(state[185][(pos1+1):pos2], " " + str(q) )

    newfile = open(path + '\\args.json', 'w+')
    newfile.writelines(state)
    newfile.close()

    return state


def peek(File, length):
    pos = File.tell()
    data = File.read(length)
    File.seek(pos)
    return data

def readAmp(name, qmax):
    f = open(name, "rb+")
    header = []
    offset = 0
    if peek(f, 1).decode('ascii') == '#':
        desc = f.read(2)
        tempdesc = desc.decode('ascii')
        if (tempdesc[1] == '@'):
            offset = np.fromfile(f, dtype=np.uint32, count=1, sep="")
            del_ = f.readline()
        elif (tempdesc[1] != '\n'):
            tmphead = f.readline()
            header.append(tmphead)

    while peek(f, 1).decode('ascii') == '#':
        header.append(f.readline())
    if offset > 0:
        f.seek(offset[0], 0)

    version_r = f.readline().rstrip()
    version = int(version_r.decode('ascii'))
    size_element_r = f.readline().rstrip()
    size_element = int(size_element_r.decode('ascii'))

    if size_element != int(2 * np.dtype(np.float64).itemsize):
        print ("error in file: " + name + "dtype is not float64\n")
        exit(1)

    tmpGridsize_r = f.readline().rstrip()
    tmpGridsize = int(tmpGridsize_r.decode('ascii'))

    tmpExtras_r = f.readline().rstrip()
    tmpExtras = int(tmpExtras_r.decode('ascii'))
    Gridsize = (tmpGridsize - tmpExtras) * 2

    # Total size of the grid
    thetaDivisions = 3
    phiDivisions = 6

    stepSize = qmax / float(Gridsize / 2)
    actualGridSize = Gridsize / 2 + tmpExtras

    i = actualGridSize
    totalsz = int((phiDivisions * i * (i + 1) * (3 + thetaDivisions + 2 * thetaDivisions * i)) / 6)
    totalsz = totalsz + 1
    totalsz = totalsz * 2

    step_size = np.fromfile(f, dtype=np.float64, count=1, sep="")

    Amplitude = np.fromfile(f, dtype=np.float64, count=totalsz, sep="")

    f.close()
    pos = 0
    header_List = []
    header_List.append(desc)
    pos = pos + len(desc)

    header_List.append(offset[0].tobytes())
    pos = pos + len(offset[0].tobytes())
    header_List.append(del_)
    pos = pos + len(del_)

    for i in header:
        header_List.append(i)
        pos = pos + len(i)
        header_List.append(del_)
        header_List.append(del_)
        pos = pos + 2 * len(del_)
    pos = np.int32(pos)
    if pos != offset[0]:
        header_List[1] = pos.tobytes()

    header_List.append(version_r + b"\n")
    header_List.append(size_element_r + b"\n")
    header_List.append(tmpGridsize_r + b"\n")
    header_List.append(tmpExtras_r + b"\n")
    header_List.append(step_size.tobytes())

    return Amplitude, header_List

def writeAmp(SubAveragedAmplitude, HeaderList_Solvent, filename):
    outputFile = open(filename, 'wb')

    for i in HeaderList_Solvent:
        outputFile.write(i)

    SubAveragedAmplitude = np.float64(SubAveragedAmplitude)
    SubAveragedAmplitude.tofile(outputFile, sep="")

    outputFile.close()

    return  filename


def averageAmps(AmpList, qmax, _filename):
    count = 0
    path = os.path.dirname(AmpList[0])
    for name in AmpList:
        if count == 0:
            Amplitude, header = readAmp(name, qmax)
        else:
            tempAmp, tempHeader = readAmp(name, qmax)
            Amplitude = Amplitude + tempAmp
        count = count + 1

    Ave_Amplitude = Amplitude / np.float64(count)

    # Exporting the averaged Amplitude file
    outputFile = path + "\\" + _filename

    return Ave_Amplitude , writeAmp(Ave_Amplitude,header,outputFile), header




def averageIntensities(ListOfFiles, outputfilename):

    data = []
    ave_I = []
    filesnames = []

    count = 1
    for name in ListOfFiles:
        f = open(name, 'r+')
        l = f.readlines()
        f.close()

        path, filename = os.path.split(os.path.abspath(name))
        filesnames.append(filename)

        fparams = fileParams(l)

        if count == 1:
            header = l[0:fparams[0, 1]]
        Intensity = np.loadtxt(l[0:len(l) - fparams[0, 2]], dtype='float', skiprows=fparams[0, 1], usecols={1})
        data.append(Intensity[0:(len(Intensity))])
        Q = np.loadtxt(l[0:len(l) - fparams[0, 2]], dtype='float', skiprows=fparams[0, 1], usecols={0})
        count = count + 1

    for v in range(0, len(Intensity)):
        temp_I = 0
        for n in range(0, len(ListOfFiles)):
            temp_I = temp_I + data[n][v] / len(ListOfFiles)
        ave_I.append(temp_I)

    newfile = open(outputfilename, 'w+')

    for filename in filesnames:
        newfile.write('#' + filename + '\n')

    for line in range(0, len(header)):
        newfile.write(str(header[line]))
    for i in range(0, len(Intensity)):
        newfile.write(str(Q[i]) + '\t' + str(ave_I[i]) + '\n')
    newfile.close()

    return Q, ave_I


def creatTree(_StateFile, _qMax, _grid):

    input = GenerateInput.load_from_state(_StateFile)
    input.state.DomainPreferences.grid_size = _grid
    input.state.DomainPreferences.q_max = _qMax

    input._create_x_vector()

    return input

def WriteFile(result, outputfile):
    new_file = open(outputfile,'w+')
    for key, value in result.headers.items():
        new_file.write(value)

    new_file.write("q [nm^-1]" + "\t" + "Intensity [e^2]" + "\n")
    for key, value in result.graph.items():
        new_file.write(str(key) + "\t" + str(value) + "\n")

    new_file.close()


def CalcPDB(Tree,Tree2, api, _path, _name):

    _name = _name[0:_name.index('.')]

    pathA = _path + '\\Amplitudes'
    pathI = _path + '\\Intensities'

    pathA = repr(pathA).replace("'", "")
    pathI = repr(pathI).replace("'", "")
    if not os.path.exists(pathA):
        os.makedirs(pathA)
    if not os.path.exists(pathI):
        os.makedirs(pathI)

    outputfile = pathI + "\\" + _name + '.out'

    Tree.state.Domain.Children[0].Children[0].use_grid = True
    Tree2.state.Domain.Children[0].Children[0].use_grid = True

    result = api.generate(Tree)
    WriteFile(result, outputfile)

    temp_Amp_dir = api._session_directory + "\\cache"

    All_Amps = os.listdir(temp_Amp_dir)
    for i in All_Amps:
        print('Moving file (1):' , i , ' to ', _name)
        shutil.move(temp_Amp_dir + "\\" + i, pathA + "\\" + _name + '.amp')
    print(_name)

    Tree2.state.Domain.Children[0].Children[0].filename = repr(pathA + "\\" + _name + '.amp').replace("'", "")
    result2 = api.generate(Tree2)

    All_Amps = os.listdir(temp_Amp_dir)
    for i in All_Amps:
        os.remove(temp_Amp_dir + "\\" + i)

    WriteFile(result2, outputfile)

    return outputfile, pathA + "\\" + _name + '.amp'


def CalculateIntensityFromAmplitud(SubAveragedAmplitudeName, AmpTree, api):
    _path, _name = os.path.split(SubAveragedAmplitudeName)
    _name = _name[0:_name.index('.')]
    _path = repr(_path).replace("'", "")

    outputfile = _name + '.out'

    AmpTree.state.Domain.Children[0].Children[0].use_grid = True
    AmpTree.state.Domain.Children[0].Children[0].filename = repr(SubAveragedAmplitudeName).replace("'", "")

    result = api.generate(AmpTree)
    WriteFile(result, outputfile)

    Intensity = []
    Qvector = []

    for key, value in result.graph.items():
        Qvector.append(key)
        Intensity.append(value)

    return Qvector, Intensity



def ExpotrIntensity(signal, Q, filename):
    # TODO construct a header stream

    FinalFile = open(filename,'w+')
    FinalFile.write("q [nm^-1]" + "\t" + "Intensity [A.U]" + "\n")

    for i in range(0,len(signal)):
        FinalFile.write(str(Q[i]) + "\t" + str(signal[i]) + "\n")
    FinalFile.close()


def averageSumAmps(AmpListProtein, AmpListSolvent, qMax):
    count = 0
    for i in range (0,len(AmpListProtein)):
        Amp_Protein, headerProtein = readAmp(AmpListProtein[i], qMax)
        Amp_Solvent, headerSolvent = readAmp(AmpListSolvent[i], qMax)
        if i == 0:
            SubAmplitude = Amp_Protein - Amp_Solvent
        else:
            SubAmplitude = SubAmplitude + (Amp_Protein - Amp_Solvent)
        count = count + 1
    SubAmplitude = SubAmplitude / np.float64(count)

    return SubAmplitude, headerSolvent

def fileParams(l):
    """
    returns the lines in the file that contain
    the relevant dada (the q and I without the comments)
    """
    k = 0
    n = 0
    for i in range(0, len(l)):
        if not is_number(l[i][0:4]) or l[i].find('nan') >= 0:
            k = k + 1
        else:
            break
    j = len(l) - 1
    while j > k and (l[j].find('nan') >= 0 or is_number(l[j][0:5].strip("'b")) == False):
        n = n + 1
        j = j - 1
    return np.matrix([len(l), k, n])

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def DeleteFiles(delete,ProteinDir,SolventDir):
        if delete:
            print("Deleting all files...")
            ProteinListDir = os.listdir(ProteinDir)
            SolventListDir = os.listdir(SolventDir)
            for f in ProteinListDir:
                if f.find(".") == -1:
                    shutil.rmtree(ProteinDir + "\\" + f,ignore_errors=True)
            for f in SolventListDir:
                if f.find(".") == -1:
                    shutil.rmtree(SolventDir + "\\" + f,ignore_errors=True)


# def averageAmps(AmpList, qmax):
#     count = 0
#     path = os.path.dirname(AmpList[0])
#     for name in AmpList:
#         f = open(name, "rb+")
#         header = []
#         offset = 0
#         if peek(f, 1).decode('ascii') == '#':
#             desc = f.read(2)
#             tempdesc = desc.decode('ascii')
#             if (tempdesc[1] == '@'):
#                 offset = np.fromfile(f, dtype=np.uint32, count=1, sep="")
#                 del_ = f.readline()
#             elif (tempdesc[1] != '\n'):
#                 tmphead = f.readline()
#                 header.append(tmphead)
#
#
#
#         while peek(f, 1).decode('ascii') == '#':
#             header.append(f.readline())
#         if offset > 0:
#             f.seek(offset[0], 0)
#
#         version_r = f.readline().rstrip()
#         version = int(version_r.decode('ascii'))
#         size_element_r = f.readline().rstrip()
#         size_element = int(size_element_r.decode('ascii'))
#
#         if size_element != int(2 * np.dtype(np.float64).itemsize):
#             print ("error in file: " + name + "dtype is not float64\n")
#             exit(1)
#
#         tmpGridsize_r = f.readline().rstrip()
#         tmpGridsize = int(tmpGridsize_r.decode('ascii'))
#
#         tmpExtras_r = f.readline().rstrip()
#         tmpExtras = int(tmpExtras_r.decode('ascii'))
#         Gridsize = (tmpGridsize - tmpExtras) * 2
#
#         # Total size of the grid
#         thetaDivisions = 3
#         phiDivisions = 6
#
#         stepSize = qmax / float(Gridsize / 2)
#         actualGridSize = Gridsize / 2 + tmpExtras
#
#         i = actualGridSize
#         totalsz = int((phiDivisions * i * (i + 1) * (3 + thetaDivisions + 2 * thetaDivisions * i)) / 6)
#         totalsz = totalsz + 1
#         totalsz = totalsz * 2
#
#         step_size = np.fromfile(f, dtype=np.float64, count=1, sep="")
#
#         if count == 0:
#             Amplitude = np.fromfile(f, dtype=np.float64, count=totalsz, sep="")
#         else:
#             tempAmp = np.fromfile(f, dtype=np.float64, count=totalsz, sep="")
#             Amplitude = Amplitude + tempAmp
#         count = count + 1
#         f.close()
#
#     Ave_Amplitude = Amplitude / np.float64(count)
#
#     # Exporting the averaged Amplitude file
#     outputFile = open(path + "\\" + "Averaged_Amp.amp", 'wb')
#
#     ##TODO - What sould the heaser include
#     header_List = []
#     outputFile.write(desc)
#     header_List.append(desc)
#
#     offset[0].tofile(outputFile, sep="")
#     header_List.append(offset[0].tobytes())
#
#     outputFile.write(del_)
#     header_List.append(del_)
#
#     for i in header:
#         outputFile.write(i)
#         header_List.append(i)
#
#     outputFile.write(del_)
#     outputFile.write(del_)
#     header_List.append(del_)
#     header_List.append(del_)
#
#     pos = np.uint32(outputFile.tell())
#     if pos != offset[0]:
#         outputFile.seek(0, 0)
#         outputFile.write(desc)
#         pos.tofile(outputFile, sep="")
#         header_List[1] = pos.tobytes()
#         outputFile.write(del_)
#
#     outputFile.seek(pos, 0)
#     outputFile.write(version_r + b"\n")
#     header_List.append(version_r + b"\n")
#
#     outputFile.write(size_element_r + b"\n")
#     header_List.append(size_element_r + b"\n")
#
#     outputFile.write(tmpGridsize_r + b"\n")
#     header_List.append(tmpGridsize_r + b"\n")
#
#     outputFile.write(tmpExtras_r + b"\n")
#     header_List.append(tmpExtras_r + b"\n")
#
#     step_size.tofile(outputFile, sep="")
#     header_List.append(step_size.tobytes())
#
#     Ave_Amplitude = np.float64(Ave_Amplitude)
#     Ave_Amplitude.tofile(outputFile, sep="")
#
#     outputFile.close()
#
#     return Ave_Amplitude , path + "\\" + "Averaged_Amp.amp", header_List
#