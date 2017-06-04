import Functions as Funcs
import PDBFunctions as PDBFuncs
import os
import numpy as np
import argparse
from dplus.CalculationRunner import LocalRunner
from dplus.CalculationInput import GenerateInput

def main():
    parser = argparse.ArgumentParser(description='Calculates the x-ray scattering intensity from simulation results.'
                                                 'Needs the simulation box of the solvent and the simulation box of the solvated protein.'
                                                 'The program will create sub folders with additional .pdb files as default.'
                                                 'The Calculation follows Eq. 10 in Chen P. and Hub J. S. "Validating Solution Ensembles'
                                                 'from Molecular Dynamics Simulation by Wide-Angle X-ray Scattering Data"'
                                                 '(2014) Biophysical Journal, 107, 435__447'
                                     )
    parser.add_argument(dest='ProteinDir', metavar='<Protein Dir>', type=str,
                        help='The folder that contains the .pdb files of the solvated proteins')
    parser.add_argument(dest='SolventDir', metavar='<Solvent Dir>', type=str,
                        help='The folder that contains the .pdb files of the solvent only')
    parser.add_argument('-d', dest='delete', action='store_const',
                        const=True, default=False,
                        help='Add -d if you wish to delete all the additional .pdb files.')
    parser.add_argument('-R', dest='Rotate', action='store_const',
                        const=True, default=False,
                        help='Add -R if you wish to align the proteins.')
    parser.add_argument('-E', dest='Envelope', action='store_const',
                        const=True, default=False,
                        help='Add -E if you wish to cut a spherical envelope around the protein. If yes you must provide Rmax.')
    parser.add_argument('--RMax', metavar='<RMax>', type=float,
                        help='The max radius of the spherical envelope.')
    parser.add_argument('-m', dest='Modify', action='store_const',
                        const=True, default=False,
                        help='Add -m if the PDB files need to be modified before the calculations. '
                             '(move to geometric center, add formal charges, remove lone pairs of the water force field)')
    parser.add_argument(dest='qMax', metavar='<qMax>', type=float,
                        help='The max q value you wish to calculate, right now the maximum q range is 10')
    parser.add_argument(dest='Grid', metavar='<Grid size>', type=float,
                        help='The amplitude grid size.')
    args = parser.parse_args()

    # Chacking the arguments for the program

    if os.path.exists(args.ProteinDir) != True:
        print("The protein path entered is invalid!")
        exit()

    if os.path.exists(args.SolventDir) != True:
        print("The Solvent path entered is invalid!")
        exit()

    if len(os.listdir(args.ProteinDir)) == 0:
        print("The protein path entered is empty!")
        exit()

    if len(os.listdir(args.SolventDir)) == 0:
        print("The Solvent path entered is empty!")
        exit()

    if args.RMax == None and args.Envelope == True:
        print("Must enter the maximum radius of the envelope")
        exit()

    if args.RMax and args.Envelope == False:
        print("For envelope cut add -E")
        exit()
    # -------------------------------------------------------------------------------------------

    # creating a list of files for the protein and the solvent PDB files
    print("Reading protein files...")
    ListPDBFiles = []
    for file in os.listdir(args.ProteinDir):
        if file.endswith(".pdb"):
            ListPDBFiles.append(args.ProteinDir + '\\' + file)

    print("Reading solvent files...")
    ListSolventFiles = []
    for file in os.listdir(args.SolventDir):
        if file.endswith(".pdb"):
            ListSolventFiles.append(args.SolventDir + '\\' + file)

    if args.Modify:
        # Converting the PDB files to D+ format and savind it
        print("Converting to D plus format ...")
        new_filenamesProtein = []
        new_filenamesSolvent = []

        for name in ListPDBFiles:
            new_filenamesProtein.append(PDBFuncs.SimToDplusFormat(name))
        for name in ListSolventFiles:
            new_filenamesSolvent.append(PDBFuncs.SimToDplusFormat(name))

        # moving the geometric center of the PDB to (0,0,0)
        print("Moving to geometric center...")
        new_filenamesProtein_GC = []
        new_filenamesSolvent_GC = []
        for name in new_filenamesProtein:
            new_filenamesProtein_GC.append(PDBFuncs.PDBtoCOM(name))
        for name in new_filenamesSolvent:
            new_filenamesSolvent_GC.append(PDBFuncs.PDBtoCOM(name))

        FinalProteinList = np.copy(new_filenamesProtein_GC)
        FinalSolventList = np.copy(new_filenamesSolvent_GC)
    else:
        FinalProteinList = np.copy(ListPDBFiles)
        FinalSolventList = np.copy(ListSolventFiles)

    # optional - if the frames need to be rotated
    if (args.Envelope == True):
        print("Cutting spherical envelope...")
        new_filenamesProtein_env = []
        new_filenamesSolvent_env = []
        errorList = []
        for name in FinalProteinList:
            xyz, out_center = PDBFuncs.ProteintoCOM(name)
            exportname, error = PDBFuncs.CutSphericalEnvelope(name, out_center, args.RMax)
            new_filenamesProtein_env.append(exportname)
            if error:
                errorList.append(name)

        if len(errorList) > 0:
            print("In the following files the protein is too close to the envelope: ")
            for i in errorList:
                print(os.path.split(i)[1])
            print("Do you wish to continue the calculation? (yes/no)")
            ans = input()
            errorcount = 0
            while ans != "no" and ans != "yes" and errorcount < 3:
                print("Please Enter yes or no.")
                ans = input()
                errorcount = errorcount + 1
            if ans == "no":
                print("exit program....")
                Funcs.DeleteFiles(True, args.ProteinDir, args.SolventDir)
                exit()
            elif ans == "yes":
                print("Continue calculations...")
            else:
                print("exit program...")
                Funcs.DeleteFiles(True, args.ProteinDir, args.SolventDir)
                exit()

        for name in FinalSolventList:
            xyz, out_center = PDBFuncs.SolventtoCOM(name)
            exportname, error = PDBFuncs.CutSphericalEnvelope(name, out_center, args.RMax)
            new_filenamesSolvent_env.append(exportname)



        FinalProteinList = np.copy(new_filenamesProtein_env)
        FinalSolventList = np.copy(new_filenamesSolvent_env)

    # optional - if the frames need to be rotated
    if (args.Rotate == True):
        print("Aligning the protein frames...")
        new_filenamesProtein_Rotated = []
        new_filenamesSolvent_Rotated = []
        Rotations = []
        filename_ref, ref_xyz, out_center = PDBFuncs.read_pdb_Protein_xyz(FinalProteinList[0])
        for name in FinalProteinList:
            outputfile, alpha, beta, gamma = PDBFuncs.ExportRotatedPDB(name, ref_xyz)
            Rotations.append([alpha, beta, gamma])
            new_filenamesProtein_Rotated.append(outputfile)
        i = 0
        for name in FinalSolventList:
            new_filenamesSolvent_Rotated.append(PDBFuncs.RotateSolvent(name, Rotations[i][0], Rotations[i][1], Rotations[i][2]))
            i = i + 1
        FinalProteinList = np.copy(new_filenamesProtein_Rotated)
        FinalSolventList = np.copy(new_filenamesSolvent_Rotated)

    #Creating a tree that will be used for the calculations of all the PDB files
    StateFile = "PDBCalculation.state"
    api = LocalRunner()
    PDBTree = Funcs.creatTree(StateFile,args.qMax, args.Grid)


    # Calculate the intensity and amplitude for each frame
    print("Starting the calculation for the frames...")
    ListOfIntensityProtein = []
    ListOfAmplitudeProtein = []

    for name in FinalProteinList:
        name_path, basename = os.path.split(name)
        PDBTree.state.Domain.Children[0].Children[0].filename = repr(name).replace("'", "")
        I, A = Funcs.CalcPDB(PDBTree, api, name_path, basename)
        ListOfAmplitudeProtein.append(A)
        ListOfIntensityProtein.append(I)


    ListOfIntensitySolvent = []
    ListOfAmplitudeSolvent = []
    for name in FinalSolventList:
        name_path, basename = os.path.split(name)
        PDBTree.state.Domain.Children[0].Children[0].filename = repr(name).replace("'", "")
        I, A = Funcs.CalcPDB(PDBTree, api, name_path, basename)
        ListOfIntensitySolvent.append(I)
        ListOfAmplitudeSolvent.append(A)

    print("The calculation is dine...")

    #Calculates the averages
    print("Calculating the averages...")
    Q_vec_Protein, Ave_intensity_Protein = Funcs.averageIntensities(ListOfIntensityProtein, "Protein_Averaged_Intensity.out")
    Q_vec_Solvent, Ave_intensity_Solvent = Funcs.averageIntensities(ListOfIntensitySolvent, "Solvent_Averaged_Intensity.out")

    if  args.Rotate == True and args.Envelope == False:
        print('This option was not completed. '
              'If you wish to align the protein please cut a sperical envelope around it.')

    else:
        Ave_Amplitude_Protein,Ave_Amplitude_Protein_Name, HeaderList_Protein  = Funcs.averageAmps(ListOfAmplitudeProtein, args.qMax,"Protein_Ave_Amps.amp")
        Ave_Amplitude_Solvent,Ave_Amplitude_Solvent_Name, HeaderList_Solvent  = Funcs.averageAmps(ListOfAmplitudeSolvent,args.qMax,"Solvent_Ave_Amps.amp")
        SubAveragedAmplitude = Ave_Amplitude_Protein - Ave_Amplitude_Solvent
        SubAveragedAmplitudeName = Funcs.writeAmp(SubAveragedAmplitude, HeaderList_Solvent, os.path.dirname(Ave_Amplitude_Solvent_Name) + "\\Sub_Ave_Amp.amp")


    #Calculate the intensities
    AmpFile = "AmpCalculation.state"
    AmpTree = Funcs.creatTree(AmpFile,args.qMax, args.Grid)

    Q_vec_Sub, SubAveragedIntensity = Funcs.CalculateIntensityFromAmplitud(SubAveragedAmplitudeName, AmpTree,api)
    Q_vec_Ave_Protein, IntensityOfTheAve_AmplituseProtein = Funcs.CalculateIntensityFromAmplitud(Ave_Amplitude_Protein_Name, AmpTree,api)
    Q_vec_Ave_Solvent, IntensityOfTheAve_AmplituseSolvent = Funcs.CalculateIntensityFromAmplitud(Ave_Amplitude_Solvent_Name, AmpTree,api)

    #Calculation of eq. 10 for the total averaged intensity
    TotalIntensity = np.asarray(SubAveragedIntensity, dtype=float) + np.asarray(Ave_intensity_Protein, dtype=float)
    TotalIntensity = TotalIntensity - np.asarray(IntensityOfTheAve_AmplituseProtein, dtype=float) - np.asarray(Ave_intensity_Solvent, dtype=float)
    TotalIntensity = TotalIntensity + np.asarray(IntensityOfTheAve_AmplituseSolvent, dtype=float)

    #export result
    Funcs.ExpotrIntensity(TotalIntensity, Q_vec_Protein)

    Funcs.DeleteFiles(args.delete, args.ProteinDir, args.SolventDir)

if __name__ == "__main__":
    main()
