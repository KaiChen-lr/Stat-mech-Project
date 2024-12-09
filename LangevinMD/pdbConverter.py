import numpy as np

def writePdb(positions, outputFile):
    with open(outputFile, 'w') as pdbFile:

        atomNumber = 1

        for i in range(len(positions)):
            x, y, z = positions[i]

            if i < 108:
                atomName = "AR"
                residueName = "ARG"
            else:
                atomName = "NE"
                residueName = "NEO"

            pdbLine = f"ATOM  {atomNumber:5}  {atomName}  {residueName} A{atomNumber:4}    {x:8.3f} {y:7.3f} {z:7.3f}  1.00 20.00           {atomName[0]}\n"

            pdbFile.write(pdbLine)

            atomNumber += 1

        pdbFile.write("END\n")


initCfg=np.loadtxt("InitialCfg.txt")
finalCfg=np.loadtxt("FinalCfg.txt")

outputPdb1 = "InitialCfg.pdb"
outputPdb2 = "FinalCfg.pdb"

writePdb(initCfg, outputPdb1)
writePdb(finalCfg, outputPdb2)
