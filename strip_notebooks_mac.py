import os, sys, time
from subprocess import Popen

# This code strips the data from notebook files in preparation for archiving them
# on github.

fromDir = "/Users/cslage/.ipython/notebooks"
toDir = "/Users/cslage/Research/LSST/code/Notebook_Keeper/mac_notebooks"

files = os.listdir(fromDir)
for file in files:
    firstCharacter = list(file)[0]
    if firstCharacter == "." or firstCharacter == '_':
        continue

    fromFile = f"{fromDir}/{file}"
    toFile = f"{toDir}/{file}"
    if os.path.isfile(fromFile):
        if fromFile.split(".")[-1] == "ipynb":
            print("Stripping ", fromFile)
            command = f"jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True --output {toFile} {fromFile}"
            strip = Popen(command, shell=True)
            Popen.wait(strip)
        else:
            command = f"cp {fromFile} {toFile}"
            copy = Popen(command, shell=True)
            Popen.wait(copy)

    elif os.path.isdir(fromFile):
        print(f"Creating directory {toFile}")
        command = f"mkdir -p {toFile}" 
        make = Popen(command, shell=True)
        Popen.wait(make)

        subFiles = os.listdir(fromFile)
        for subFile in subFiles:
            firstCharacter = list(subFile)[0]
            if firstCharacter == "." or firstCharacter == '_':
                continue

            subFromFile = f"{fromDir}/{file}/{subFile}"
            subToFile = f"{toDir}/{file}/{subFile}"        
            #print(subFromFile)
            #print(subToFile)
            #sys.exit()

            if os.path.isfile(subFromFile):
                if subFromFile.split(".")[-1] == "ipynb":
                    print("Stripping ", subFromFile)
                    command = f"jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True --output {subToFile} {subFromFile}"
                    strip = Popen(command, shell=True)
                    Popen.wait(strip)
                else:
                    command = f"cp {subFromFile} {subToFile}"
                    copy = Popen(command, shell=True)
                    Popen.wait(copy)

            else:
                continue

print("Done")
