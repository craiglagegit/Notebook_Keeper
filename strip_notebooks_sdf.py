import os, time
from subprocess import Popen

# This code strips the data from notebook files in preparation for archiving them
# on github.

dirs = os.listdir("./")
for dir in dirs:
        try:
                files = os.listdir(dir)
        except NotADirectoryError:
                continue
        for file in files:
                thisFile = f"{dir}/{file}"
                if file.split(".")[-1] == "ipynb":
                        print("Stripping ", thisFile)
                        command = "jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True --inplace %s"%thisFile
                        strip = Popen(command, shell=True)
                        Popen.wait(strip)
                else:
                        print("Skipping ", thisFile)

print("Done")
