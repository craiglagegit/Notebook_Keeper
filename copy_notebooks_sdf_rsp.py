import sys, os, time
from pathlib import Path
from subprocess import Popen
# This code strips the data from the notebooks, and copies them over for pushing to github.

dirs = [["/sdf/group/rubin/u/cslage/BOT_LSSTCam/notebooks", "/home/c/cslage/WORK/Notebook_Keeper/bot_notebooks/"], \
        ["/sdf/group/rubin/u/cslage/AuxTel/notebooks", "/home/c/cslage/WORK/Notebook_Keeper/auxtel_notebooks/"], \
        ["/sdf/group/rubin/u/cslage/ComCam/notebooks", "/home/c/cslage/WORK/Notebook_Keeper/comcam_notebooks/"], \
        ["/sdf/group/rubin/u/cslage/BOT_LSSTCam/notebooks/reca", "/home/c/cslage/WORK/Notebook_Keeper/reca_notebooks/"]]
for [get_dir, put_dir] in dirs:
        files = os.listdir(get_dir)
        filesToStore = []
        for file in files:
                if file.split(".")[-1] == "ipynb":
                        thisFile = f"{get_dir}/{file}"
                        storedFile = f"{put_dir}/{file}"
                        thisPath = Path(thisFile)
                        thisLastModified = thisPath.stat().st_mtime
                        try:
                                storedPath = Path(storedFile)
                                storedLastModified = storedPath.stat().st_mtime
                        except:
                                filesToStore.append(thisFile)
                        if thisLastModified > storedLastModified:
                                filesToStore.append(thisFile)
        for thisFile in filesToStore:
                print("Stripping and copying ", thisFile)
                command = f"jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True \
                --output-dir={put_dir} {thisFile}"
                strip_and_copy = Popen(command, shell=True)
                Popen.wait(strip_and_copy)
print("Done")

