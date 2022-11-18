import sys, os, time
from subprocess import Popen
# This code strips the data from the notebooks, and copies them over for pushing to github.

dirs = [["/sdf/group/rubin/u/cslage/BOT_LSSTCam/notebooks", "/sdf/home/c/cslage/WORK/Notebook_Keeper/bot_notebooks/"], \
        ["/sdf/group/rubin/u/cslage/AuxTel/notebooks", "/sdf/home/c/cslage/WORK/Notebook_Keeper/auxtel_notebooks/"], \
        ["/sdf/group/rubin/u/cslage/ComCam/notebooks", "/sdf/home/c/cslage/WORK/Notebook_Keeper/comcam_notebooks/"], \
        ["/sdf/group/rubin/u/cslage/BOT_LSSTCam/notebooks/reca", "/sdf/home/c/cslage/WORK/Notebook_Keeper/reca_notebooks/"]]
for [get_dir, put_dir] in dirs:
        files = os.listdir(get_dir)
        for file in files:
                thisFile = f"{get_dir}/{file}"
                if file.split(".")[-1] == "ipynb":
                        print("Stripping and copying ", thisFile)
                        command = f"jupyter nbconvert --to notebook --ClearOutputPreprocessor.enabled=True \
                        --output-dir={put_dir} {thisFile}"
                        strip_and_copy = Popen(command, shell=True)
                        Popen.wait(strip_and_copy)
                else:
                        print("Skipping ", thisFile)

print("Done")
