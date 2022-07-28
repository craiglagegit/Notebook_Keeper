import os, time
from subprocess import Popen
# This code syncs the notebooks and copies the markdown files over for pushing to github.

dirs = [["/project/cslage/BOT_LSSTCam/notebooks", "/home/cslage/alternate_branches/Notebook_Keeper/bot_markdown/"], \
        ["/project/cslage/AuxTel/notebooks", "/home/cslage/alternate_branches/Notebook_Keeper/auxtel_markdown/"], \
        ["/project/cslage/ComCam/notebooks", "/home/cslage/alternate_branches/Notebook_Keeper/comcam_markdown/"], \
        ["/project/cslage/BOT_LSSTCam/notebooks/reca", "/home/cslage/alternate_branches/Notebook_Keeper/reca_markdown/"]]
for [get_dir, put_dir] in dirs:
        cmd2 = "jupytext --sync %s/*.ipynb"%get_dir
        ex_cmd2 = Popen(cmd2, shell=True)
        Popen.wait(ex_cmd2)
        cmd3 = "cp %s/markdown/* %s/"%(get_dir, put_dir)
        ex_cmd3 = Popen(cmd3, shell=True)
        Popen.wait(ex_cmd3)

print("Done")
