import subprocess
import os

def BuildBackground(fastapath):
    homerpath = '/kb/module/work/tmp/homer_background.fa'
    backgroundCommand = 'mv ' + fastapath + ' ' + homerpath

    try:
        out_txt = subprocess.check_output(backgroundCommand,shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print('************GET BACKGROUND ERROR************\n')
        print(e.returncode)
        exit(0)
    assert os.path.isfile(homerpath)
