import sys
import os
import json
import subprocess

def build_homer_motif_command(inputFilePath,min,max,background):
    if not os.path.exists('/kb/module/work/tmp/homer_out'):
        os.mkdir('/kb/module/work/tmp/homer_out')
    outputDirPath = '/kb/module/work/tmp/homer_out'
    #outputDirPath = './temp/homer_out'
    #command = 'findMotifs.pl ' + inputFilePath + ' fasta ' + outputDirPath
    #/kb/module/work/homer/bin/
    #command = 'findMotifs.pl ' + inputFilePath + ' fasta ' + outputDirPath +' -basic'
    command = '/kb/deployment/bin/findMotifs.pl ' + inputFilePath + ' fasta ' + outputDirPath +' -basic'
    homerpath = '/kb/module/work/tmp/homer_background.fa'
    if background == 1:
        command += ' -fasta ' + homerpath
    command += ' -len '
    numLens = 0
    for i in range(min,max+1,2):
        numLens += 1
        command += str(i)
        if i <= max:
            command += ','
    command += ' -S 5'
    return command

def build_homer_location_command(inputFilePath):
    outputDirPath = '/kb/module/work/tmp/homer_out'
    #outputDirPath = './temp/homer_out'
    outputFilePath = outputDirPath + '/homerMotifs.all.motifs'
    outputTo = outputDirPath + '/homer_locations.txt'
    command = '/kb/deployment/bin/scanMotifGenomeWide.pl ' + outputFilePath + ' ' + inputFilePath + ' > ' + outputTo
    #command = 'scanMotifGenomeWide.pl ' + outputFilePath + ' ' + inputFilePath + ' > ' + outputTo
    return command


def run_homer_command(command):
    try:
        subprocess.check_output(command,shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print('*******HOMER ERROR******** : ' + command)
        print(e.returncode)
    #os.system(command)

def parse_homer_output(path,location):
    #outputDirPath = '/kb/module/work/tmp/homer_out'
    #outputDirPath = './temp/homer_out'
    #comment no cache
    #outputFilePath = outputDirPath + '/homerMotifs.all.motifs'
    #locationFilePath = outputDirPath + '/homer_locations.txt'
    outputFilePath = path
    locationFilePath = location
    homerFile = open(outputFilePath,'r')
    motifList = []
    motifDict = {}
    pwmList = []
    for line in homerFile:
        if '>' in line:
            if len(motifDict) != 0:
                motifDict['pwm'] = pwmList
                pwmList = []
                motifDict['Locations'] = []
                motifList.append(motifDict)
                motifDict = {}
            elems = line.split()
            motif = elems[0].replace('>','')
            motifDict['Iupac_signature'] = motif
            p_val = float(elems[5].split(',')[2].split(':')[1])
            motifDict['p-value'] = p_val
        else:
            elems = line.split()
            rowList = []
            rowList.append(('A',float(elems[0])))
            rowList.append(('C',float(elems[1])))
            rowList.append(('G',float(elems[2])))
            rowList.append(('T',float(elems[3])))
            pwmList.append(rowList)

    locationFile = open(locationFilePath,'r')
    for line in locationFile:
        if len(line.split()) == 7:
            elems = line.split()
            motif = elems[0].split('-')[1]
            for m in motifList:
                if m['Iupac_signature'] == motif:
                    locList = []
                    locList.append(elems[1])
                    locList.append(elems[2])
                    locList.append(elems[3])
                    locList.append(elems[4])
                    m['Locations'].append(locList)
                    break
    jsonFilePath = outputDirPath + '/homer.json'
    with open(jsonFilePath,'w') as jsonFile:
        json.dump(motifList,jsonFile)
    return motifList
