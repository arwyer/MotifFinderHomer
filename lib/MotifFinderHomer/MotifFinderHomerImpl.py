# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import uuid
from installed_clients.KBaseReportClient import KBaseReport
import json
from Bio import SeqIO
from pprint import pprint, pformat
from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.SequenceSetUtilsClient import SequenceSetUtils
from installed_clients.MotifUtilsClient import MotifUtils
import subprocess
import re
from pprint import pprint, pformat
from datetime import datetime
import uuid
from MotifFinderHomer.Utils.MotifSetUtil import MotifSetUtil
#import MotifFinderHomer.Utils.HomerUtil as HU
from MotifFinderHomer.Utils.HomerUtil import HomerUtil
from MotifFinderHomer.Utils.MakeNewReport import MakeNewReport
import subprocess
from biokbase.workspace.client import Workspace
#from MotifFinderHomer.Utils.FastaUtils import RemoveRepeats
from MotifFinderHomer.Utils.FastaUtils import FastaUtils
from MotifFinderHomer.Utils.BackgroundUtils import BackgroundUtils
from installed_clients.AssemblyUtilClient import AssemblyUtil
from MotifFinderHomer.Utils.TestUtils import TestUtils
#END_HEADER


class MotifFinderHomer:
    '''
    Module Name:
    MotifFinderHomer

    Module Description:
    A KBase module: MotifFinderHomer
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        #END_CONSTRUCTOR
        pass


    def find_motifs(self, ctx, params):
        """
        :param params: instance of type "find_motifs_params" (SS_ref -
           optional, used for exact genome locations if possible) ->
           structure: parameter "workspace_name" of String, parameter
           "fastapath" of String, parameter "motif_min_length" of Long,
           parameter "motif_max_length" of Long, parameter "SS_ref" of
           String, parameter "obj_name" of String
        :returns: instance of type "extract_output_params" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN find_motifs
        if 'motif_min_length' not in params:
            params['motif_min_length'] = 8
        if 'motif_max_length' not in params:
            params['motif_max_length'] = 16
        motMin = params['motif_min_length']
        motMax = params['motif_max_length']

        promoterFastaFilePath = params['fastapath']

        HU=HomerUtil()
        HomerMotifCommand = HU.build_homer_motif_command(promoterFastaFilePath,motMin,motMax,params['background'])
        HomerLocationCommand = HU.build_homer_location_command(promoterFastaFilePath)
        HU.run_homer_command(HomerMotifCommand)
        HU.run_homer_command(HomerLocationCommand)
        homer_out_path = '/kb/module/work/tmp/homer_out'
        #homer_out_path = os.path.join(self.shared_folder, str(uuid.uuid4()))
        homer_params = {'ws_name' : params['workspace_name'], 'path' : homer_out_path + '/homerMotifs.all.motifs','location_path' : homer_out_path + '/homer_locations.txt','obj_name' : params['obj_name']}
        MOU = MotifUtils(self.callback_url)
        dfu = DataFileUtil(self.callback_url)
        locDict = {}
        if 'SS_ref' in params:
            get_ss_params = {'object_refs' : [params['SS_ref']]}
            SS = dfu.get_objects(get_ss_params)['data'][0]['data']
            for s in SS['sequences']:
                if s['source']['assembly_id'] != '' and s['source']['location'] != []:
                    locDict['sequence_id'] = {'contig' : s['source']['location'][0][0],'start':str(s['source']['location'][0][1])}
        if len(locDict.keys()) > 0:
            homer_params['absolute_locations'] = locDict
        homer_params['min_len'] = motMin
        homer_params['max_len'] = motMax
        obj_ref = MOU.UploadFromHomer(homer_params)['obj_ref']
      
        HU.write_obj_ref(homer_out_path, obj_ref)
        
        #file = open(homer_out_path+"/homer_obj.txt","w") 
        #file.write(obj_ref) 
        #file.close() 
        #HERE:
        #we've got object ref
        #we've got html building functions
        #build report, setup return,
        #make report and return it


        timestamp = int((datetime.utcnow() - datetime.utcfromtimestamp(0)).total_seconds()*1000)
        timestamp = str(timestamp)
        htmlDir = self.shared_folder + '/html' +  timestamp
        os.mkdir(htmlDir)
        lineCount = 0
        with open(promoterFastaFilePath,'r') as pFile:
            for line in pFile:
                lineCount += 1
        numFeat = lineCount/2
        with open(promoterFastaFilePath,'r') as pFile:
            fileStr = pFile.read()
        promHtmlStr = '<html><body> '  + fileStr + ' </body></html>'
        with open(htmlDir + '/promoters.html','w') as promHTML:
            promHTML.write(promHtmlStr)
        JsonPath = '/kb/module/work/tmp'

        dfu = DataFileUtil(self.callback_url)
        get_obj_params = {'object_refs' : [obj_ref]}
        homerMotifSet = dfu.get_objects(get_obj_params)['data'][0]['data']
        mr=MakeNewReport()
        mr.MakeReport(htmlDir,homerMotifSet)


        try:
            html_upload_ret = dfu.file_to_shock({'file_path': htmlDir ,'make_handle': 0, 'pack': 'zip'})
        except:
            raise ValueError ('error uploading HTML file to shock')


        reportName = 'HomerMotifFinder_report_'+str(uuid.uuid4())

        reportObj = {'objects_created': [{'ref' : obj_ref, 'description' : 'Motif Set generated by Homer'}],
                     'message': '',
                     'direct_html': None,
                     'direct_html_link_index': 0,
                     'file_links': [],
                     'html_links': [],
                     'html_window_height': 220,
                     'workspace_name': params['workspace_name'],
                     'report_object_name': reportName
                     }


        # attach to report obj
        reportObj['direct_html'] = ''
        reportObj['direct_html_link_index'] = 0
        reportObj['html_links'] = [{'shock_id': html_upload_ret['shock_id'],
                                    #'name': 'promoter_download.zip',
                                    'name': 'index.html',
                                    'label': 'Save promoter_download.zip'
                                    }
                                   ]


        report = KBaseReport(self.callback_url, token=ctx['token'])
        report_info = report.create_extended_report(reportObj)
        output = { 'report_name': report_info['name'], 'report_ref': report_info['ref'] }

        #END find_motifs

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method find_motifs return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def BuildFastaFromSequenceSet(self, ctx, params):
        """
        :param params: instance of type "BuildSeqIn" -> structure: parameter
           "workspace_name" of String, parameter "SequenceSetRef" of String,
           parameter "fasta_outpath" of String
        :returns: instance of type "BuildSeqOut" -> structure: parameter
           "fasta_outpath" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN BuildFastaFromSequenceSet
        dfu = DataFileUtil(self.callback_url)

        bu = BackgroundUtils()
        TU = TestUtils()
        if params['TESTFLAG'] and params['background']:
            targetpath = '/kb/module/work/tmp/testgenome.fa'
            TU.GetGenome(targetpath)
            bu.BuildBackground(targetpath)
        elif params['background']:

            ws = Workspace('https://appdev.kbase.us/services/ws')
            subset = ws.get_object_subset([{
                                         'included':['/features/[*]/location', '/features/[*]/id','/assembly_ref'],
    'ref':params['genome_ref']}])
            aref = subset[0]['data']['assembly_ref']
            assembly_ref = {'ref': aref}
            print('Downloading Assembly data as a Fasta file.')
            assemblyUtil = AssemblyUtil(self.callback_url)
            fasta_file = assemblyUtil.get_assembly_as_fasta(assembly_ref)['path']
            bu.BuildBackground(fasta_file)


        get_objects_params = {'object_refs' : [params['SequenceSetRef']]}

        SeqSet = dfu.get_objects(get_objects_params)['data'][0]['data']
        outFile = open(params['fasta_outpath'],'w')
        for s in SeqSet['sequences']:
            sname = '>' + s['sequence_id'] + '\n'
            outFile.write(sname)
            sseq = s['sequence'] + '\n'
            outFile.write(sseq)
        outFile.close()

        fu=FastaUtils()
        if params['mask_repeats']:
            fu.RemoveRepeats(params['fasta_outpath'],params['fasta_outpath'])

        output = {'fasta_outpath' : params['fasta_outpath']}
        #END BuildFastaFromSequenceSet

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method BuildFastaFromSequenceSet return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def ExtractPromotersFromFeatureSetandDiscoverMotifs(self, ctx, params):
        """
        :param params: instance of type "extract_input" -> structure:
           parameter "workspace_name" of String, parameter "genome_ref" of
           String, parameter "featureSet_ref" of String, parameter
           "promoter_length" of Long, parameter "motif_min_length" of Long,
           parameter "motif_max_length" of Long, parameter "obj_name" of
           String
        :returns: instance of type "extract_output_params" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN ExtractPromotersFromFeatureSetandDiscoverMotifs
        SSU = SequenceSetUtils(self.callback_url)

        BuildParams = {'ws_name' : params['workspace_name'], 'FeatureSet_ref' : params['featureSet_ref'], 'genome_ref' : params['genome_ref'], 'upstream_length' : params['promoter_length']}
        SSret =  SSU.buildFromFeatureSet(BuildParams)
        SSref = SSret['SequenceSet_ref']
        fastapath = '/kb/module/work/tmp/tmpSeqSet.fa'
        newfastapath = '/kb/module/work/tmp/SeqSet.fa'
        fastapath = newfastapath
        FastaParams = {'workspace_name' : params['workspace_name'] , 'SequenceSetRef' : SSref , 'fasta_outpath' : fastapath}

        output = self.BuildFastaFromSequenceSet(ctx,FastaParams)

        findmotifsparams= {'workspace_name' : params['workspace_name'],'fastapath':fastapath,'motif_min_length':params['motif_min_length'],'motif_max_length':params['motif_max_length'],'SS_ref':SSref,'obj_name' :params['obj_name']}

        output = self.find_motifs(ctx,findmotifsparams)[0]
        #END ExtractPromotersFromFeatureSetandDiscoverMotifs

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method ExtractPromotersFromFeatureSetandDiscoverMotifs return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def DiscoverMotifsFromFasta(self, ctx, params):
        """
        :param params: instance of type "discover_fasta_input" -> structure:
           parameter "workspace_name" of String, parameter "fasta_path" of
           String
        :returns: instance of type "extract_output_params" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN DiscoverMotifsFromFasta
        #END DiscoverMotifsFromFasta

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method DiscoverMotifsFromFasta return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def DiscoverMotifsFromSequenceSet(self, ctx, params):
        """
        :param params: instance of type "discover_seq_input" -> structure:
           parameter "workspace_name" of String, parameter "genome_ref" of
           String, parameter "SS_ref" of String, parameter "promoter_length"
           of Long, parameter "motif_min_length" of Long, parameter
           "motif_max_length" of Long, parameter "obj_name" of String,
           parameter "background" of Long, parameter "mask_repeats" of Long,
           parameter "background_group" of mapping from String to String
        :returns: instance of type "extract_output_params" -> structure:
           parameter "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN DiscoverMotifsFromSequenceSet
        fastapath = '/kb/module/work/tmp/tmpSeqSet.fa'
        newfastapath = '/kb/module/work/tmp/SeqSet.fa'
        fastapath = newfastapath

        if params['background_group'] == None:
            params['background_group'] = {'background': 0}

        FastaParams = {
            'workspace_name': params['workspace_name'],
            'SequenceSetRef' : params['SS_ref'],
            'fasta_outpath' : fastapath,
            'background':params['background_group']['background'],
            'mask_repeats':params['mask_repeats']
        }

        if params['background_group']['background'] == 1:
            FastaParams['genome_ref'] = params['background_group']['genome_ref']
        else:
            FastaParams['genome_ref'] = 'NULL'
        if 'TESTFLAG' in params:
            FastaParams['TESTFLAG'] = params['TESTFLAG']
        else:
            FastaParams['TESTFLAG'] = 0
        output = self.BuildFastaFromSequenceSet(ctx,FastaParams)

        #RemoveRepeats(fastapath,newfastapath)
        findmotifsparams= {'workspace_name' : params['workspace_name'],'fastapath':fastapath,'motif_min_length':params['motif_min_length'],'motif_max_length':params['motif_max_length'],'SS_ref':params['SS_ref'],'obj_name':params['obj_name']}
        if params['background_group']['background'] == 1:
            findmotifsparams['background'] = 1
        else:
            findmotifsparams['background'] = 0

        output = self.find_motifs(ctx,findmotifsparams)[0]
        #END DiscoverMotifsFromSequenceSet

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method DiscoverMotifsFromSequenceSet return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
