import errno
import json
import multiprocessing
import os
import re
import subprocess
import sys
import time
import traceback
import uuid
import zipfile

from pathos.multiprocessing import ProcessingPool as Pool

from DataFileUtil.DataFileUtilClient import DataFileUtil
from KBaseReport.KBaseReportClient import KBaseReport
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from SetAPI.SetAPIServiceClient import SetAPI
from Workspace.WorkspaceClient import Workspace as Workspace
from kb_Bowtie2.kb_Bowtie2Client import kb_Bowtie2
from kb_QualiMap.kb_QualiMapClient import kb_QualiMap


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))  


class TopHatUtil:

    TOPHAT2_TOOLKIT_PATH = '/kb/deployment/bin/TopHat2'

    OPTIONS_MAP = {'read_mismatches': '--read-mismatches',
                   'read_gap_length': '--read-gap-length',
                   'read_edit_dist': '--read-edit-dist',
                   'min_intron_length': '--min-intron-length',
                   'max_intron_length': '--max-intron-length',
                   'min_anchor_length': '--min-anchor-length',
                   'report_secondary_alignments': '--report-secondary-alignments',
                   'no_coverage_search': '--no-coverage-search',
                   'library_type': '--library-type',
                   'num_threads': '--num-threads'
                   }

    BOOLEAN_OPTIONS = ['report_secondary_alignments', 'no_coverage_search']

    @staticmethod
    def _mkdir_p(path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    @staticmethod
    def _validate_run_tophat2_app_params(params):
        """
        _validate_run_tophat2_app_params:
                validates params passed to run_tophat2_app method
        """

        log('start validating run_tophat2_app params')

        # check for required parameters
        for p in ['input_ref', 'assembly_or_genome_ref', 'workspace_name', 
                  'alignment_suffix']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    @staticmethod
    def _run_command(command):
        """
        _run_command: run command and print result
        """

        log('start executing command:\n{}'.format(command))

        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output = pipe.communicate()[0]
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed commend:\n{}\n'.format(command) +
                'Exit Code: {}\nOutput:\n{}'.format(exitCode, output))
        else:
            error_msg = 'Error running commend:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}'.format(exitCode, output)
            raise ValueError(error_msg)

    def _get_bowtie_index(self, result_directory, assembly_or_genome_ref, workspace_name):
        """
        _get_bowtie_index: gets genome index file using kb_Bowtie2
        """

        log('start generating genome index file using kb_Bowtie2')

        output_dir = os.path.join(result_directory, 'bowtie2_index_' + str(int(time.time() * 100)))
        get_bowtie_index_params = {'ref': assembly_or_genome_ref,
                                   'output_dir': output_dir,
                                   'ws_for_cache': workspace_name}

        genome_index_file_dir = self.bt.get_bowtie2_index(get_bowtie_index_params)['output_dir']

        return genome_index_file_dir

    @staticmethod
    def _get_type_from_obj_info(info):
        """
        _get_type_from_obj_info: parses object type from object info
        """
        return info[2].split('-')[0]

    def _get_input_object_info(self, input_ref):
        """
        _get_input_object_info: gets input object data type and info
        """

        info = self.ws.get_object_info3({'objects': [{'ref': input_ref}]})['infos'][0]

        obj_type = self._get_type_from_obj_info(info)

        if obj_type in ['KBaseAssembly.PairedEndLibrary', 'KBaseAssembly.SingleEndLibrary',
                        'KBaseFile.PairedEndLibrary', 'KBaseFile.SingleEndLibrary']:
            return {'run_mode': 'single_library', 'info': info, 'ref': input_ref}
        elif obj_type == 'KBaseRNASeq.RNASeqSampleSet':
            return {'run_mode': 'sample_set', 'info': info, 'ref': input_ref}
        elif obj_type == 'KBaseSets.ReadsSet':
            return {'run_mode': 'sample_set', 'info': info, 'ref': input_ref}
        else:
            raise ValueError('Object type of input_ref is not valid, was: ' + str(obj_type))

    def _get_reads_file(self, reads_ref, reads_type, result_directory):
        """
        _get_reads_file: gets reads file from Single/Paired End Libiary
        """

        reads_file_paths = []

        download_reads_params = {'read_libraries': [reads_ref],
                                 'interleaved': 'false',
                                 'gzipped': None}

        reads_files = self.ru.download_reads(download_reads_params)['files'][reads_ref]['files']

        reads_file_dir = os.path.join(result_directory, 
                                      'reads_file_' + str(int(time.time() * 100)))
        self._mkdir_p(reads_file_dir)

        if reads_type.split('.')[1] == 'SingleEndLibrary':
            se_file_path = os.path.join(reads_file_dir, 'SE_reads.fastq')
            os.rename(reads_files['fwd'], se_file_path)
            reads_file_paths.append(se_file_path)
        elif reads_type.split('.')[1] == 'PairedEndLibrary':
            pe_fwd_file_path = os.path.join(reads_file_dir, 'PE_reads_1.fastq')
            pe_rev_file_path = os.path.join(reads_file_dir, 'PE_reads_2.fastq')
            os.rename(reads_files['fwd'], pe_fwd_file_path)
            os.rename(reads_files['rev'], pe_rev_file_path)
            reads_file_paths.append(pe_fwd_file_path)
            reads_file_paths.append(pe_rev_file_path)

        return reads_file_paths

    def _generate_command(self, genome_index_base, reads_files, 
                          tophat_result_dir, cli_option_params):
        """
        _generate_command: generate tophat2 command
        """

        command = self.TOPHAT2_TOOLKIT_PATH + '/tophat '

        command += '-o {} '.format(tophat_result_dir)

        for key, option in self.OPTIONS_MAP.items():
            option_value = cli_option_params.get(key)
            if key in self.BOOLEAN_OPTIONS and option_value:
                option_value = ' '
            if option_value:
                command += '{} {} '.format(option, option_value)

        preset_options = cli_option_params.get('preset_options')
        if preset_options:
            command += '{} '.format('--' + preset_options)

        command += '{} {}'.format(genome_index_base, ' '.join(reads_files))

        log('generated TopHat2 command: {}'.format(command))

        return command

    def _save_alignment(self, tophat_result_dir, alignment_name, reads_ref,
                        assembly_or_genome_ref, workspace_name, reads_condition):
        """
        _save_alignment: upload Alignment object
        """

        log('starting saving ReadsAlignment object')

        bam_file_path = self._merge_bam_files(tophat_result_dir)
        destination_ref = workspace_name + '/' + alignment_name
        if reads_condition:
            condition = reads_condition
        else:
            condition = 'unspecified'

        upload_alignment_params = {'file_path': bam_file_path,
                                   'destination_ref': destination_ref,
                                   'read_library_ref': reads_ref,
                                   'assembly_or_genome_ref': assembly_or_genome_ref,
                                   'aligned_using': 'tophat2',
                                   'aligner_version': '2.1.1',
                                   'condition': condition}

        reads_alignment_object_ref = self.rau.upload_alignment(upload_alignment_params)['obj_ref']

        return reads_alignment_object_ref

    def _merge_bam_files(self, tophat_result_dir, merged_file_name="merged_hits.bam"):
        """
        Tophat splits results into a mapped file and unmapped file while the alignment
        upload expects these to be in one file (like hisat and bowtie produces). This uses
        samtools to merge the files.
        """
        command = 'samtools merge {}/{} {}/accepted_hits.bam {}/unmapped.bam'.format(
            tophat_result_dir, merged_file_name, tophat_result_dir, tophat_result_dir)
        self._run_command(command)
        return os.path.join(tophat_result_dir, merged_file_name)

    def _save_alignment_set(self, reads_alignment_object_refs, workspace_name, alignment_set_name,
                            conditions):
        """
        _save_alignment_set: upload AlignmentSet object
        """

        log('starting saving ReadsAlignmentSet object')

        items = []
        for reads_alignment_object_ref, condition in zip(reads_alignment_object_refs, conditions):
            items.append({'ref': reads_alignment_object_ref, 'label': condition})
        alignment_set_data = {'description': 'Alignments using TopHat2', 
                              'items': items}

        alignment_set_save_params = {'data': alignment_set_data,
                                     'workspace': workspace_name,
                                     'output_object_name': alignment_set_name}

        save_result = self.set_client.save_reads_alignment_set_v1(alignment_set_save_params)
        alignment_set_object_ref = save_result['set_ref']

        return alignment_set_object_ref

    def _process_single_reads_library(self, input_object_info, genome_index_base, 
                                      result_directory, cli_option_params):
        """
        _process_single_reads_library: process single reads library
        """
        try:
            reads_obj_type = self._get_type_from_obj_info(input_object_info['info'])
            reads_obj_name = input_object_info['info'][1]
            reads_files = self._get_reads_file(input_object_info['ref'], 
                                               reads_obj_type, 
                                               result_directory)
            
            tophat_result_dir = os.path.join(result_directory, 
                                             'tophat2_result_' + reads_obj_name + 
                                             '_' + str(int(time.time() * 100)))
            command = self._generate_command(genome_index_base, reads_files, 
                                             tophat_result_dir, cli_option_params)
            self._run_command(command)

            alignment_object_name = reads_obj_name + cli_option_params.get('alignment_suffix')
            assembly_or_genome_ref = cli_option_params.get('assembly_or_genome_ref')
            reads_alignment_object_ref = self._save_alignment(tophat_result_dir,
                                                              alignment_object_name,
                                                              input_object_info['ref'],
                                                              assembly_or_genome_ref,
                                                              cli_option_params.get('workspace_name'),
                                                              cli_option_params.get('reads_condition'))
        except:
            log('caught exception in worker')
            e = sys.exc_info()[0]

            error_msg = 'ERROR -- {}: {}'.format(e, ''.join(traceback.format_stack()))

            reads_alignment_object_ref = error_msg
        finally:
            return reads_alignment_object_ref

    def _generate_report_single_library(self, reads_alignment_object_ref, result_directory, 
                                        workspace_name):
        """
        _generate_report_single_library: generate summary report for single library
        """

        log('start creating report')

        output_files = self._generate_output_file_list_single_library(result_directory)
        output_html_files = self._generate_html_report(reads_alignment_object_ref)

        description = 'Alignment generated by TopHat2'
        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'objects_created': [{'ref': reads_alignment_object_ref,
                                              'description': description}],
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 333,
                         'report_object_name': 'kb_tophat2_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_report_sets_library(self, reads_alignment_object_ref, result_directory, 
                                      workspace_name):
        """
        _generate_report_sets_library: generate summary report for sample sets
        """

        objects_created = [{'ref': reads_alignment_object_ref,
                            'description': 'AlignmentSet generated by TopHat2'}]
        alignment_set_data = self.ws.get_objects2({'objects': 
                                                  [{'ref':
                                                   reads_alignment_object_ref}]})['data'][0]
        alignment_refs = alignment_set_data['data'].get('items')
        for alignment_ref in alignment_refs:
            objects_created.append({'ref': alignment_ref['ref'],
                                    'description': 'Alignment generated by TopHat2'})

        output_files = self._generate_output_file_list_sets_library(result_directory)
        output_html_files = self._generate_html_report(reads_alignment_object_ref)

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'objects_created': objects_created,
                         'html_links': output_html_files,
                         'direct_html_link_index': 0,
                         'html_window_height': 333,
                         'report_object_name': 'kb_tophat2_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output

    def _generate_html_report(self, reads_alignment_object_ref):
        """
        _generate_html_report: generate html summary report
        """

        log('start generating html report')

        # running qualimap
        qualimap_report = self.qualimap.run_bamqc({'input_ref': reads_alignment_object_ref})
        qc_result_zip_info = qualimap_report['qc_result_zip_info']

        html_report = list()

        html_report.append({'shock_id': qc_result_zip_info['shock_id'],
                            'name': qc_result_zip_info['index_html_file_name'],
                            'label': qc_result_zip_info['name']})

        return html_report

    @staticmethod
    def _generate_output_file_list_single_library(result_directory):
        """
        _generate_output_file_list_single_library: zip result files and generate file_links 
                                                   for report
        """

        log('start packing result files')

        output_files = list()
        result_file = os.path.join(result_directory, 'TopHat2_result.zip')

        result_dirs = os.listdir(result_directory)
        tophat2_result_dir_name = filter(re.compile('tophat2_result_*').match, result_dirs)[0]
        tophat2_result_dir = os.path.join(result_directory, tophat2_result_dir_name)

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(tophat2_result_dir):
                for file in files:
                    if not (file.endswith('.DS_Store')):
                        zip_file.write(os.path.join(root, file), file)

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File generated by TopHat2 App'})

        return output_files

    @staticmethod
    def _generate_output_file_list_sets_library(result_directory):
        """
        _generate_output_file_list_sets_library: zip result files and generate file_links 
                                                 for report
        """

        log('start packing result files')

        output_files = list()
        result_dirs = os.listdir(result_directory)
        tophat2_result_dir_names = filter(re.compile('tophat2_result_*').match, result_dirs)
        result_files = list()

        for tophat2_result_dir_name in tophat2_result_dir_names:
            tophat2_result_dir = os.path.join(result_directory, tophat2_result_dir_name)
            file_name = tophat2_result_dir_name.split('tophat2_result_')[1].rsplit('_', 1)[0]
            result_file = result_file = os.path.join(result_directory, 
                                                     '{}.zip'.format(file_name))
            with zipfile.ZipFile(result_file, 'w',
                                 zipfile.ZIP_DEFLATED,
                                 allowZip64=True) as zip_file:
                for root, dirs, files in os.walk(tophat2_result_dir):
                    for file in files:
                        if not (file.endswith('.DS_Store')):
                            zip_file.write(os.path.join(root, file), file)
            result_files.append(result_file)

        for result_file in result_files:
            output_files.append({'path': result_file,
                                 'name': os.path.basename(result_file),
                                 'label': os.path.basename(result_file),
                                 'description': 'File generated by TopHat2 App'})

        return output_files

    def fetch_reads_refs_from_sampleset(self, ref, info):
        """
        Note: adapted from kbaseapps/kb_hisat2 - file_util.py

        From the given object ref, return a list of all reads objects that are a part of that
        object. E.g., if ref is a ReadsSet, return a list of all PairedEndLibrary or 
        SingleEndLibrary refs that are a member of that ReadsSet. 
        This is returned as a list of dictionaries as follows:
        {
            "ref": reads object reference,
            "condition": condition string associated with that reads object
        }
        The only one required is "ref", all other keys may or may not be present, based on the 
        reads object or object type in initial ref variable. E.g. a RNASeqSampleSet might have 
        condition info for each reads object, but a single PairedEndLibrary may not have that info.
        If ref is already a Reads library, just returns a list with ref as a single element.
        """
        obj_type = self._get_type_from_obj_info(info)
        refs = list()
        if "KBaseSets.ReadsSet" in obj_type or "KBaseRNASeq.RNASeqSampleSet" in obj_type:
            print("Looking up reads references in ReadsSet object")
            reads_set = self.set_client.get_reads_set_v1({'ref': ref,
                                                          'include_item_info': 0,
                                                          'include_set_item_ref_paths': 1})
            for reads in reads_set["data"]["items"]:
                refs.append({'ref': reads['ref_path'],
                             'condition': reads['label']
                             })
        else:
            raise ValueError("Unable to fetch reads reference from object {} "
                             "which is a {}".format(ref, obj_type))

        return refs

    def _process_set_reads_library(self, input_object_info, genome_index_base, 
                                   result_directory, cli_option_params):
        """
        _process_set_reads_library: process set reads library
        """

        reads_refs = self.fetch_reads_refs_from_sampleset(input_object_info['ref'],
                                                          input_object_info['info'])

        set_object_name = input_object_info['info'][1]
        alignment_set_name = set_object_name + cli_option_params['alignment_set_suffix']

        arg_1 = []
        arg_2 = [genome_index_base] * len(reads_refs)
        arg_3 = [result_directory] * len(reads_refs)
        arg_4 = []
        conditions = []
        for reads_ref in reads_refs:
            reads_input_object_info = self._get_input_object_info(reads_ref['ref'])
            option_params = cli_option_params.copy()
            option_params['reads_condition'] = reads_ref['condition']
            conditions.append(reads_ref['condition'])
            arg_1.append(reads_input_object_info)
            arg_4.append(option_params)

        cpus = min(cli_option_params.get('num_threads'), multiprocessing.cpu_count())
        pool = Pool(ncpus=cpus)
        log('running _process_alignment_object with {} cpus'.format(cpus))

        reads_alignment_object_refs = pool.map(self._process_single_reads_library, 
                                               arg_1, arg_2, arg_3, arg_4)

        for reads_alignment_object_ref in reads_alignment_object_refs:
            if reads_alignment_object_ref.startswith('ERROR'):
                error_msg = 'Caught exception in worker\n'
                error_msg += '{}'.format(reads_alignment_object_ref)
                raise ValueError(error_msg)

        workspace_name = cli_option_params['workspace_name']
        reads_alignment_set_object_ref = self._save_alignment_set(reads_alignment_object_refs,
                                                                  workspace_name,
                                                                  alignment_set_name,
                                                                  conditions)

        return reads_alignment_set_object_ref

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.scratch = config['scratch']
        self.srv_wiz_url = config['srv-wiz-url']
        self.ws = Workspace(self.ws_url, token=self.token)
        self.bt = kb_Bowtie2(self.callback_url)
        self.rau = ReadsAlignmentUtils(self.callback_url)
        self.qualimap = kb_QualiMap(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)
        self.set_client = SetAPI(self.srv_wiz_url)

    def run_tophat2_app(self, params):
        """
        run_tophat2_app: run TopHat2 app
        (https://ccb.jhu.edu/software/tophat/manual.shtml)

        required params:
        input_ref: input reads object (Single/Paired_reads, reads_set, sample_set)
        assembly_or_genome_ref: ref to Assembly, ContigSet, or Genome
        workspace_name: the name of the workspace it gets saved to
        alignment_set_suffix: suffix append to alignment set object name
        alignment_suffix: suffix append to alignment object name

        optional params:
        reads_condition: condition associated with the input reads objec (ignored for sets of 
                         samples)
        num_threads: number of processing threads
        read_mismatches: read mismatch cutoff
        read_gap_length: read gap cutoff
        read_edit_dist: read edit cutoff
        min_intron_length: minimum intron length
        max_intron_length: maximum intron length
        min_anchor_length: minimum anchor length
        report_secondary_alignments: use this option to output secondary alignments
        no_coverage_search: use this option to disable the coverage-based search for junctions
        library_type: library type (fr-unstranded, fr-firststrand, fr-secondstrand)
        preset_options: alignment preset options (b2-very-fast, b2-fast, b2-sensitive, 
                                                  b2-very-sensitive)

        return:
        result_directory: folder path that holds all files generated by run_tophat2_app
        reads_alignment_object_ref: generated Alignment/AlignmentSet object reference
        report_name: report name generated by KBaseReport
        report_ref: report reference generated by KBaseReport
        """

        log('--->\nrunning TopHatUtil.run_tophat2_app\n' +
            'params:\n{}'.format(json.dumps(params, indent=1)))

        self._validate_run_tophat2_app_params(params)

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(result_directory)

        genome_index_file_dir = self._get_bowtie_index(result_directory, 
                                                       params.get('assembly_or_genome_ref'),
                                                       params.get('workspace_name'))
        genome_index_files = os.listdir(genome_index_file_dir)

        log('generated genome index files: {}'.format(genome_index_files)) 
        genome_index_file = filter(re.compile(".*\.\d\..*").match, genome_index_files)[0]
        if re.match('.*\.rev\.\d\..*', genome_index_file):
            genome_index_file_prefix = genome_index_file.split('rev')[0][:-1]
        else:
            genome_index_file_prefix = ''
            for prefix in genome_index_file.split('.'):
                if prefix.isdigit():
                    break
                else:
                    genome_index_file_prefix += '.' + prefix
            genome_index_file_prefix = genome_index_file_prefix[1:]

        genome_index_base = genome_index_file_dir + '/' + genome_index_file_prefix

        input_object_info = self._get_input_object_info(params.get('input_ref'))

        if input_object_info['run_mode'] == 'single_library':
            reads_alignment_object_ref = self._process_single_reads_library(input_object_info, 
                                                                            genome_index_base,
                                                                            result_directory,
                                                                            params)
            report_output = self._generate_report_single_library(reads_alignment_object_ref,
                                                                 result_directory,
                                                                 params.get('workspace_name'))
        elif input_object_info['run_mode'] == 'sample_set':
            reads_alignment_object_ref = self._process_set_reads_library(input_object_info,
                                                                         genome_index_base,
                                                                         result_directory,
                                                                         params)
            report_output = self._generate_report_sets_library(reads_alignment_object_ref,
                                                               result_directory,
                                                               params.get('workspace_name'))

        returnVal = {'result_directory': result_directory,
                     'reads_alignment_object_ref': reads_alignment_object_ref}

        returnVal.update(report_output)

        return returnVal
