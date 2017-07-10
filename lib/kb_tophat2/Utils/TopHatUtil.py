import time
import json
import os
import uuid
import errno
import subprocess
import re

from Workspace.WorkspaceClient import Workspace as Workspace
from kb_Bowtie2.kb_Bowtie2Client import kb_Bowtie2
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from ReadsUtils.ReadsUtilsClient import ReadsUtils


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

    def _mkdir_p(self, path):
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

    def _validate_run_tophat2_app_params(self, params):
        """
        _validate_run_tophat2_app_params:
                validates params passed to run_tophat2_app method
        """

        log('start validating run_tophat2_app params')

        # check for required parameters
        for p in ['input_ref', 'assembly_or_genome_ref', 'workspace_name', 
                  'alignment_object_name']:
            if p not in params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _run_command(self, command):
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

    def _get_type_from_obj_info(self, info):
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

        bam_file_path = os.path.join(tophat_result_dir, 'accepted_hits.bam')
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

    def _process_single_reads_library(self, input_object_info, genome_index_base, 
                                      result_directory, cli_option_params):
        """
        _process_single_reads_library: process single reads library
        """

        obj_type = self._get_type_from_obj_info(input_object_info['info'])
        reads_files = self._get_reads_file(input_object_info['ref'], obj_type, result_directory)
        
        tophat_result_dir = os.path.join(result_directory, 
                                         'tophat2_result_' + str(int(time.time() * 100)))
        command = self._generate_command(genome_index_base, reads_files, 
                                         tophat_result_dir, cli_option_params)
        self._run_command(command)

        reads_alignment_object_ref = self._save_alignment(tophat_result_dir,
                                                          cli_option_params.get('alignment_object_name'),
                                                          input_object_info['ref'],
                                                          cli_option_params.get('assembly_or_genome_ref'),
                                                          cli_option_params.get('workspace_name'),
                                                          cli_option_params.get('reads_condition'))

        return reads_alignment_object_ref

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.scratch = config['scratch']
        self.ws = Workspace(self.ws_url, token=self.token)

        self.bt = kb_Bowtie2(self.callback_url, service_ver='dev')
        self.ru = ReadsUtils(self.callback_url)
        self.rau = ReadsAlignmentUtils(self.callback_url, service_ver='dev')

    def run_tophat2_app(self, params):
        """
        run_tophat2_app: run TopHat2 app
        (https://ccb.jhu.edu/software/tophat/manual.shtml)

        required params:
        input_ref: input reads object (Single/Paired_reads, reads_set, sample_set)
        assembly_or_genome_ref: ref to Assembly, ContigSet, or Genome
        workspace_name: the name of the workspace it gets saved to
        alignment_object_name: output Alignment or AlignmentSet object name

        optional params:
        reads_condition: condition associated with the input reads objec (ignored for sets of samples)
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
        preset_options: alignment preset options (b2-very-fast, b2-fast, b2-sensitive, b2-very-sensitive)

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
        genome_index_file_prefix = filter(re.compile(".*\.\d\..*").match, 
                                          genome_index_files)[0].split('.')[0]
        genome_index_base = genome_index_file_dir + '/' + genome_index_file_prefix

        input_object_info = self._get_input_object_info(params.get('input_ref'))

        if input_object_info['run_mode'] == 'single_library':
            reads_alignment_object_ref = self._process_single_reads_library(input_object_info, 
                                                                            genome_index_base,
                                                                            result_directory,
                                                                            params)
        elif input_object_info['run_mode'] == 'sample_set':
            reads_alignment_object_ref = self._process_set_reads_library()

        report_output = self._generate_report(reads_alignment_object_ref,
                                              result_directory,
                                              params.get('workspace_name'))

        returnVal = {'result_directory': result_directory,
                     'reads_alignment_object_ref': reads_alignment_object_ref}

        returnVal.update(report_output)

        return returnVal
