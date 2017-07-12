import time
import json
import os
import uuid
import errno
import subprocess
import re
import zipfile
from pathos.multiprocessing import ProcessingPool as Pool
import multiprocessing

from Workspace.WorkspaceClient import Workspace as Workspace
from kb_Bowtie2.kb_Bowtie2Client import kb_Bowtie2
from ReadsAlignmentUtils.ReadsAlignmentUtilsClient import ReadsAlignmentUtils
from ReadsUtils.ReadsUtilsClient import ReadsUtils
from KBaseReport.KBaseReportClient import KBaseReport
from SetAPI.SetAPIServiceClient import SetAPI
from DataFileUtil.DataFileUtilClient import DataFileUtil


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
                  'alignment_suffix']:
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

    def _generate_alignment_set_data(self, reads_alignment_object_refs, sampleset_id):
        """
        _generate_alignment_set_data: generate AlignmentSet data from alignment objects
        """
        alignment_set_data = {}
        mapped_alignments_ids = []
        mapped_rnaseq_alignments = []
        read_sample_ids = []
        sample_alignments = []

        for reads_alignment_object_ref in reads_alignment_object_refs:
            alignment_data = self.ws.get_objects2({'objects': 
                                                  [{'ref':
                                                   reads_alignment_object_ref}]})['data'][0]['data']
            reads_ref = alignment_data.get('read_sample_id')
            read_sample_ids.append(reads_ref)
            sample_alignments.append(reads_alignment_object_ref)

            alignment_name = self.ws.get_object_info3({'objects': 
                                                      [{'ref': 
                                                       reads_alignment_object_ref}]})['infos'][0][1]

            reads_name = self.ws.get_object_info3({'objects': [{'ref': reads_ref}]})['infos'][0][1]

            mapped_alignments_ids.append({reads_ref: reads_alignment_object_ref})
            mapped_rnaseq_alignments.append({reads_name: alignment_name})

            genome_id = alignment_data.get('genome_id')

            if not alignment_set_data.get('genome_id'):
                alignment_set_data.update({'genome_id': genome_id})

        alignment_set_data.update({'mapped_alignments_ids': mapped_alignments_ids})
        alignment_set_data.update({'mapped_rnaseq_alignments': mapped_rnaseq_alignments})
        alignment_set_data.update({'read_sample_ids': read_sample_ids})
        alignment_set_data.update({'sample_alignments': sample_alignments})
        alignment_set_data.update({'sampleset_id': sampleset_id})

        return alignment_set_data

    def _save_alignment_set(self, reads_alignment_object_refs, workspace_name, alignment_set_name,
                            sampleset_id):
        """
        _save_alignment_set: upload AlignmentSet object
        """

        log('starting saving ReadsAlignmentSet object')

        if isinstance(workspace_name, int) or workspace_name.isdigit():
            workspace_id = workspace_name
        else:
            workspace_id = self.dfu.ws_name_to_id(workspace_name)

        object_type = 'KBaseRNASeq.RNASeqAlignmentSet'

        alignment_set_data = self._generate_alignment_set_data(reads_alignment_object_refs, 
                                                               sampleset_id)

        save_object_params = {
            'id': workspace_id,
            'objects': [{'type': object_type,
                         'data': alignment_set_data,
                         'name': alignment_set_name}]
        }

        dfu_oi = self.dfu.save_objects(save_object_params)[0]
        alignment_set_object_ref = str(dfu_oi[6]) + '/' + str(dfu_oi[0]) + '/' + str(dfu_oi[4])

        return alignment_set_object_ref

    def _process_single_reads_library(self, input_object_info, genome_index_base, 
                                      result_directory, cli_option_params):
        """
        _process_single_reads_library: process single reads library
        """

        reads_obj_type = self._get_type_from_obj_info(input_object_info['info'])
        reads_obj_name = input_object_info['info'][1]
        reads_files = self._get_reads_file(input_object_info['ref'], reads_obj_type, result_directory)
        
        tophat_result_dir = os.path.join(result_directory, 
                                         'tophat2_result_' + reads_obj_name + 
                                         '_' + str(int(time.time() * 100)))
        command = self._generate_command(genome_index_base, reads_files, 
                                         tophat_result_dir, cli_option_params)
        self._run_command(command)

        alignment_object_name = reads_obj_name + cli_option_params.get('alignment_suffix')
        reads_alignment_object_ref = self._save_alignment(tophat_result_dir,
                                                          alignment_object_name,
                                                          input_object_info['ref'],
                                                          cli_option_params.get('assembly_or_genome_ref'),
                                                          cli_option_params.get('workspace_name'),
                                                          cli_option_params.get('reads_condition'))

        return reads_alignment_object_ref

    def _generate_report_single_library(self, reads_alignment_object_ref, result_directory, 
                                        workspace_name):
        """
        _generate_report_single_library: generate summary report for single library
        """

        log('start creating report')

        output_files = self._generate_output_file_list_single_library(result_directory)
        output_html_files = self._generate_html_report_single_library(reads_alignment_object_ref, 
                                                                      result_directory)

        report_params = {'message': '',
                         'workspace_name': workspace_name,
                         'file_links': output_files,
                         'objects_created': [{'ref': reads_alignment_object_ref,
                                              'description': 'RNASeq Alignment generated by TopHat2'}],
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
                            'description': 'RNASeq Alignment Set generated by TopHat2'}]
        alignment_set_data = self.ws.get_objects2({'objects': 
                                                  [{'ref':
                                                   reads_alignment_object_ref}]})['data'][0]
        alignment_refs = alignment_set_data['data'].get('sample_alignments')
        for alignment_ref in alignment_refs:
            objects_created.append({'ref': alignment_ref,
                                    'description': 'RNASeq Alignment generated by TopHat2'})

        output_files = self._generate_output_file_list_sets_library(result_directory)
        output_html_files = self._generate_html_report_sets_library(reads_alignment_object_ref, 
                                                                    result_directory)

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

    def _generate_html_report_sets_library(self, reads_alignment_object_ref, result_directory):
        """
        _generate_html_report_sets_library: generate html summary report
        """

        log('start generating html report')

        html_report = list()
        result_file_path = os.path.join(result_directory, 'report.html')

        alignment_set_data = self.ws.get_objects2({'objects': 
                                                  [{'ref':
                                                   reads_alignment_object_ref}]})['data'][0]['data']
        alignment_set_info = self.ws.get_objects2({'objects': 
                                                  [{'ref':
                                                   reads_alignment_object_ref}]})['data'][0]['info']
        sample_alignments = alignment_set_data['sample_alignments']

        Overview_Content = ''
        Overview_Content += '<br/><table><tr><th>Generated AlignmentSet Object</th></tr>'
        Overview_Content += '<tr><td>{} ({})</td></tr></table>'.format(alignment_set_info[1],
                                                                       reads_alignment_object_ref)
        Overview_Content += '<p><br/></p>'
        Overview_Content += '<table><tr><th>Generated Alignment Objects</th><th></th><th></th>'
        Overview_Content += '<th></th><th></th><th></th><th></th></tr>'
        Overview_Content += '<tr><th>Alignment Name</th><th>Condition</th><th>Total Reads</th>'
        Overview_Content += '<th>Mapped Reads</th><th>Unmapped Reads</th><th>Singletons</th>'
        Overview_Content += '<th>Alignment Rate</th></tr>'
        
        for sample_alignment in sample_alignments:
            alignment_data = self.ws.get_objects2({'objects': 
                                                  [{'ref':
                                                   sample_alignment}]})['data'][0]['data']

            alignment_info = self.ws.get_objects2({'objects': 
                                                  [{'ref':
                                                   sample_alignment}]})['data'][0]['info']
            alignment_stats = alignment_data['alignment_stats']

            Overview_Content += '<tr><td>{} ({})</td><td>{}</td>'.format(alignment_info[1],
                                                                         sample_alignment,
                                                                         alignment_data['condition'])
            Overview_Content += '<td>{}</td><td>{}</td>'.format(alignment_stats['total_reads'],
                                                                alignment_stats['mapped_reads'])
            Overview_Content += '<td>{}</td><td>{}</td>'.format(alignment_stats['unmapped_reads'],
                                                                alignment_stats['singletons'])
            Overview_Content += '<td>{}%</td></tr>'.format(alignment_stats['alignment_rate'])
        Overview_Content += '</table>'
        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for TopHat2 App'})
        return html_report

    def _generate_html_report_single_library(self, reads_alignment_object_ref, result_directory):
        """
        _generate_html_report_single_library: generate html summary report
        """

        log('start generating html report')

        html_report = list()
        result_file_path = os.path.join(result_directory, 'report.html')

        alignment_data = self.ws.get_objects2({'objects': 
                                              [{'ref':
                                                reads_alignment_object_ref}]})['data'][0]['data']
        alignment_info = self.ws.get_objects2({'objects': 
                                              [{'ref':
                                                reads_alignment_object_ref}]})['data'][0]['info']
        alignment_stats = alignment_data['alignment_stats']

        Overview_Content = ''
        Overview_Content += '<p>Generated Alignment Object:</p>'
        Overview_Content += '<p>{} ({})</p>'.format(alignment_info[1], reads_alignment_object_ref)
        Overview_Content += '<br />'
        Overview_Content += '<p>Library Type: {}</p>'.format(alignment_data['library_type'])
        Overview_Content += '<p>Condition: {}</p>'.format(alignment_data['condition'])
        Overview_Content += '<p>Total Reads: {}</p>'.format(alignment_stats['total_reads'])
        Overview_Content += '<p>Mapped Reads: {}</p>'.format(alignment_stats['mapped_reads'])
        Overview_Content += '<p>Unmapped Reads: {}</p>'.format(alignment_stats['unmapped_reads'])
        Overview_Content += '<p>Singletons: {}</p>'.format(alignment_stats['singletons'])
        Overview_Content += '<p>Alignment Rate: {}%</p>'.format(alignment_stats['alignment_rate'])

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                result_file.write(report_template)

        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for TopHat2 App'})
        return html_report

    def _generate_output_file_list_single_library(self, result_directory):
        """
        _generate_output_file_list_single_library: zip result files and generate file_links for report
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

    def _generate_output_file_list_sets_library(self, result_directory):
        """
        _generate_output_file_list_sets_library: zip result files and generate file_links for report
        """

        log('start packing result files')

        output_files = list()
        result_dirs = os.listdir(result_directory)
        tophat2_result_dir_names = filter(re.compile('tophat2_result_*').match, result_dirs)
        result_files = list()

        for tophat2_result_dir_name in tophat2_result_dir_names:
            tophat2_result_dir = os.path.join(result_directory, tophat2_result_dir_name)
            result_file = result_file = os.path.join(result_directory, 
                                                     '{}.zip'.format(tophat2_result_dir_name.split('tophat2_result_')[1].rsplit('_', 1)[0]))
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
        object. E.g., if ref is a ReadsSet, return a list of all PairedEndLibrary or SingleEndLibrary
        refs that are a member of that ReadsSet. This is returned as a list of dictionaries as follows:
        {
            "ref": reads object reference,
            "condition": condition string associated with that reads object
        }
        The only one required is "ref", all other keys may or may not be present, based on the reads
        object or object type in initial ref variable. E.g. a RNASeqSampleSet might have condition info
        for each reads object, but a single PairedEndLibrary may not have that info.
        If ref is already a Reads library, just returns a list with ref as a single element.
        """
        obj_type = self._get_type_from_obj_info(info)
        refs = list()
        if "KBaseSets.ReadsSet" in obj_type:
            print("Looking up reads references in ReadsSet object")
            set_client = SetAPI(self.srv_wiz_url)
            reads_set = set_client.get_reads_set_v1({'ref': ref,
                                                     'include_item_info': 0
                                                     })
            for reads in reads_set["data"]["items"]:
                refs.append({'ref': reads['ref'],
                             'condition': reads['label']
                             })
        elif "KBaseRNASeq.RNASeqSampleSet" in obj_type:
            print("Looking up reads references in RNASeqSampleSet object")
            sample_set = self.ws.get_objects2({"objects": [{"ref": ref}]})["data"][0]["data"]
            for i in range(len(sample_set["sample_ids"])):
                refs.append({'ref': sample_set["sample_ids"][i],
                             'condition': sample_set["condition"][i]
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
        arg_2 = [genome_index_base] * 4
        arg_3 = [result_directory] * 4
        arg_4 = []
        for reads_ref in reads_refs:
            reads_input_object_info = self._get_input_object_info(reads_ref['ref'])
            cli_option_params['reads_condition'] = reads_ref['condition']
            arg_1.append(reads_input_object_info)
            arg_4.append(cli_option_params)

        cpus = min(cli_option_params.get('num_threads'), multiprocessing.cpu_count())
        pool = Pool(ncpus=cpus)
        log('running _process_alignment_object with {} cpus'.format(cpus))
        reads_alignment_object_refs = pool.map(self._process_single_reads_library, 
                                               arg_1, arg_2, arg_3, arg_4)

        reads_alignment_set_object_ref = self._save_alignment_set(reads_alignment_object_refs,
                                                                  cli_option_params['workspace_name'],
                                                                  alignment_set_name,
                                                                  input_object_info['ref'])

        return reads_alignment_set_object_ref

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['KB_AUTH_TOKEN']
        self.shock_url = config['shock-url']
        self.scratch = config['scratch']
        self.srv_wiz_url = config['srv-wiz-url']
        self.ws = Workspace(self.ws_url, token=self.token)

        self.bt = kb_Bowtie2(self.callback_url, service_ver='dev')
        self.rau = ReadsAlignmentUtils(self.callback_url, service_ver='dev')
        self.ru = ReadsUtils(self.callback_url)
        self.dfu = DataFileUtil(self.callback_url)

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
