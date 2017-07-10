# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests  # noqa: F401

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from biokbase.workspace.client import Workspace as workspaceService
from kb_tophat2.kb_tophat2Impl import kb_tophat2
from kb_tophat2.kb_tophat2Server import MethodContext
from kb_tophat2.authclient import KBaseAuth as _KBaseAuth
from kb_tophat2.Utils.TopHatUtil import TopHatUtil


class kb_tophat2Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_tophat2'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_tophat2',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = kb_tophat2(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_tophat2_" + str(suffix)
        cls.wsClient.create_workspace({'workspace': cls.wsName})

        cls.tophat_runner = TopHatUtil(cls.cfg)

        cls.prepare_data()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    @classmethod
    def prepare_data(cls):
        pass

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.__class__.wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_bad_run_tophat2_app_params(self):
        invalidate_input_params = {
            'missing_reads_ref': 'reads_ref',
            'bowtie_index': 'bowtie_index',
            'workspace_name': 'workspace_name',
            'alignment_object_name': 'alignment_object_name'
        }
        with self.assertRaisesRegexp(
                ValueError, '"reads_ref" parameter is required, but missing'):
            self.getImpl().run_tophat2_app(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
            'reads_ref': 'reads_ref',
            'missing_bowtie_index': 'bowtie_index',
            'workspace_name': 'workspace_name',
            'alignment_object_name': 'alignment_object_name'
        }
        with self.assertRaisesRegexp(
                ValueError, '"bowtie_index" parameter is required, but missing'):
            self.getImpl().run_tophat2_app(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
            'reads_ref': 'reads_ref',
            'bowtie_index': 'bowtie_index',
            'missing_workspace_name': 'workspace_name',
            'alignment_object_name': 'alignment_object_name'
        }
        with self.assertRaisesRegexp(
                ValueError, '"workspace_name" parameter is required, but missing'):
            self.getImpl().run_tophat2_app(self.getContext(), invalidate_input_params)

        invalidate_input_params = {
            'reads_ref': 'reads_ref',
            'bowtie_index': 'bowtie_index',
            'workspace_name': 'workspace_name',
            'missing_alignment_object_name': 'alignment_object_name'
        }
        with self.assertRaisesRegexp(
                ValueError, '"alignment_object_name" parameter is required, but missing'):
            self.getImpl().run_tophat2_app(self.getContext(), invalidate_input_params)
