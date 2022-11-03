#!/usr/bin/env python3
# cloud_workstation 2.0.0
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# See https://documentation.dnanexus.com/developer for documentation and
# tutorials on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import datetime
import sys

import dxpy

# in case PYTHONPATH is not set correctly
sys.path.append("/usr/lib/python3/")

import asset_builder_tools
import dx_utils

TIMEOUT_FILE = '/home/dnanexus/.dx.timeout'
TIME_FORMAT = '%Y %m %d %H %M %S'

SERVER_READY_TAG = 'server_running'


def _get_timeout():
    with open(TIMEOUT_FILE, 'r') as fh:
        timeout = fh.read().strip()
    return datetime.datetime.strptime(timeout, TIME_FORMAT)


@dxpy.entry_point('main')
def main(**job_inputs):
    asset_builder_tools.create_before_file_list()

    if 'snapshot' in job_inputs:
        cmd = 'dx cat {0} | sudo tar -C / -zxvf - '.format(job_inputs['snapshot']['$dnanexus_link'])
        dx_utils.run_cmd(cmd)
    if 'fids' in job_inputs:
        for fid in job_inputs['fids']:
            cmd = 'dx download {0}'.format(fid['$dnanexus_link'])
            dx_utils.run_cmd(cmd)

    # Loop until it's time to quit
    dx_utils.run_cmd('dx-set-timeout {0} '.format(job_inputs['max_session_length']))
    dx_utils.run_cmd('sudo chmod 666 {0} '.format(TIMEOUT_FILE))

    cmd = ['dx', 'tag', dxpy.JOB_ID, SERVER_READY_TAG]
    dx_utils.run_cmd(cmd)

    # interactive apps should run sub-jobs detached automatically.
    with open("/home/dnanexus/.bashrc", "a") as bashrc:
        bashrc.write("\nexport DX_RUN_DETACH=1\n")

    # Loop until it's time to quit
    while True:
        dx_utils.run_cmd(["sleep", "10"])
        timeout = _get_timeout()
        if datetime.datetime.now() > timeout:
            break

    output = {}

    return output


dxpy.run()
