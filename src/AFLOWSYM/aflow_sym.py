import json
import subprocess
import os


class Symmetry:

    def __init__(self, aflow_executable='aflow'):
        self.aflow_executable = aflow_executable


    def aflow_commmand(self, cmd):
        try:
            return subprocess.check_output(
                self.aflow_executable + cmd,
                shell=True
            )
        except subprocess.CalledProcessError:
            print(("Error aflow executable not found at: " + self.aflow_executable))
    
    
    def get_symmetry(self, input_file, tol=None):
        fpath = os.path.realpath(input_file.name)
        output = ""

        
        if tol is None:
            output = self.aflow_commmand(
                ' --aflowSYM --print=json --screen_only' + ' < ' + fpath
            )
        else:
            output = self.aflow_commmand(
                ' --aflowSYM=' + str(tol) + ' --print=json --screen_only' + ' < ' + fpath
            )
        res_json = json.loads(output)
        return res_json


    def get_edata(self, input_file, tol=None):
        fpath = os.path.realpath(input_file.name)
        output = ""
        
        if tol is None:
            output = self.aflow_commmand(
                ' --edata --print=json' + ' < ' + fpath
            )
        else:
            output = self.aflow_commmand(
                ' --edata=' + str(tol) + ' --print=json' + ' < ' + fpath
            )
        res_json = json.loads(output)
        return res_json


    def get_sgdata(self, input_file, tol=None):
        fpath = os.path.realpath(input_file.name)
        output = ""
        
        if tol is None:
            output = self.aflow_commmand(
                ' --sgdata --print=json' + ' < ' + fpath
            )
        else:
            output = self.aflow_commmand(
                ' --sgdata=' + str(tol) + ' --print=json' + ' < ' + fpath
            )
        res_json = json.loads(output)
        return res_json

