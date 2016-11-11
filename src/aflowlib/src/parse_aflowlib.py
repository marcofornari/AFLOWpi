import AFLOWpi
import urllib2
import ast
import copy 
import re
class parser():



    def __init__(self):
        self.search_base = 'http://aflowlib.duke.edu/search/API/?'
        self.search_parameters = {}
        self.auid_list         = []
        self.search_vals       = []
        self.results           = {}
#    def by_icsd(self,icsd_number,attribute=None):
        #    connection = urllib2.urlopen('http://aflowlib.org/material.php?id=%s'%icsd_number)                        


    def add_value(self,parameter):
        self.search_vals.append(parameter)


    def remove_value(self,parameter):
        pass

    def remove_condition(self,parameter):
        condition_string=''.join(condition_string.split())
        self.search_parameters.remove(condition_string)

    def add_condition(self,parameter,condition):
#        condition_string=''.join(condition_string.split()

        self.search_parameters[parameter]=condition
    def get_file(self,filename):
#        relax_file_str = 
        for auid,entry in self.results.iteritems():
            search_str = 'http://'+entry['aurl'].replace(':','/')+'/'

            try:
                connection = urllib2.urlopen(search_str+filename)
                relax_file_str = connection.read()
                
                page_str = '!# http://aflowlib.org/material.php?id=aflow:'+auid+'\n'
                   

                self.results[auid][filename]=page_str+relax_file_str
            except:
                print '%s could not be extracted from AFLOWlib for entry %s'%(filename,auid)
                self.results[auid][filename]=''



    
#    def _transform_condition(parameter,condition):

        

    def search(self,):
        
        search_string=''
        for parameter,condition in self.search_parameters.iteritems():
            search_string+=parameter+'('+condition+')'+','
        for value in self.search_vals:
            search_string+=value+','
#            search_string+=self._transform_condition(parameter,condition)+','
        #truncate the tail comma

        null_to_none   = re.compile('null')
        true_to_True   = re.compile('true')
        false_to_False = re.compile('false')

        paging=1
        return_dict={}
        # while True:
        #     search = self.search_base+search_string+'paging(%d)'%paging
        #     connection = urllib2.urlopen(search)            

        #     res_str = connection.read()

        #     #end of the paging
        #     if res_str.strip()=='[]':
        #         break


        #     #change some json stuff to python form before using ast.literal
        #     res_str = true_to_True.sub('True',res_str)
        #     res_str = false_to_False.sub('False',res_str)
        #     res_str = null_to_none.sub('None',res_str)
            
        #     temp_dict=ast.literal_eval(res_str).values()     

        #     for v in temp_dict:
        #         ID=v['auid'].split(':')[-1]
        #         v_copy = copy.deepcopy(v)
        #         del v_copy['auid']
        #         return_dict[ID]=v_copy




        search = self.search_base+search_string+'paging(%d)'%paging
        connection = urllib2.urlopen(search)            
        
        res_str = connection.read()
        
        #end of the paging        

        #change some json stuff to python form before using ast.literal
        res_str = true_to_True.sub('True',res_str)
        res_str = false_to_False.sub('False',res_str)
        res_str = null_to_none.sub('None',res_str)
        
        temp_dict=ast.literal_eval(res_str).values()     

        for v in temp_dict:
            ID=v['auid'].split(':')[-1]
            v_copy = copy.deepcopy(v)
            del v_copy['auid']
            return_dict[ID]=v_copy


        paging+=1
        
        print len(return_dict.keys())

        self.results = return_dict
        return return_dict
#        return connection

    def display(self,parameter_list=None):
        pass
