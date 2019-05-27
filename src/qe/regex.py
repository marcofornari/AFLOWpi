import re
import AFLOWpi






def cell_parameters(string,return_which='content',regex_or_string='string'):
    modifier_regex = re.compile(r'(?:CELL_PARAMETERS)\s*(?:\s*[\[\(\{])*\s*([\w]*)\s*(?:[\]\)\}])*\s*\n*',re.MULTILINE)
    try:
        modifier = modifier_regex.findall(string)[-1]
    except:
        alatSearch = re.compile(r'(?:CELL_PARAMETERS)\s*(?:\s*[\[\(\{])*\s*(?:alat\s*=\s*([0-9.]*))\s*(?:[\]\)\}])*\s\*\n*')
        if len(alatSearch.findall(string))==0:
            modifier = ''
        else:
            modifier= '{bohr}'
    if return_which=='modifier':
        if regex_or_string=='regex':
            return modifier_regex
        else:
            return modifier

    if return_which=='content':
        content_return = re.compile(r'(?:CELL_PARAMETERS.*\n)(\s*(?:\s*(?:[-0-9.E+]+)\s+(?:[-0-9.E+]+)\s+(?:[-0-9.E+]+)\s*)+)',re.MULTILINE)    

    if regex_or_string=='string':
        try:
            content_return =  content_return.findall(string)[-1]
        except Exception,e: 
            return


    return content_return

def k_points(string,return_which='content',regex_or_string='string'):
    modifier_regex = re.compile(r'(?:K_POINTS)\s*(?:(?:\s*[\[\(\{])*\s*([A-Za-z_]+)\s*(?:[\]\)\}])*\s*\n*|(?:\s*))',re.MULTILINE)

    if return_which=='modifier':
        if regex_or_string=='regex':
            return modifier_regex
        else:
            try:
                returnVal = modifier_regex.findall(string)[-1]
                return returnVal
            except:
                return ''



    if return_which=='modifier':
        try:
            if regex_or_string=='regex':
                return modifier_regex
            else:

                return modifier

        except:
            return ''


    content_return = re.compile(r'(?:K_POINTS).*\n*((?:(?:\s*[A-Za-z]*\s*[-.0-9]+\n*)+))(?=^[\w]+|)')    
    modifier_regex = re.compile(r'(?:K_POINTS)\s*(?:(?:\s*[\[\(\{])*\s*([A-Za-z_]+)\s*(?:[\]\)\}])*\s*\n*|(?:\s*))',re.MULTILINE)
    try:
        modifier =  modifier_regex.findall(string)[-1]
        if modifier.strip().lower()=='crystal' or modifier.strip().lower()=='crystal_b':
            content_return = re.compile(r'(?:K_POINTS)(?:(?:\s*[\[\(\{]*)\s*(?:\w*)\s*(?:[\]\)\}]*))\s*\n([0-9\s]*\n\s*(?:\s*(?:[-0-9.]+)\s+(?:[-0-9.]+)\s+(?:[-0-9.]+)\s*[\s\d]*(?:[A-Za-z0-9!\s]*)\s*\n*)+)(?=(?:[A-Z|_|\s]+\n)|)')

    except:
        pass

    
#    content_return = re.compile(r'(?:K_POINTS).*\n*((?:(?:\s*[A-Za-z]*\s*[-.0-9]+\n*)+))(?=^[\w]+|)')    

#    content_return = re.compile(r'(?:K_POINTS).*\n*((?:[(0-9.\s]+\n*)+)(?=^[\w]+|)',re.MULTILINE)    
#    if return_which=='content':
#        if modifier=='automatic' or modifier=='' or modifier==None:

    if regex_or_string=='string':
        try:
            content_return =  content_return.findall(string)[-1]
        except Exception,e: 
            return ''
    return content_return





def atomic_positions(string,return_which='content',regex_or_string='string'):

    modifier_regex = re.compile(r'(?:ATOMIC_POSITIONS)\s*(?:\s*[\[\(\{])*\s*(\w*)\s*(?:[\]\)\}])*\s*\n*',re.MULTILINE)
    if return_which=='modifier':
        if regex_or_string=='regex':
            return modifier_regex
        else:
            try:
                returnVal = modifier_regex.findall(string)[-1]
                return returnVal
            except:
                return ''



    if return_which=='content':
        content_return = re.compile(r'(?:ATOMIC_POSITIONS)(?:(?:\s*[\[\(\{]*)\s*(?:\w*)\s*(?:[\]\)\}]*))\s*\n(\s*(?:(?:[A-Za-z0-9]+)\s+(?:[-0-9.]+)\s+(?:[-0-9.]+)\s+(?:[-0-9.]+)\s*(?:[-0-9.]+)\s*[\s\d]*\n*)+)(?=(?:[A-Z|_|\s]+\n)|)',re.MULTILINE)    
        if regex_or_string=='regex':
            return content_return

        try:
            modifier = modifier_regex.findall(string)[-1]
        except Exception,e:
            return ''




    if regex_or_string=='string':
        try:
            content_return =  content_return.findall(string)[-1]
        except Exception,e: 
            return
    return content_return


    

    
