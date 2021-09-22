#!/usr/bin/python3

"""
// This is a standalone QoMMMa Python file, which allows for automatic generation of a documentation file in .txt or .html format. //
"""

# Global imports.
import glob, os
import sys
import re

# Local imports.
import gauutil as gauutil
import jagutil as jagutil
import molutil as molutil
import orcautil as orcautil
import qomfreq as qomfreq
import qommma as qommma
import qomutil as qomutil
import qomdoc as qomdoc

def all_files(pattern, search_path):
    """
    
    // Function which searches for all available files matching a certain pattern in a given directory. //
    
    Arguments
    ----------
    pattern : string
         This is '*.py' as it is representative of all Python files.
    search_path : string
        The directory which is to be searched for files matching the defined pattern.

    """
    
    pathsep = os.pathsep
    for path in search_path.split(pathsep):
        for match in glob.glob(os.path.join(path, pattern)):
            yield match

def docwrite(py_list, mode, fdoc):
    """
    
    // Function which generates documentation for all functions contained within Python files. //
    
    Arguments
    ----------
    py_list : list
        A list of all the Python files in the working directory for which documentation is to be written.
    mode : string
        If the output file is desired simply as text, then 'text', and for html, it is 'html'.
    fn string : string
        The name of the output documentation file.

    """
  
    # Depending on whether the mode set to text or html file, the formatting of headings and main-body text is changed.
    if mode == 'text':
        h1 = '%s\n\n'
        h2 = '%s\n'
        h3 = '%s'
        doc = '%s\n\n'
    elif mode == 'html':
        h1 = '<h1>%s</h1>\n'
        h2 = '<h2>%s</h2>\n'
        h3 = '<h3>%s</h3>\n'
        doc = '<pre>%s</pre>\n'
    else:
        print('Unknown format ' + mode)
        sys.exit()
    
    for py_file in py_list:
        # The given Python file is opened, the documentation file is created.
        f_in = open(py_file, 'r')
        f_doc = open(fn, 'a')
        
        # Using a regex search, all the functions and their arguments defined in the given Python file are added to the lists function_list and args_list.
        function_list= []
        args_list = []
        for line in f_in:
            line = line.rstrip()
            if re.search('def' + ' ', line):
                temp_func = line.split("(")[0]
                temp_args = line.split("(")[1]
                args = '(' + temp_args
                func = temp_func.split(" ")
                function_list.append(func[1])
                args_list.append(args)
    
        # The title of the Python file is written to the documentation file.
        title = py_file.rstrip('.py')
        f_doc.write(h2%title)
        
        # Using a series of if statements, the docstring for the given Python file are written to the documentation file.
        # If new Python files are added to QoMMMa, then they will need to be added here since it is not possible to dynamically allocate imports to __doc__.
        if py_file == 'gauutil.py':
            docstring = gauutil.__doc__
        elif py_file == 'jagutil.py':
            docstring = jagutil.__doc__  
        elif py_file == 'molutil.py':
            docstring = molutil.__doc__
        elif py_file == 'orcautil.py':
            docstring = orcautil.__doc__  
        elif py_file == 'qomfreq.py':
            docstring = qomfreq.__doc__   
        elif py_file == 'qomutil.py':
            docstring = qomutil.__doc__                 
        elif py_file == 'qommma.py':
            docstring = qommma.__doc__
        elif py_file == 'qomdoc.py':
            docstring = qomdoc.__doc__
        else:
            pass
        f_doc.write(doc%docstring) 
        f_doc.write(doc%"// Within this file, the following functions are defined... //")
        
        # Using a series of if statements, the functions name and the docstrings relating to them for the given Python file are written to the documentation file.
        # Again, if new Python files are added to QoMMMa, then they will need to be added here since it is not possible to dynamically allocate imports to getattr.
        for num in range(len(function_list)):
            if py_file == 'gauutil.py':
                docstring = getattr(gauutil, function_list[num]).__doc__
            elif py_file == 'jagutil.py':
                docstring = getattr(jagutil, function_list[num]).__doc__  
            elif py_file == 'molutil.py':
                docstring = getattr(molutil, function_list[num]).__doc__
            elif py_file == 'orcautil.py':
                docstring = getattr(orcautil, function_list[num]).__doc__  
            elif py_file == 'qomfreq.py':
                docstring = getattr(qomfreq, function_list[num]).__doc__   
            elif py_file == 'qomutil.py':
                docstring = getattr(qomutil, function_list[num]).__doc__                 
            elif py_file == 'qommma.py':
                docstring = getattr(qommma, function_list[num]).__doc__
            elif py_file == 'qomdoc.py':
                docstring = getattr(qomdoc, function_list[num]).__doc__
            else:
                pass
            f_doc.write(h3%(function_list[num] + args_list[num]))
            f_doc.write(doc%docstring)  
        f_in.close()
    f_doc.close()

if __name__ == "__main__":
    # The mode is read from the second command line argument.
    try:
        mode = sys.argv[1]
    except:
        print('Usage: qomdoc.py [text|html]')
      
    # Depending on whether the mode set to text or html file, the name of the output file is set.
    if mode == 'html':
        fn = 'qomdoc.html'
    elif mode == 'text':
        fn = 'qomdoc.txt'
    else:
        print('Unknown format ' + mode)
        sys.exit()
    
    # If the output file already exists, then it is removed so the new one can be made.
    if os.path.exists(fn):
        os.remove(fn)
    
    # There is specific header and footer formatting required for a html file, so these variables are initialised.
    fd = open(fn, 'w')
    if mode == 'html':
        header = '''<!DOCTYPE html
    PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
    "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
    <html>
     <head>
      <title>QoMMMa programmers\' manual</title>
      <link rel="stylesheet" type="text/css" href="qommma.css" />
      <meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1" />
    <h1>QoMMMa programmers\' manual</h1>
    <h3>// QoMMMa is currently in version 8.07 and is written in Python 3. //<h3>
    <h3>// The Python based QoMMMa files are designed to couple and run QM, MM, and QoMMMa Fortran codes. //<h3>
    <h3>// All Python files are written following (largely) PEP 8 style guidelines. //<h3>
    '''
        footer = '''</body>
    </html>
    '''
        fd.write(header)
    fd.close()
    
    # The name of every Python file is added to the list py_list.
    # This list is used in docwrite to create the output documentation file containing docstrings from all python files in the working directory.
    py_list = []
    for match in all_files('*.py', os.getcwd()):
        fname = os.path.split(match)[1]
        py_list.append(fname)
    docwrite(py_list, mode, fd) 
    
    # For formatting, the footer of the html file is written.
    fd = open(fn, 'a')
    if mode == 'html':
        fd.write(footer)
    fd.close()