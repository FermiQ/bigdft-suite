"""
This file contains some low-level useful functions
"""

from __future__ import print_function

def find_files(regexp, archive=None):
    """
    Returns a list of all the paths to files that follow the regular expression
    regexp. They are searched from the current working directory or from an archive 
    given as optional argument.


    :param regexp: A regular expression
    :type regexp: string
    :param archive: an opened tarfile archive (optional)
    :type archive: 
    :returns: a list of all the paths that agree with the regexp
    :rtype: list of strings
    :raises: ValueError if the regexp does not find a single path.


    Example::

        #Find all python files in the current working directory
        find_files('*py') 

        #An exmple outside of the current working directory
        find_files('*/log-*.yaml') 

        #Example using a tarfile
        import tarfile
        my_archive = tarfile,open('archive.tar.gz')
        find_files('*/*/log-*.yaml', archive=my_archive)
    """

    if archive:
        return _find_files_from_archive(regexp, archive) 
    else:
        return _find_files_from_current_folder(regexp)


def _find_files_from_archive(re, archive):
    """
    This function retrieves the list of Logfiles instances 
    from the file archived satisfying a regular expression.   
    #function to identify an archive out of its regexp, 
    #solves the bug in re for '*' (solved in Python 2.7.6)
    """
    #Maybe some
    #from ../../../../bigdft/src/python/BigDFT import Logfiles as lf 
    #from futile import Yaml
    result=[]
    for f in archive.getnames():
        isthere=True
        for pattern in re.split('*'):
            isthere = isthere and (pattern in f)
        if isthere: 
            #print(f)
            filetar = archive.extractfile(f)
            content = filetar.read()
            dictionary = Yaml.yaml.load(stream=content)
            result.append(lf.Logfile(dictionary=dictionary))
    return result


def _find_files_from_current_folder(regexp):
    """
    Returns a list of all the files paths that follow the regular expression
    regexp in the current working directory.
    """
    import os
    paths = os.popen('ls '+regexp).read().splitlines()
    #print(type(paths), paths)
    if paths == []:
        raise ValueError('The regexp "{}" leads to no file. '\
                         'Consider using another one.'.format(regexp))
    else:
        return paths




if __name__ == '__main__':
    #Tests of the find_files function
    print(find_files("*py"))
    print(find_files("Utils.py"))
    try:
        find_files('*html')
    except ValueError as e:
        print('This raised the following ValueError:')
        print(e)
    try:
        find_files('*py', archive='file.tar.gz')
    except Exception as e:
        #print(dir(e))
        print('This raised the following Exception:')
        print(e)
    import os
    fname = 'file.tar.gz'
    os.system('find * -name "*py" | tar -zcvf '+fname+' -T -')
    os.system('ls '+fname)
    import tarfile
    my_archive = tarfile.open(fname)
    ##not working at present, need to import Logfile and Yaml
    ##in _find_files_from_archive
    #find_files('*py', archive=my_archive)
    os.system('rm '+fname)
