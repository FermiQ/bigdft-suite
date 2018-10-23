"""
This file contains some low-level useful functions
"""

from __future__ import print_function

def write(*args,**kwargs):
    """
    Wrapper for print function or print to ensure compatibility with python 2
    The arguments are used similarly as the print_function
    They can also be generalized to python 2 cases
    """
    return print(*args,**kwargs)

def push_path(inp,*keys):
    """
    Follow in the dictionary inp the path indicated by the keys.
    If this path does not exists creates it.

    Args:
       inp (dict): dictionary
       keys (str): keys of the path to follow

    Returns:
       (``branch``,``key``) tuple, where

       * ``branch`` (dict): the dictionary of the second-last item of the path
       * ``key`` (str): the last item of the path
    """
    tmp=inp
    for i,key in enumerate(keys):
        k=key
        if i==len(keys)-1: break
        if key not in tmp: tmp[key]={}
        tmp=tmp[key]
    return tmp,k


def dict_set(inp,*subfields):
    """Ensure the provided fields and set the value

    Provide a entry point to the dictionary.
    Useful to define a key in a dictionary that may not have the
    previous keys already defined.

    Arguments:
       inp (dict): the top-level dictionary
       subfields (str,object): keys, ordered by level, that have to be retrieved from topmost level of ``inp``.
          The last item correspond to the value to be set .

    Example:

       >>> inp={}
       >>> dict_set(inp,'dft','nspin','mpol',2)
       >>> print (inp)
       {'dft': {'nspin': {'mpol': 2}}}

    """
    if len(subfields) <= 1:
        raise ValueError('invalid subfields')
    keys=subfields[:-1]
    tmp,key=push_path(inp,*keys)
    #tmp=inp
    #for i,key in enumerate(subfields):
    #    k=key
    #    if i==len(subfields)-2: break
    #    tmp[key]={}
    #    tmp=tmp[key]
    tmp[key]=subfields[-1]


def dict_merge(dest, src):
    """ Recursive dict merge. Inspired by :meth:`dict.update`, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``src`` is merged into
    ``dest``.  From :ref:`angstwad/dict-merge.py <https://gist.github.com/angstwad/bf22d1822c38a92ec0a9>`

    Arguments:
       dest (dict): dict onto which the merge is executed
       src (dict): dict merged into dest

    """
    import collections
    for k, v in src.iteritems():
        if (k in dest and isinstance(dest[k], dict)
                and isinstance(src[k], collections.Mapping)):
            dict_merge(dest[k], src[k])
        else:
            dest[k] = src[k]

def file_time(filename):
    """
    Gives the date of the creation of the file, if exists.

    :param str filename: name of the file
    :returns: if the file exists, the date of the filename as per os.path.getmtime.
     Otherwise it returns 0
    """
    import os
    if os.path.isfile(filename):
        return os.path.getmtime(filename)
    else:
        return 0

def make_dict(inp):
    """
    Transform the instance ``inp`` into a python dictionary. If inp is already a dictionary, it perfroms a copy.

    Args:
       inp (dict): a instance of a Class which inherits from dict

    Returns:
       dict: the copy of the class, converted as a dictionary
    """
    import copy
    local_tmp=copy.deepcopy(inp)
    local_input={}
    local_input.update(local_tmp)
    return local_input

def function_signature_regenerator(target_kwargs_function,fun_name='',fun_docstring='',**kwargs):
    '''
    Generate the function of the name provided by `fun_name`, with signature provided by the
    kwargs dictionary.

    Args:
       target_kwargs_function (func): keyword arguments function that will be used for the generated function.
       fun_name (str): name of the regenerated function. If empty it will be the ``target_kwargs_functon.__name__`` prefixed by ``regenerated``, which will be copied in the docstring of the regenerated function.
       fun_docstring (str): docstring of the generated function, if empty it will take the docstring from ``target_kwargs_function``.
       **kwargs: keyword arguments which will represent the signature of the generated function.

    Example:
        >>> def write_kwargs(**kwargs):
        >>>     """
        >>>     Convert keyword arguments into a string
        >>>     """
        >>>     return str(kwargs)
        >>> write_opts=function_signature_regenerator(write_kwargs,fun_name='write_opts',opt1='default1',opt2='default2')
        >>> help(write_opts)
        >>> print (write_opts())
        Help on function write_opts:

        write_opts(opt1='default1', opt2='default2')
              Convert keyword arguments into a string

        {'opt1': 'default1', 'opt2': 'default2'}

    '''
    signature=option_line_generator(',',**kwargs).lstrip(',')
    docstring=target_kwargs_function.__doc__ if not fun_docstring else fun_docstring
    if docstring is None: docstring="Automatically generated function from the target function '"+target_kwargs_function.__name__+"'"
    docstring='   """\n'+docstring+'\n   """'
    fname="regenerated_"+target_kwargs_function.__name__ if  not fun_name else fun_name
    function="def %s(%s):\n%s\n   return target_function(**locals())" % (fname,signature,docstring)
    gen_locals={}
    gen_object=compile(function,'generated_fun','exec')
    eval(gen_object,{'target_function':target_kwargs_function},gen_locals)
    return gen_locals[fname]

def option_line_generator(separator='--',**kwargs):
    """
    Associate to each of the keyword arguments a command line argument.

    Args:
       separator (str): The string needed to separate the options.
       Might be '--' for command-line arguments, but also ',' for function signatures.

    Warning:
        The separator comes **before** the first argument therefore pay attention to
        lstrip it in case you want to use it as a function signature string.

    Example:
        >>> option_line_generator(arg1='val1',arg2='val2')
        '--arg1=val1 --arg2=val2'
    """
    command=''
    for option,value in kwargs.items():
        command+=separator+option+'="'+str(value)+'" '
    return command

def option_line_generator(separator='--',**kwargs):
    """
    Associate to each of the keyword arguments a command line argument.

    Args:
       separator (str): The string needed to separate the options.
       Might be '--' for command-line arguments, but also ',' for function signatures.

    Warning:
        The separator comes **before** the first argument therefore pay attention to
        lstrip it in case you want to use it as a function signature string.

    Example:
        >>> option_line_generator(arg1='val1',arg2='val2')
        '--arg1=val1 --arg2=val2'
    """
    command=''
    for option,value in kwargs.items():
        command+=separator+option+'="'+str(value)+'" '
    return command

def kw_pop(*args,**kwargs):
    """
    Treatment of kwargs. Eliminate from kwargs the tuple in args.
    """
    arg=kwargs.copy()
    key,default=args
    if key in arg:
        return arg,arg.pop(key)
    else:
        return arg,default


def find_files(regexp, archive=None):
    """
    Returns a list of the paths to the files that follow the regular expression
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
        my_archive = tarfile.open('archive.tar.gz')
        find_files('*/*/log-*.yaml', archive=my_archive)
    """
    import os

    #Get a list of all paths to files satisfying the regexp
    if archive is not None:
        paths = _find_files_from_archive(regexp, archive)
    else:
        paths = os.popen('ls '+regexp).read().splitlines()

    #Test that the regexp found files
    if paths == []:
        raise ValueError('The regexp "{}" leads to no file. '\
                         'Consider using another one.'.format(regexp))
    else:
        return paths


def _find_files_from_archive(re, archive):
    """
    This function retrieves the list of Logfiles instances
    from the file archived satisfying a regular expression.
    #function to identify an archive out of its regexp,
    #solves the bug in re for '*' (solved in Python 2.7.6)
    """
    import tarfile

    #Open the archive
    with tarfile.open(archive, 'r') as arch:
    #Return paths to logfiles satisfying the regexp
        return [f for f in arch.getnames()
                if all(pattern in f for pattern in re.split('*'))]

def ensure_copy(src,dest):
    """Copy src into dest.

    Guarantees that the file indicated by ``dest`` is a copy of the file ``src``

    Args:
      src (str): path of the source file. Should be valid.
      dest (src): path of the destination file

    Returns:
      bool: ``True`` if the file needed to be copied, ``False`` if ``src`` and ``dest`` are identical
    """
    import shutil,os
    copied=False
    if (os.path.isfile(dest) and os.stat(dest) != os.stat(src)) or not os.path.isfile(dest):
        shutil.copy2(src,os.path.dirname(dest))
        copied=True
    return copied

def ensure_dir(file_path):
    """
    Guarantees the existance on the directory given by the (relative) file_path

    Args:
       file_path (str): path of thh directory to be created

    Returns:
       bool: True if the directory needed to be created, False if it existed already
    """
    import os
    directory = file_path
    created=False
    if not os.path.exists(directory):
        os.makedirs(directory)
        created=True
    return created

if __name__ == '__main__':
    import os

    #Tests of the find_files function
    #
    print("Test finding all python files in this directory")
    print(find_files("*py"))
    print()

    #
    print("Test finding the Utils.py file in this directory")
    print(find_files("Utils.py"))
    print()

    #
    print("Test raising a ValueError because the regexp leads to no files")
    try:
        find_files('*html')
    except ValueError as e:
        print('This raised the following ValueError:')
        print(e)
    print()

    #
    print("Test raising an exception because there is no such archive.")
    fname = 'file.tar.gz'
    if fname in os.popen('ls'): os.system('rm '+fname)
    #os.system('rm '+fname)
    try:
        find_files('*py', archive=fname)
    except Exception as e:
        #print(dir(e))
        print('This raised the following Exception:')
        print(e)
    print()

    #
    print("Test without error using an archive")
    os.system('find * -name "*py" | tar -zcvf '+fname+' -T -')
    os.system('ls '+fname)
    find_files('*py', archive=fname)
    os.system('rm '+fname)
