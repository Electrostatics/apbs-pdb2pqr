"""SCons.Tool.cpython

Tool-specific initialization for Python binary builder.

There normally shouldn't be any need to import this module directly.
It will usually be imported through the generic SCons.Tool.Tool()
selection method.

"""

#
# Copyright (c) 2001-7,2010 The SCons Foundation
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__revision__ = "src/engine/SCons/Tool/python.py 3263 2008/07/31 13:50:51 MatiGruca"

import os

import SCons.Action
import SCons.Builder
import SCons.Errors


##########################################################################
#  Create Python builder
def createPythonBuilder(env):
    """This is a utility function that creates the InstallPython
    Builder in an Environment if it's not there already.

    If it's there already, we return the existing one.

    This builder is based on Install/InstallAs methods. It makes use
    of those builder's functions: installFunc(), and
    add_targets_to_INSTALLED_FILES()).
    """

    try:
        PythonInstallBuilder = env['BUILDERS']['InstallPython']
    except KeyError:
        from SCons.Tool.install import installFunc, add_targets_to_INSTALLED_FILES
        installpython_action = SCons.Action.Action(installFunc, '$CPYTHON_PYCOMSTR')
        PythonInstallBuilder = SCons.Builder.Builder(
                                 action = installpython_action,
                                 src_suffix = '$CPYTHON_SUFFIX',
                                 target_factory = env.fs.Entry,
                                 source_factory = env.fs.Entry,
                                 multi = 1,
                                 emitter = [ add_targets_to_INSTALLED_FILES ],
                                 name = 'InstallPythonBuilder')

    return PythonInstallBuilder

def InstallPython(env, target=None, source=None, dir=None, **kw):
    """InstallPython creates .pyc or .pyo files for .py source files
    and adds them to the list of targets along with the source files.
    They are later copied to the destination (target) directory.

    InstallPython takes a target (destination) directory as its first
    argument and a list of source files/directories as a second argument.

    InstallPython returns the list of target files to copy to the
    target directory.
    """

    if target and dir:
        raise SCons.Errors.UserError, "Both target and dir defined for InstallPython(), only one may be defined."
    if not dir:
        dir=target
    
    try:
        dnodes = env.arg2nodes(dir, env.fs.Dir)
    except TypeError:
        raise SCons.Errors.UserError, "Target `%s' of Install() is a file, but should be a directory.  Perhaps you have the InstallPython() arguments backwards?" % str(dir)

    sources = env.arg2nodes(source, env.fs.Entry)
    tgt = []

    try:
        import py_compile
    except ImportError:
        raise SCons.Errors.InternalError, "Couldn't import py_compile module"

    # import `compileall` module only if there is a dir in sources list
    import SCons.Node
    dir_in_sources = [isinstance(i, SCons.Node.FS.Dir) for i in sources]
    if True in dir_in_sources:
        try:
            import compileall
        except ImportError:
            raise SCons.Errors.InternalError, "Couldn't import compileall module"
        import glob

    compile_pyc = True
    try:
        if int(env.subst('$CPYTHON_PYC')) == 0:
            compile_pyc = False
    except ValueError:
        pass

    if compile_pyc:
        CPYTHON_TARGETSUFFIX = 'c'
    else:
        CPYTHON_TARGETSUFFIX = 'o'

    PIB = PythonInstallBuilder

    for dnode in dnodes:
        for src in sources:
            # add *.py and *.pyc files from a directory to tgt list
            if isinstance(src, SCons.Node.FS.Dir) and compile_pyc:
                compileall.compile_dir(str(src), maxlevels = 0, quiet = 1)
                globpath = src.path + os.sep + '*.py'
                py_and_pycs = glob.glob(globpath) + glob.glob(globpath+'c')
                for filename in py_and_pycs:
                    target = env.fs.Entry('.'+os.sep+filename, dnode)
                    tgt.extend(apply(PIB, (env, target, filename), kw))
            # add *.py and *.pyo files from a directory to tgt list
            elif isinstance(src, SCons.Node.FS.Dir):
                to_compile = []
                py_files = glob.glob(src.path + os.sep + '*.py')
                for py_file in py_files:
                    to_compile.append(py_file)
                    target_path = '.' + os.sep + py_file

                    # add '.py' file to tgt list
                    py_src = env.fs.Entry(py_file)
                    py_tgt = env.fs.Entry(target_path, dnode)
                    tgt.extend(apply(PIB, (env, py_tgt, py_src), kw))

                    # add '.pyo' file to tgt list
                    pyo_src = env.fs.Entry(py_file + CPYTHON_TARGETSUFFIX)
                    pyo_tgt = env.fs.Entry(target_path + CPYTHON_TARGETSUFFIX, dnode)
                    tgt.extend(apply(PIB, (env, pyo_tgt, pyo_src), kw))
                act = SCons.Action.CommandAction('@$CPYTHON_PYCOM %s' % (' '.join(to_compile)))
                act([], [], env)
            # add single '.py' and '.pyc' or '.pyo' file to tgt list
            else:
                # add '.py' file to tgt list
                target = env.fs.Entry('.'+os.sep+src.path, dnode)
                tgt.extend(apply(PIB, (env, target, src), kw))

                # .pyc or .pyo source and target files
                pyco_src = env.fs.Entry(src.path + CPYTHON_TARGETSUFFIX)
                pyco_tgt = env.fs.Entry(target.path + CPYTHON_TARGETSUFFIX)

                if compile_pyc:
                    py_compile.compile(src.path)
                else:
                    act = SCons.Action.CommandAction('@$CPYTHON_PYCOM %s' % (src.path))
                    act([], [], env)

                # add '.pyc' or '.pyo' file to tgt list
                tgt.extend(apply(PIB, (env, pyco_tgt, pyco_src), kw))

    return tgt

def generate(env):
    from SCons.Tool.install import copyFunc

    try:
        env['INSTALL']
    except KeyError:
        env['INSTALL'] = copyFunc

    global PythonInstallBuilder
    PythonInstallBuilder = createPythonBuilder(env)

    env['CPYTHON_PYC'] = 1 # generate '.pyc' files by default

    env['CPYTHON_EXE'] = 'python'
    env['CPYTHON_PYO_FLAGS'] = '-O'
    env['CPYTHON_PYO_CMD'] = "-c 'import sys,py_compile; [py_compile.compile(i) for i in sys.argv[1:]]'"
    env['CPYTHON_PYCOM'] = '$CPYTHON_EXE $CPYTHON_PYO_FLAGS $CPYTHON_PYO_CMD'
    env['CPYTHON_PYCOMSTR'] = 'Install file: "$SOURCE" as "$TARGET"'

    env['CPYTHON_SUFFIX'] = '.py' # extension for Python source files
    
    try:
        env.AddMethod(InstallPython, "InstallPython")
    except AttributeError:
        # Looks like we use a pre-0.98 version of SCons...
        from SCons.Script.SConscript import SConsEnvironment
        SConsEnvironment.InstallPython = InstallPython

def exists(env):
    return 1
