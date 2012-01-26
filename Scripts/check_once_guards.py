#!/usr/bin/env python
import os
import re
import sys
from contextlib import contextmanager
 
@contextmanager
def WorkingDirectory(dst_dir):
    current_dir = os.getcwd()
    os.chdir(dst_dir)
    try:
        yield dst_dir
    finally:
        os.chdir(current_dir)

includeFinder = re.compile(r'\A#include\s+"(.+?)"')
def FileNonSystemIncludedFileGenerator(filename):
    """Given a filename, return a generator that will yield the line
    numbers and paths included by non-system include statements in the
    file.
    """
    # Want to match: start of line; '#include'; one or more space;
    # double quote; anything (but non-greedily); double quote.
    
    codefile = file(filename)
    for iLine, line in enumerate(codefile):
        match = includeFinder.match(line)
        if match:
            yield (iLine + 1), match.group(1)


extensions = set(('.cc', '.hpp', '.h'))
def Walk(codeDir):
    for dirpath, dirnames, filenames in os.walk(codeDir):
        if 'build' in dirnames:
            dirnames.remove('build')
            pass
        
        for name in filenames:
            base, ext = os.path.splitext(name)
            if ext in extensions or (ext == '.in' and os.path.splitext(base)[1] in extensions):
                yield os.path.normpath(os.path.join(dirpath, name))

def CheckIncludePaths(sourceFile):
    errors = False
    for lineNumber, include in FileNonSystemIncludedFileGenerator(sourceFile):
        if include in ignoredIncludes:
            continue
        
        if not os.path.exists(include) and not os.path.exists(include + '.in'):
            sys.stderr.write(
                '{file}:{line} Bad include path "{dodgy}"\n'.format(file=sourceFile,
                                                                    line=lineNumber,
                                                                    dodgy=include)
                )
            errors = True
    return errors

def GetGuardLines(filename):
    f = file(filename)
    lines = [f.readline(), f.readline()]
    
    for line in f:
        continue
    lines.append(line)
    return lines

def SplitAll(path):
    head, tail = os.path.split(path)
    if head == '':
        ans = []
    else:
        ans = SplitAll(head)
        
    ans.append(tail)
    return ans

def CheckGuardErrors(sourceFile):
    lines= GetGuardLines(sourceFile)
    
    parts = SplitAll(sourceFile)
    parts[-1] = parts[-1].replace('.', '_')
    define = ('HEMELB_' + '_'.join(parts)).upper()

    error = False
    
    line0 = '#ifndef {define}\n'.format(define=define)
    if lines[0] != line0:
        sys.stderr.write(
                '{file}:1 Bad include guard; must be {required!r}\n'.format(file=sourceFile,
                                                                            required=line0)
            )
        error = True

    line1 = '#define {define}\n'.format(define=define)
    if lines[1] != line1:
        sys.stderr.write(
                '{file}:2 Bad include guard; must be {required!r}\n'.format(file=sourceFile,
                                                                            required=line1)
            )
        error = True
        
    return error

if __name__ == '__main__':
    script = os.path.abspath(sys.argv[0])
    scriptDir = os.path.dirname(script)
    codeDir = os.path.normpath(os.path.join(scriptDir, os.pardir, 'Code'))

    ignoredIncludes = set((
        'tinyxml.h',
        'parmetis.h',
        'ctemplate/template.h',
        ))

    errors = False
    with WorkingDirectory(codeDir):
        
        for sourceFile in Walk('.'):
            includeErrors = CheckIncludePaths(sourceFile)
            
            base, ext = os.path.splitext(sourceFile)
            if ext == '.h' or ext == '.hpp' or ext == '.in':
                guardErrors = CheckGuardErrors(sourceFile)
            else:
                guardErrors = False
            errors = errors or (includeErrors or guardErrors)
            
    if errors:
        raise SystemExit(1)
    
        
