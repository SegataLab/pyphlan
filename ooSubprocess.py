#!/usr/bin/env python
# Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy


import subprocess
import os
import multiprocessing
from multiprocessing.pool import ThreadPool
import sys
import cStringIO


class ooSubprocess:

    def __init__(self, tmp_dir='tmp/'):
        self.chain_cmds = []
        self.tmp_dir = tmp_dir
        mkdir(tmp_dir)

    def ex(
            self,
            prog,
            args=[],
            get_output=False,
            get_out_pipe=False,
            out_fn=None,
            in_pipe=None,
            error_pipe=None,
            current_dir=None,
            verbose=True):

        if isinstance(args, str):
            args = args.split()

        if not isinstance(args, list):
            args = [args]

        cmd = [prog] + args
        print_cmd = 'ooSubprocess: ' + ' '.join(cmd)
        if verbose and out_fn and (not get_output):
            print_stderr(print_cmd + ' > ' + out_fn)
        elif verbose:
            print_stderr(print_cmd)

        if get_output:
            result = subprocess.check_output(
                                                cmd,
                                                cwd=current_dir,
                                                stdin=in_pipe,
                                                stderr=error_pipe)
        elif get_out_pipe:
            result = subprocess.Popen(
                                        cmd, 
                                        stdin=in_pipe, 
                                        stderr=error_pipe,
                                        stdout=subprocess.PIPE, 
                                        cwd=current_dir).stdout
        elif out_fn:
            ofile = open(out_fn, 'w') if out_fn else None
            result = subprocess.check_call(
                                            cmd, 
                                            stdin=in_pipe, 
                                            stderr=error_pipe,
                                            stdout=ofile, 
                                            cwd=current_dir)
            ofile.close()
        else:
            result = subprocess.check_call(
                                            cmd, 
                                            stdin=in_pipe, 
                                            stderr=error_pipe,
                                            cwd=current_dir)
        return result

    def chain(
            self,
            prog,
            args=[],
            stop=False,
            in_pipe=None,
            error_pipe=None,
            get_output=False,
            get_out_pipe=False,
            out_fn=None,
            current_dir=None,
            verbose=True):

        if in_pipe is None and self.chain_cmds != []:
            print_stderr(
                'Error: the pipeline was not stopped before creating a new one!')
            print_stderr('In cach: %s' % (' | '.join(self.chain_cmds)))
            exit(1)
        if out_fn and stop == False:
            print_stderr(
                'Error: out_fn (output_file_name) is only specified when stop = True!')
            exit(1)

        if isinstance(args, str):
            args = args.split()

        if not isinstance(args, list):
            args = [args]
        cmd = [prog] + args

        print_cmd = ' '.join(cmd)
        if out_fn and (not get_output):
            print_cmd += ' > ' + out_fn
        self.chain_cmds.append(print_cmd)

        if stop:
            if in_pipe is None:
                print_stderr('Error: No input process to create a pipeline!')
                exit(1)

            if verbose:
                print_stderr('ooSubprocess: ' + ' | '.join(self.chain_cmds))

            self.chain_cmds = []
            if get_output:
                result = subprocess.check_output(
                                                    cmd,
                                                    stdin=in_pipe,
                                                    stderr=error_pipe,
                                                    cwd=current_dir)
            elif get_out_pipe:
                result = subprocess.Popen(
                                            cmd,
                                            stdin=in_pipe,
                                            stderr=error_pipe,
                                            stdout=subprocess.PIPE,
                                            cwd=current_dir).stdout
            elif out_fn:
                ofile = open(out_fn, 'w')
                result = subprocess.check_call(
                                                cmd,
                                                stdin=in_pipe,
                                                stderr=error_pipe,
                                                stdout=ofile,
                                                cwd=current_dir)
                ofile.close()
            else:
                result = subprocess.check_call(
                                                cmd,
                                                stdin=in_pipe,
                                                stderr=error_pipe,
                                                cwd=current_dir)
        else:
            result = subprocess.Popen(
                                        cmd,
                                        stdin=in_pipe,
                                        stderr=error_pipe,
                                        stdout=subprocess.PIPE,
                                        cwd=current_dir).stdout
        return result

    def ftmp(self, ifn):
        return os.path.join(self.tmp_dir, os.path.basename(ifn))


def fdir(dir, ifn):
    return os.path.join(dir, os.path.basename(ifn))


def mkdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
    elif not os.path.isdir(dir):
        print_stderr('Error: %s is not a directory!' % dir)
        exit(1)


def replace_ext(ifn, old_ext, new_ext):
    # if not os.path.isfile(ifn):
    #    print_stderr('Error: file %s does not exist!'%(ifn))
    #    exit(1)
    if ifn[len(ifn) - len(old_ext):] != old_ext:
        # print_stderr('Error: the old file extension %s does not match!'%old_ext)
        # exit(1)
        new_ifn = ifn + new_ext
    else:
        new_ifn = ifn[:len(ifn) - len(old_ext)] + new_ext
    return new_ifn


def splitext(ifn):
    base, ext = os.path.splitext(ifn)
    base1, ext1 = os.path.splitext(base)
    if ext1 in ['.tar', '.fastq', '.fasta', '.sam']:
        base = base1
        ext = ext1 + ext
    return base, ext


def splitext2(ifn):
    basename = os.path.basename(ifn)
    base, ext = splitext(basename)
    return base, ext


def parallelize(func, args, nprocs=1, use_threads=False):
    if nprocs > 1:
        if use_threads:
            pool = ThreadPool(nprocs)
        else:
            pool = multiprocessing.Pool(nprocs)
        results = pool.map(func, args)
        pool.close()
        pool.join()
    else:
        results = serialize(func, args)
    return results


def serialize(func, args):
    results = []
    for arg in args:
        results.append(func(arg))
    return results


def print_stderr(*args):
    sys.stderr.write(' '.join(map(str, args)) + '\n')
    sys.stderr.flush()


def print_stdout(*args):
    sys.stdout.write(' '.join(map(str, args)) + '\n')
    sys.stdout.flush()
