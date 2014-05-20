#!/usr/bin/env python
import subprocess
import os

class ooSubprocess:
	def __init__(self, a_tmp_dir = 'tmp/'):
		self.chain_cmds = []
		self.tmp_dir = a_tmp_dir
		self.mkdir(a_tmp_dir)
		
	def mkdir(self, dir):
		if not os.path.exists(dir):
			os.makedirs(dir)
		elif not os.path.isdir(dir):
			print 'Error: %s is not a directory!'%dir
			exit(1)

	def ex(self, prog, args = [], out_fn = None, get_output = False, a_cwd = None, verbose = True):
		if type(args) is str:
			args = args.split()
			
		if type(args) is not list:
			args = [args]

		cmd = [prog] + args
		print_cmd = ' '.join(cmd)
		if out_fn:
			if verbose:
				print 'ooSubprocess: ' + print_cmd + ' > ' + out_fn
			ofile = open(out_fn, 'w')
			if get_output:
				return subprocess.check_output(cmd, stdout = ofile, cwd = a_cwd)
			else:
				return subprocess.call(cmd, stdout = ofile, cwd = a_cwd)
		else:
			if verbose:
				print 'ooSubprocess: ' + print_cmd
			if get_output:
				return subprocess.check_output(cmd, cwd = a_cwd)
			else:
				return subprocess.check_call(cmd, cwd = a_cwd)


	def chain(self, prog, args = [], out_fn = None, stop = False, in_proc = None, get_output = False, a_cwd = None, verbose = True):
		if in_proc is None and self.chain_cmds != []:
			print 'Error: the pipeline was not stopped before creating a new one!'
			print 'In cach: %s'%(' | '.join(self.chain_cmds))
			exit(1)
		if out_fn and stop == False:
			print 'Error: out_fn (output_file_name) is only specified when stop = True!'
			exit(1)

		if type(args) is str:
			args = args.split()
			
		if type(args) is not list:
			args = [args]
		cmd = [prog] + args

		print_cmd = ' '.join(cmd)
		if out_fn: print_cmd += ' > ' + out_fn
		self.chain_cmds.append(print_cmd)

		in_pipe = in_proc.stdout if in_proc else None
		if stop:
			if in_pipe is None:
				print 'Error: No input process to create a pipeline!'
				exit(1)

			if verbose:
				print 'ooSubprocess: ' + ' | '.join(self.chain_cmds)
				
			self.chain_cmds = []
			ofile = open(out_fn, 'w') if out_fn else None
			if get_output:
				return subprocess.check_output(cmd, stdout = ofile, stdin = in_pipe, cwd = a_cwd)
			else:
				return subprocess.check_call(cmd, stdout = ofile, stdin = in_pipe, cwd = a_cwd)
		else:
			return subprocess.Popen(cmd, stdout = subprocess.PIPE, stdin = in_pipe, cwd = a_cwd)
	
	def change_ext(self, ifn, old_ext, new_ext):
		#if not os.path.isfile(ifn):
		#	print 'Error: file %s does not exist!'%(ifn)
		#	exit(1)
		if ifn[len(ifn) - len(old_ext):] != old_ext:
			#print 'Error: the old file extension %s does not match!'%old_ext
			#exit(1)
			new_ifn = ifn + new_ext
		else:
			new_ifn = ifn[:len(ifn) - len(old_ext)] + new_ext
		return new_ifn

	def ftmp(self, ifn):
		return os.path.join(self.tmp_dir, os.path.basename(ifn))
