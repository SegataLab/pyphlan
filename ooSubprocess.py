#!/usr/bin/env python
import subprocess

class ooSubprocess:
	def __init__(self):
		self.chain_cmds = []

	def ex(self, prog, args = [], out_fn = None, verbose = True):
		if type(args) is not list:
			args = [args]
		cmd = [prog] + args
		print_cmd = ' '.join(cmd)
		if out_fn:
			if verbose:
				print 'ooSubprocess: ' + print_cmd + ' > ' + out_fn
			ofile = open(out_fn, 'w')
			return subprocess.call(cmd, stdout = ofile)
		else:
			if verbose:
				print 'ooSubprocess: ' + print_cmd
			return subprocess.call(cmd)


	def chain(self, prog, args = [], out_fn = None, stop = False, in_proc = None, verbose = True):
		if in_proc is None and self.chain_cmds != []:
			print 'Error: the pipeline was not stopped before creating a new one!'
			print 'In cach: %s'%(' | '.join(self.chain_cmds))
			exit(1)
		if out_fn and stop == False:
			print 'Error: out_fn (output_file_name) is only specified when stop = True!'
			exit(1)

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
			return subprocess.call(cmd, stdout = ofile, stdin = in_pipe)
		else:
			return subprocess.Popen(cmd, stdout = subprocess.PIPE, stdin = in_pipe)

