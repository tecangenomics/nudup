#!/usr/bin/env python2.7
"""Marks/removes PCR introduced duplicate molecules based on the molecular tagging
technology used in NuGEN products.

For SINGLE END reads, duplicates are marked if they fulfill the following
criteria: a) start at the same genomic coordinate b) have the same strand
orientation c) have the same molecular tag sequence. The read with the
highest mapping quality is kept as the non-duplicate read.

For PAIRED END reads, duplicates are marked if they fulfill the following
criteria: a) start at the same genomic coordinate b) have the same template
length c) have the same molecular tag sequence. The read pair with the highest
mapping quality is kept as the non-duplicate read.

Here are the two cases for running this tool:

- CASE 1 (Standard): User supplies two input files,
 1) SAM/BAM file that a) ends with .sam or .bam extension b) contains unique
    alignments only
 2) FASTQ file (ie: the Index FASTQ) that contains the molecular tag sequence
    for each read name in the corresponding SAM/BAM file as either a) the read
    sequence or b) in the FASTQ entry header name. If the index FASTQ read
    length is 6, 8, 12, 14, or 16nt long as expected for NuGEN products, the
    molecular tag sequence to be extracted from the read according to -s and -l
    parameters, otherwise the molecular tag will be extracted from the header
    of the FASTQ entry.

- CASE 2 (Runtime Optimized): User supplies only one input file,
 1) SAM/BAM file that a) ends with .sam or .bam extension b) contains unique
    alignments only c) is sorted d) has a fixed length sequence containing the
    molecular tag appended to each read name.

Author: {author}
Contact: {company}, {email}
"""

import contextlib
from itertools import chain, groupby, islice
import logging
import os
import re
import shlex
import shutil
import subprocess as sp
import sys
import tempfile
import time

__author__ = 'Anand Patel'
__company__ = "NuGEN Technologies Inc."
__email__ = "techserv@nugen.com"
__version__ = '2.3'
__copyright__ = """Copyright (C) 2015 NuGEN Technologies Inc.
   
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>."""


IUPAC = 'ATCGRYMKSWHBVDN'
ALLOWED_FASTQ = ['.fq','.fastq.gz']

MAX_MEMORY = 4*(2**30) # in bytes

logger = logging.getLogger('ote')

# Define Logging parameters"
logger.setLevel(logging.DEBUG)

# Log to Console
ch = logging.StreamHandler(sys.stdout)
#ch.setLevel(logging.INFO)
ch.setLevel(logging.DEBUG)

# Specifies format of log
ch.setFormatter(logging.Formatter('%(asctime)s [%(levelname) 9s] - %(message)s'))
logger.addHandler(ch)

@contextlib.contextmanager
def make_named_pipe(*args, **kwargs):
	temp_dir = tempfile.mkdtemp(*args, **kwargs)
	try:
		path = os.path.join(temp_dir, 'named_pipe')
		os.mkfifo(path)
		yield path
	finally:
		shutil.rmtree(temp_dir)

class SubprocessChain(object):
	""" Creates sequential Popens, where output for the ith Popen is redirected to
	 the i+1th Popen. 

	  First cmd reads from filepath, and the final output must be to a file. If the last cmd does not have a '{out}', then 'tee {out}' is added as the final command. Each command must finish reading/writing from pipe otherwise errors are thrown.

	WARNING: By default, std error are redirected to devnull. 
	"""
	def __init__(self, cmds, final_out=None, suppress_stderr=False, stdout=None):
		""" Init
		
		Arguments:
		  cmds -- iterable of basic parsing commands run in shell (str)
		  final_out -- file path to write stdout of last process (str)
		"""
		self._cmds = list(cmds)
		if suppress_stderr:
			self._ferr = open(os.devnull, 'w')
		else:
			self._ferr = None
		self._popen = []

		if final_out is None and stdout is None:
			raise Exception('Both final_out and stdout cannot be None')
		self._out_path = final_out
		self._stdout = stdout
		if self._stdout is None:
			self._fout = open(os.devnull, 'w')
		else:
			self._fout = self._stdout

		self._POLL_TIME = 2 # Number of seconds to wait between polls of process finishes
		self.errors = [] # error return code for commands in the process		

	def __str__(self):
		return ' | '.join(self._cmds)

	def _flush(self):		
		""" Flush out stdout of previous processes """
		for p in self._popen[:-1]:
			p.stdout.close()

	def __enter__(self):
		if self._stdout is None and not '{out}' in self._cmds[-1]:
			self._cmds.append('tee {out}')
		self._cmds[-1] = self._cmds[-1].format(out=self._out_path)
		for idx, cmd in enumerate(self._cmds):
			# PIPE all output except quench the last popen
			stdout = sp.PIPE if idx<len(self._cmds)-1 else self._fout
			# First Popen has no stdin otherwise connect previous Popen
			stdin = None if idx==0 else self._popen[-1].stdout

			self._popen.append(sp.Popen(shlex.split(cmd),stdin=stdin, stdout=stdout, stderr=self._ferr))
	
		return self
		
	def __exit__(self, type, value, tb):
		self._flush()		
		# Loop for processes to finish
		returncode = dict((p,p.poll()) for p in self._popen)

		success = None
		while success is None:
			returncode = dict((p,p.poll()) for p in self._popen)
			for ret in returncode.itervalues():
				if ret>0:
					success = False
			returned = list((v for v in returncode.itervalues() if not v is None))
			if len(returned)==len(self._popen) and sum(returned)==0:
				success = True

			if success is None:
				time.sleep(self._POLL_TIME)

		# All processes finished
		popen_idx = dict(zip(self._popen,xrange(len(self._popen))))
		if not success:
			for p, ret in returncode.iteritems():
				if p.poll() is None:
					p.terminate()
				elif ret>0:
					self.errors.append((popen_idx[p],returncode[p]))
		
		if not self._ferr is None:
			self._ferr.close()	
		if self._stdout is None:
			self._fout.close()

	def get_last_process(self):
		return self._popen[-1]

class FQFilePath(object):
	""" Manages type of FastQ file 

	Attributes:
	  fq -- name of fq file (str)
	  fq_type -- kind of fq file a value matching ALLOWED_FASTQ
	  fq_dir -- directory fq file is found in
	  fq_base -- basename of fq file

	"""
	def __init__(self, fq_file):
		"""Init

		Checks if fq file conforms to an ALLOWED_FASTQ file
		If it conforms, strip the ALLOWED_FASTQ extension from the basename

		Arguments:
			fq_file -- fastq file (str)
		"""
		self.fq = '' if fq_file is None else fq_file
		
		fq_type = filter(self.fq.endswith, ALLOWED_FASTQ)
		self.fq_type = fq_type[0] if len(fq_type)>0 else None

		self.fq_dir, self.fq_base = os.path.split(self.fq)
		self.fq_dir = os.path.normpath(self.fq_dir)

		if not self.fq_type is None:
			self.fq_base = self.fq_base[:-1*len(self.fq_type)]	

	def is_allowed(self):
		""" Checks if file conforms to an ALLOWED_FASTQ type see __init__ """
		return False if self.fq_type is None else True
			
	def is_compressed(self):
		""" Checks if file has a gzip extension """
		return (not self.fq_type is None) and self.fq_type.endswith('.gz')

	def uncompress(self, new_fpath):
		if not self.is_compressed():
			os.link(self.fq, new_fpath)		
		else:
			gzip_cmd = 'gzip -c -d {in_fq}'.format(in_fq=self.fq)
			print gzip_cmd
			with open(new_fpath, 'wb') as f:
				gzip_p = sp.check_call(shlex.split(gzip_cmd), stdout=f)
		new_fq = FQFilePath(new_fpath)
		assert not new_fq.is_compressed()
		return new_fq	

	def head(self, n=64):
		""" Get first n lines of a fastq file """
		if self.is_compressed():
			gzip_cmd = 'gzip -c -d {in_fq}'.format(in_fq=self.fq)
			head_cmd = 'head -n {0:d}'.format(n)
			#print gzip_cmd + ' | ' + head_cmd
			with open(os.devnull, 'w') as fnull:
				gzip_p = sp.Popen(shlex.split(gzip_cmd), stdout=sp.PIPE, stderr=fnull)
				head_p = sp.Popen(shlex.split(head_cmd), stdin=gzip_p.stdout, stdout=sp.PIPE)
				gzip_p.stdout.close()
				lines = head_p.communicate()[0].split('\n')		
		else:
			with open(self.fq, 'rb') as f:
				lines = islice(f, 0,n)
		return lines

class UMIStripFunc(object):
	def __init__(self, parent):
		self._umi_trim_length = -(parent._umi.length+1)
	def __call__(self, sam_row):
		if not sam_row[0].startswith('@'):
			sam_row[0] = sam_row[0][:self._umi_trim_length]
		return sam_row

class AbstractWriter(object):
	""" Abstract writer for sam tokens (reads) to sam files 
	
		sam_row_func - performs a manipulation on sam rows
	"""
	def __init__(self, sam_row_func=None):
		self._func = sam_row_func
		self._fnull = open(os.devnull, 'wb')

	def _get_sam_str(self, sam_row):
		if self._func is None:
			return '\t'.join(sam_row)
		else:
			return '\t'.join(self._func(sam_row))

	def write(self, sam_row, to_rm_flag=True):
		raise NotImplementedError()

	def close(self):
		self._fnull.close()
	

class RmdupWriter(AbstractWriter):
	""" By Default writes markdup and rmdup files using samtools view converter """
	def __init__(self, rmdup_path, sam_row_func=None):
		AbstractWriter.__init__(self,sam_row_func)
		
		self.rmdup_path = rmdup_path
		#self._fnull = None
		# open pipes to send data to 
		samtobam_fmt = 'samtools view -bS -o {out_bam} -'
		
		torm_cmd =samtobam_fmt.format(out_bam=rmdup_path) 
		logger.debug(torm_cmd)
		self._torm = sp.Popen(shlex.split(torm_cmd), stdin=sp.PIPE, stderr=self._fnull)
	
	def _write_str(self, sam_str, to_rm_flag=True):
		""" Writes toRM open processes """
		if to_rm_flag:
			self._torm.stdin.write(sam_str)

	def write(self, sam_row, to_rm_flag=True):
		""" Writes toRM open processes """
		if to_rm_flag:
			sam_str = self._get_sam_str(sam_row)
			self._torm.stdin.write(self._get_sam_str(sam_row))
	def close(self):
		""" Flushes all input to subprocess and waits """
		self._torm.stdin.close()
		self._torm.wait()
		AbstractWriter.close(self)


class MarkdupWriter(AbstractWriter):
	""" Writes markdup file using samtools view converter """
	def __init__(self, markdup_path, sam_row_func=None):
		AbstractWriter.__init__(self,sam_row_func)
		self.markdup_path = markdup_path
		
		# open pipes to send data to 
		samtobam_fmt = 'samtools view -bS -o {out_bam} -'
		
		tomark_cmd =samtobam_fmt.format(out_bam=markdup_path) 
		logger.debug(tomark_cmd)
		self._tomark = sp.Popen(shlex.split(tomark_cmd), stdin=sp.PIPE, stderr=self._fnull)
	
	def _write_str(self, sam_str, to_rm_flag=True):
		""" Writes toMark open processes """
		self._tomark.stdin.write(sam_str)

	def write(self, sam_row, to_rm_flag=True):
		""" Writes toMark open processes """
		self._tomark.stdin.write(self._get_sam_str(sam_row))

	def close(self):
		""" Flushes all input to subprocess and waits """
		self._tomark.stdin.close()
		self._tomark.wait()
		AbstractWriter.close(self)

class RmdupMarkdupWriter(AbstractWriter):
	""" Delegates write functions to MarkdupWriter and RmdupWriter """
	def __init__(self, rmdup_path, markdup_path, sam_row_func=None):
		AbstractWriter.__init__(self, sam_row_func)
		self._torm_writer = RmdupWriter(rmdup_path)
		self._tomark_writer = MarkdupWriter(markdup_path)

	def write(self, sam_row, to_rm_flag=True):
		""" Writes to both rmdup and markdup files """
		sam_str = self._get_sam_str(sam_row)
		self._torm_writer._write_str(sam_str, to_rm_flag)
		self._tomark_writer._write_str(sam_str, to_rm_flag)

	def close(self):
		self._torm_writer.close()
		self._tomark_writer.close()
		AbstractWriter.close(self)

class UMISeq(object):
	""" Extracts UMI sequence from SAM row."""
	def __init__(self, length):
		self.length = length
	def get_umi_seq(self, sam_row):	
		try:
			return sam_row[0][-self.length:]  
		except TypeError, IndexError:
			raise Exception('Need to set molecular tag sequence length, to get umi from read name suffix')

class ValidUMISeq(UMISeq):
	def __init__(self, length):
		UMISeq.__init__(self,length)
		self.umi_re = None
	def get_umi_seq(self, sam_row):
		""" Validates the UMI seq belongs to the IUPAC alphabet """
		try:
			umi_seq = sam_row[0][-self.length:]  
		except TypeError, IndexError:
			raise Exception('Need to set molecular tag sequence length, to get molecular tag sequence from read name suffix')

		# Validates the UMI code matches a string of fixed umi length
		try:
			if not self.umi_re.match(umi_seq):
				raise AssertionError('molecular tag sequence, {0} does not match IUPAC format. Is molecular tag sequence in SAM read name?'.format(umi_seq))
		except AttributeError:
			# If no umi pattern is compiled, then sets adn returns the new umi_seq
			logger.debug("Setting molecular tag sequence pattern for validation")
			self.umi_re = re.compile('^[%s]{%d}$'%(IUPAC,self.length))
			return self.get_umi_seq(sam_row)
				
		return umi_seq

class MarkRmDups(object):
	""" Makes two bam files, one with PCR duplicates marked, one with duplicates removed.

	The input must be a sorted SAM file where each read has the last UMI_LENGTH letters representing a unique molecular identifier. A PCR duplicate is a read, which maps to the same contig, positions, orientation, and has the same molecular identifier as another read with greater MAPQ score.  
	"""

	def __init__(self, out_prefix):
		self._out_mark_suffix = '.sorted.markdup.bam'
		self._out_rm_suffix = '.sorted.dedup.bam'
		self._out_prefix = os.path.normpath(out_prefix)
		logger.debug("Outputting to %s", self._out_prefix)
		self.aligned_count = None
		self.unaligned_count = None
		self.umi_dup_count = None
		self.dup_count = None
			
		self._umi = None
		self._writer = None

		if os.sep in self._out_prefix and not os.path.isdir(os.path.dirname(self._out_prefix)):
			raise IOError('directory path not found %s'%self._out_prefix)

	def set_writer(self, writer):
		""" Sets the writer after allowing modifications to output paths """
		if not self._writer is None:
			raise RuntimeError("Setting writer multiple times")
			return
		self._writer = writer

	def set_umi_length(self, length=6, validator=False):
		if validator:
			self._umi = ValidUMISeq(length)
		else:
			self._umi = UMISeq(length)

	def iter_unique_position(self, sorted_sam_fh):
		""" Generator that yields on groups of sam rows, where contig and pos are identical """
		contig = None
		pos = None
		self.aligned_count = 0
		self.unaligned_count = 0
		grp = []
		for line in sorted_sam_fh:
			#logger.debug(line)
			if line.startswith('@'):
				# just passes header lines back to 
				yield [[line,],]
				continue

			# Splits the sam row into a list
			row = line.split('\t')
			try:
				n_contig = row[2]
				n_pos = int(row[3])	
			except KeyError:
				raise KeyError('Not a valid SAM entry.')
			#logger.debug('Parsing... %s, %d, %s', n_contig, n_pos, row)
			
			if n_contig == '*':
				# unaligned reads returned, line by line
				if len(grp)>0:
					yield grp
					grp = []
				contig = None
				pos = None
				yield [row,]
				self.unaligned_count += 1
			elif pos==n_pos and contig==n_contig:
				grp.append(row)
				self.aligned_count += 1
			elif pos<n_pos or contig!=n_contig or contig is None:
				# handles initial call, if contigs and positions are different 
				if len(grp)>0:
					yield grp
				grp = [row,]
				contig = n_contig
				pos = n_pos
				self.aligned_count += 1
			else:
				#  Only reaches here if pos>n_pos and contig==l_contig
				self.get_unique_molecule_id(row) # Catch UMI Seq validation errors
				raise Exception('Sam file is not sorted.')

		if len(grp)>0:
			yield grp

	def _check_is_forward(self, sam_row):
		# Get the 0x10 (5th) bit of the SAM flag and check if it is 0
		return int(sam_row[1]) & 0x10 == 0
	
	def get_unique_molecule_id(self, sam_row):
		""" Gets a key from the sam_row that identifies unique molecules 
			
			Unique molecules have the same mapping orientation, and 
			same umi_length string found at the end of a read name.

			Returns: (read is forward, umi sequence) """
		return (int(sam_row[1]) & 0x10, self._umi.get_umi_seq(sam_row))

	def set_flag_and_write(self, sam_row, is_dup):
		flag = int(sam_row[1])
		if is_dup:
			# mark the duplicate
			if flag & 0x400 == 0:
				flag += 0x400
			sam_row[1] = "%d"%flag
			self._writer.write(sam_row, to_rm_flag=False)
		else:
			if flag & 0x400 > 0:
				flag -= 0x400
			sam_row[1] = "%d"%flag
			self._writer.write(sam_row)

	def set_markdup_path(self, mark_path):
		self._out_mark_suffix = None
		self.markdup_path = mark_path
	def set_rmdup_path(self, rmdup_path):
		self._out_rm_suffix = None
		self.rmdup_path = rmdup_path

	def get_markdup_path(self):
		if self._out_mark_suffix is None:
			return self.markdup_path
		return self._out_prefix + self._out_mark_suffix
	def get_rmdup_path(self):
		if self._out_rm_suffix is None:
			return self.rmdup_path
		return self._out_prefix + self._out_rm_suffix

	def mark_from_sorted_sam_with_umi_in_header(self, sorted_sam_fh):
		""" Marks UMI duplicates in one bam, and Removes UMI duplicates in another bam """
		#self.set_writer()
		self.umi_dup_count = 0
		self.dup_count = 0
		
		for sam_rows in self.iter_unique_position(sorted_sam_fh):
			#logger.debug(sam_rows)
			if len(sam_rows)==1:
				# output singletons (includes both @ lines and unique_positions)
				self._writer.write(sam_rows[0])
			elif len(sam_rows)>1:
				#logger.debug('Count of sam rows in group, %s', len(sam_rows))
				self.dup_count += len(sam_rows)-1
		
				for key, umi_sam_rows in groupby(sorted(sam_rows, key=self.get_unique_molecule_id), self.get_unique_molecule_id):
					# rows matching a unique molecule
					is_best = True
					#umi_dup_count = 0 
					for sam_row in sorted(umi_sam_rows, key=lambda x:int(x[4]), reverse=True):
						# sorts rows by maximum MAPq score
						flag = int(sam_row[1])
						#if is_best and sam_row[4]!='255':
						if is_best:
							# executed exactly once.
							self.set_flag_and_write(sam_row, False)
							is_best = False
						else:
							# mark the duplicate
							self.set_flag_and_write(sam_row, True)
							self.umi_dup_count += 1 
							#umi_dup_count += 1
					#logger.debug('Key %s, Unique Dup Count, %s: Len(rows)=%s', key, umi_dup_count, len(sam_rows))

		self._writer.close()

class MarkRmDupsPairedEnd(MarkRmDups):
	""" Makes two bam files, one with PCR duplicates marked, one with duplicates removed from Paired End data

	The input must be a sorted SAM file where each read has the last UMI_LENGTH letters representing a unique molecular identifier. A PCR duplicate is a read, which maps to the same contig, positions, orientation, and has the same molecular identifier as another read with greater MAPQ score.  
	"""

	def __init__(self, out_prefix):
		MarkRmDups.__init__(self, out_prefix)
		
	def _iter_grp_and_post(self, grp, pgrp):
		if len(grp) > 0:
			yield grp
			for row in pgrp:
				yield [row,]
		grp[:] = []
		pgrp[:] = []

	def iter_unique_position(self, sorted_sam_fh):
		""" Generator that yields on groups of forward sam rows, where contig and pos are identical 
		  Forward reads that are properly aligned (SAM Flag 0x2 unset) with the same starting position
		  will be returned as a group. Additional criteria, Forward reads must be first in template 
		  (SAM Flag 0x40 set) and not reverse complemented (SAM Flag 0x10 unset). Aligned_count counts only forward
		  reads with concordant f/r paired mappings. 
		"""
		contig = None
		pos = None
		self.aligned_count = 0
		self.unaligned_count = 0
		grp = []
		interleaved_with_grp = []

		for line in sorted_sam_fh:
			#logger.debug(line)
			if line.startswith('@'):
				# just passes header lines back to 
				yield [[line,],]
				continue

			# Splits the sam row into a list
			row = line.split('\t')
			try:
				n_contig = row[2]
				n_pos = int(row[3])
				n_flag = int(row[1])
				n_contig_next = row[6]
				n_pos_next = int(row[7])
			except KeyError:
				raise KeyError('Not a valid SAM entry: {0}.'.format(line))
			#logger.debug('Parsing... %s, %d, %s', n_contig, n_pos, row)
			
			if n_contig == '*':
				# unaligned reads returned, line by line
				for pgrp in self._iter_grp_and_post(grp,interleaved_with_grp):
						yield pgrp
				#grp = []
				#interleaved_with_grp = []
				contig = None
				pos = None
				yield [row,]
				self.unaligned_count += 1
			elif n_flag & 0x1 == 0 or n_flag & 0x2 == 0 or n_flag & 0x10 > 0 or n_contig_next!='=' or n_pos_next<n_pos:
				# Do not group the following cases:
				#  - template does not have multiple segments in sequence
				#  - Unproperly aligned according to aligner
				#  - Is a reverse complemented (if paired with a groupable end, then will be catched by parent function
				#  - template length is 0
				if len(grp)>0:
					# Finish the group, and then 
					interleaved_with_grp.append(row)
				else:
					yield [row,]
			elif n_flag & 0x10 == 0 and n_flag & 0x20 == 0x20 and  pos==n_pos and contig==n_contig:
				# consider only Forward/reverse mapping reads, with this being the first forward read.
				self.aligned_count += 1
				grp.append(row)
			elif n_flag & 0x10 == 0 and n_flag & 0x20 == 0x20 and pos<n_pos or contig!=n_contig or contig is None:
				# handles initial call, if contigs and positions are different 
				for pgrp in self._iter_grp_and_post(grp,interleaved_with_grp):
						yield pgrp
				#interleaved_with_grp = []
				self.aligned_count += 1
				grp = [row,]
				contig = n_contig
				pos = n_pos
			else:
				#  Only reaches here if pos>n_pos and contig==l_contig
				self.get_unique_molecule_id(row) # Catch UMI Seq validation errors
				raise Exception('Sam file is not sorted. or Unexpected read with FLAG %d. Expects 0x1, 0x2, 0x10, and 0x20 to be set.', n_flag)

		for pgrp in self._iter_grp_and_post(grp,interleaved_with_grp):
				yield pgrp

	def get_unique_molecule_id(self, sam_row):
		""" Gets a key from the sam_row that identifies unique molecules 
			
			Unique molecules have the same mapping orientation,
			same umi_length string found at the end of a read name, and
			same template length.

			Returns: (umi sequence, tlen) """
		return (self._umi.get_umi_seq(sam_row), int(sam_row[8]))

	def mark_from_sorted_sam_with_umi_in_header(self, sorted_sam_fh):
		""" Marks UMI duplicates in one bam, and Removes UMI duplicates in another bam 
			Since paired end, the last segment in template read must be marked or removed depending on
			the grouping of the first properly mapped read.
		"""
		#self.set_writer()
		self.umi_dup_count = 0
		self.dup_count = 0
	
		is_reverse_read_a_dup = {} 
	
		for sam_rows in self.iter_unique_position(sorted_sam_fh):
			#logger.debug(sam_rows)
			if len(sam_rows)==1:
				# output singletons (includes both @ lines and unique_positions)
				qname = sam_rows[0][0]
				if qname in is_reverse_read_a_dup:
					self.set_flag_and_write(sam_rows[0], is_reverse_read_a_dup.pop(qname))
					
				else:	
					self._writer.write(sam_rows[0])
			elif len(sam_rows)>1:
				#logger.debug('Count of sam rows in group, %s', len(sam_rows))
				self.dup_count += len(sam_rows)-1

				for key, umi_sam_rows in groupby(sorted(sam_rows, key=self.get_unique_molecule_id), self.get_unique_molecule_id):
					# rows matching a unique molecule
					is_best = True
					for sam_row in sorted(umi_sam_rows, key=lambda x:int(x[4]), reverse=True):
						# sorts rows by maximum MAPq score
						flag = int(sam_row[1])
						#if is_best and sam_row[4]!='255':
						if is_best:
							# executed exactly once.
							is_reverse_read_a_dup[sam_row[0]] = False
							self.set_flag_and_write(sam_row, False)
							is_best = False
						else:
							is_reverse_read_a_dup[sam_row[0]] = True
							self.set_flag_and_write(sam_row, True)
							self.umi_dup_count += 1 
					#logger.debug('Key %s, Unique Dup Count, %s: Len(rows)=%s', key, umi_dup_count, len(sam_rows))

		if len(is_reverse_read_a_dup)>0:
			logger.error('There are %d reverse end reads, that were not caught and set properly.', len(is_reverse_read_a_dup))
		for k,v in is_reverse_read_a_dup.iteritems():
			logger.error('READ NAME: %s IsDup: %d', k,v)
		self._writer.close()
	
class IndexFinder(object):
	def __init__(self, fq_file):
		self._fq = FQFilePath(fq_file)
		self.index_seq_length = None # Length of index sequences
		self.index_seq_count = None

	def _check_valid_sequences_by_count(self, seq_lines, seq_words, seq_chars, max_seq, fq_check):
		""" Checks if all index sequences have the same length """
		if seq_lines!=seq_words:
			logger.info('Not valid %s FASTQ File, ill-formatted Index sequences, possibly spaces', fq_check)
			return False
		elif (seq_chars/(max_seq+1))!=seq_lines:
			logger.info('Not valid %s FASTQ File, differing lengths of Index sequences', fq_check)
			return False

		self.index_seq_length =  max_seq
		return True
	def _get_check_index_cmds(self):
		cmds = ['cat {fq}'.format(fq=self._fq.fq)]
		if self._fq.is_compressed():
			cmds.append('gzip -c -d')
		cmds.append('sed -n 2~4p')
		return cmds

	def _get_index_cmds(self):
		cmds = ['cat {fq}'.format(fq=self._fq.fq)]
		if self._fq.is_compressed():
			cmds.append("gzip -c -d")

		# Grabs the first word of the read title and tabs it with the FASTQ sequence (index sequence)
		sed_pair_reformat = r"$!N;s/^@\?\([^ \t]\+\).*\n/\1\t/"
		cmds.append("sed -e '{re}'".format(re=sed_pair_reformat))
		cmds.append("sed -n 1~2p")
		return cmds

	def _get_check_read_title_cmds(self):
		cmds = ['cat {fq}'.format(fq=self._fq.fq)]
		if self._fq.is_compressed():
			cmds.append('gzip -c -d')
		cmds.append('sed -n 1~4p')
		cmds.append('grep -o "[ACGTN]\{6,\}"')
		return cmds

	def _get_read_title_cmds(self):
		cmds = ['cat {fq}'.format(fq=self._fq.fq)]
		if self._fq.is_compressed():
			cmds.append("gzip -c -d")

		cmds.append('sed -n 1~4p')
		sed_reformat = r"s/^@\([^ \t]\+\).*\([ACGTN]\{6,\}\).*/\1\t\2/"
		cmds.append("sed -e '{re}'".format(re=sed_reformat))
		return cmds

	def is_index_fastq(self):
		""" Index FASTQ has sequence lengths 6,8,14,12,16 """ 
		cmds = self._get_check_index_cmds()
		# Always ordered: newline, word, character,  byte,  maximum  line
		cmds.append('wc -lwcmL')
		with SubprocessChain(cmds, stdout=sp.PIPE) as pchain:
			logger.debug('Checking FASTQ is Index fastq, %s', pchain)
			stdout, stderr = pchain.get_last_process().communicate()
			seq_lines, seq_words, seq_chars, seq_byte, max_seq = map(int,stdout.strip().split())
		self.index_seq_count = seq_lines
		if not max_seq in [6,8,12,14,16]:
			return False

		return self._check_valid_sequences_by_count(seq_lines, seq_words, seq_chars, max_seq, 'Index')
	def is_index_in_read_title(self):
		""" Check if [ACGTN]6+ in read title """
		# Requires knowing the number of Records in FASTQ file to be valid
		assert not self.index_seq_count is None 
		cmds = self._get_check_read_title_cmds()
		# Always ordered: newline, word, character,  byte,  maximum  line
		cmds.append('wc -lwcmL')
		with SubprocessChain(cmds, stdout=sp.PIPE) as pchain:
			logger.debug('Checking FASTQ is Read fastq, %s', pchain)
			stdout, stderr = pchain.get_last_process().communicate()
			seq_lines, seq_words, seq_chars, seq_byte, max_seq= map(int,stdout.strip().split())
		return self._check_valid_sequences_by_count(seq_lines, self.index_seq_count, seq_chars, max_seq, 'Read')

	def get_cmds_for_parsing_fastq(self):
		if self.is_index_fastq():
			logger.info('Using molecular tag sequence from Index FASTQ read')
			return self._get_index_cmds()
		elif self.is_index_in_read_title():
			logger.info('Using molecular tag sequence from FASTQ header name')
			return self._get_read_title_cmds()
		else:
			logger.error('No Valid molecular tag sequence information found in FASTQ header name, please provide a valid Index FASTQ file')
			sys.exit(1)

class PrepDeDup(object):
	""" Autodetects modular flow of file processing required for marking and removing deduplicated reads
	in sorted bam files.

	Depends on subprocesses ran in a BASH shell. Requires GNU-coreutils (cut,diff,grep,gzip,sed) and samtools. 
	Creates BASH scripts, and runs them for sanity checks and processing files on the stream.	
	"""

	def __init__(self, sam_file, tmp_prefix, fq_file=None, out_prefix='', old_samtools=False, rmdup_only=False):
		""" Init 

		Arguments:
		  sam_file -- absolute path to sam/bam sorted/unsorted file (str)
		  tmp_prefix -- prefix for reading/writing temp files and named pipes (str)
		  fq_file -- absolute path to fastq or gzipped fastq (str)
		  out_prefix -- absolute prefix of path to write sorted BAM files to (str)
		  old_samtools -- flag for using old samtools sort style (boolean)
		  rmdup_only -- only reports a sam file without duplicated reads (boolean)
		"""
		
		self._sam = sam_file
		self._tmp_prefix = tmp_prefix		

		self._out_prefix = out_prefix

		self._old_samtools = old_samtools
		self._rmdup_only = rmdup_only

		# Sets how to grab sam regardless of bam or sam files	
		if self.is_bam():
			self._nohead_sam_cmd = "samtools view {sam}".format(sam=self._sam)
		else:
			self._nohead_sam_cmd = "grep -v '^@' {sam}".format(sam=self._sam)

		# series of piped processes
		
		self._index = None if fq_file is None else FQFilePath(fq_file)

		# how to obtain UMI from fq sequence
		self._umi_from_fq_fmt = r"""awk '{{print $1 "\t" substr($2,length($2)-{s:d}+1,{l:d});}}'"""
		self.DupMain = MarkRmDups
		self.type_str = 'single'


	def is_bam(self):
		return self._sam[-4:]=='.bam'

	def _set_dup_writer(self, w, sam_row_func=None):
		""" Factory method for Writer """
		if self._rmdup_only:
			w.set_writer(RmdupWriter(w.get_rmdup_path(), sam_row_func))
		else:
			w.set_writer(RmdupMarkdupWriter(w.get_rmdup_path(), w.get_markdup_path(), sam_row_func))
		
	def process_unsynced_sam(self, umi_start, umi_length):
		""" Processing sam that is not synced. Also performs Sanity Checks on RNAME key"""
		if umi_length is None:
			raise Exception('Requires a molecular tag sequence length')
		if umi_start is None:
			raise Exception('Requires a molecular tag sequence start')
		if self._index is None:
			raise Exception('Requires an index FASTQ')

		# Sets how to parse fastq files, must pass through
		self._index_parser = IndexFinder(self._index.fq)
		self._fq_cmds = self._index_parser.get_cmds_for_parsing_fastq()
		
		if self._index_parser.index_seq_length < umi_start:
			logger.error("Specified start {0:d} is before the length of Index sequence {1:d}".format(umi_start, self._index_parser.index_seq_length))
			

		logger.info("Appending molecular tag sequence to SAM/BAM read name")
		w = self.DupMain(self._out_prefix)
		w.set_umi_length(umi_length)
		self._set_dup_writer(w, UMIStripFunc(w))

		with tempfile.NamedTemporaryFile(prefix=self._tmp_prefix, suffix='_sorted.fq') as fq_path, tempfile.NamedTemporaryFile(prefix=self._tmp_prefix, suffix='_sorted.sam') as sam_path, tempfile.NamedTemporaryFile(prefix=self._tmp_prefix, suffix='_umi_sorted.bam') as sorted_umi_sam:

			sort_by_name = 'sort -k 1b,1 -S {mem}b'.format(mem=MAX_MEMORY/2)
			# Parse FQ and SAM, add UMI, sort by name
			fq_chain = chain(self._fq_cmds, [self._umi_from_fq_fmt.format(s=umi_start, l=umi_length), sort_by_name])
			with SubprocessChain(fq_chain, fq_path.name) as fq_p, SubprocessChain([self._nohead_sam_cmd,sort_by_name], sam_path.name, suppress_stderr=True) as sam_p:
					logger.debug("Shell FQ: %s", fq_p)
					logger.debug("Shell SAM: %s", sam_p)
		
	
			# Checks if the Index has unique rnames
			uniq_cmd = 'sort -k 1b,1 -c -u {fq}'.format(fq=fq_path.name)
			logger.debug(uniq_cmd)
			index_uniq_out = sp.check_output(shlex.split(uniq_cmd))
			if len(index_uniq_out.strip())>0:
				logger.error("Unique index check failed: %s", index_uniq_out[:1024])
				raise Exception('Index is not unique failed')

			# Checks if each RNAME in SAM is in the FQ
			uniq_to_sam_cmd = "join -t '	' -j 1 -v 2 {fq} {sam}".format(fq=fq_path.name, sam=sam_path.name)
			logger.debug(uniq_to_sam_cmd)
			sam_uniq_out = sp.check_output(shlex.split(uniq_to_sam_cmd))	
			if len(sam_uniq_out.strip())>0:
				logger.error("SAM read names did not match read names in Index\n %s", sam_uniq_out[:1024])
				raise Exception('SAM and Index do not match check failed')

			# A little not obvious. Requires thinking of a pipe, backwards
			# Define the reader, before you start writing.

			# Declare all the pipes we will use first
			with make_named_pipe(prefix=self._tmp_prefix) as umi_sam_path, make_named_pipe(prefix=self._tmp_prefix) as sam_head_path, make_named_pipe(prefix=self._tmp_prefix) as umi_join_path:
	
				# SAMwUMIwHead --> BAMwUMI --> Sorted BAMwUMI --> Sorted SAMwUMIwHead
				samtobam = 'samtools view -bhS {sam}'.format(sam=umi_sam_path)
				if self._old_samtools:
					sort_bam_cmd = 'samtools sort -m {mem} - {{out}}'.format(mem=MAX_MEMORY)
				else:			
					sort_bam_cmd = 'samtools sort -m {mem} -T {prefix} -o {{out}}.bam'.format(prefix=self._tmp_prefix+'sort', mem=MAX_MEMORY)

				#bamtosam = 'samtools view -h -'
				#with SubprocessChain([samtobam, sort_bam_cmd, bamtosam], sorted_umi_sam.name, suppress_stderr=True) as sort_sam:
				with SubprocessChain([samtobam, sort_bam_cmd], sorted_umi_sam.name[:-4], suppress_stderr=True) as sort_sam:
					logger.debug('Sorting %s',sort_sam)

					# Redirection for header SAMwUMIwHead <-- SAM/BAMwHead <-- SAMwUMI
					cat_cmd = "cat {head} {umi_join}".format(head=sam_head_path, umi_join=umi_join_path)
					with SubprocessChain([cat_cmd,], umi_sam_path) as umi_sam:
						logger.debug('Cat %s', umi_sam)
						# Adds SAM header for sorting with bam
						sam_opt = '' if self.is_bam() else 'S'
						add_header_cmd = 'samtools view -H{opt} -o {head} {sam}'.format(opt=sam_opt, sam=self._sam, head=sam_head_path)
						logger.debug('Add header %s',add_header_cmd)
						with open(os.devnull, 'w') as fnull:
							p = sp.check_call(shlex.split(add_header_cmd), stderr=fnull)
						if p>0:
							logger.error('Could not find header in SAM file')
							raise Exception('Could not find head of SAM file. ')

						# Add SAM file with UMI in header
						join_cmd = "join -t '	' -j 1 {fq} {sam}".format(fq=fq_path.name, sam=sam_path.name)
						add_umi_cmd = r"sed -e 's/^\([^\t]\+\)\t\([^\t]\+\)/\1:\2/'"	
						with SubprocessChain([join_cmd, add_umi_cmd], umi_join_path) as umi_nohead_sam:
							logger.debug('Add molecular tag sequence to sam: %s', umi_nohead_sam)
			#raw_input('continue? [enter] \n')
			self._sam = sorted_umi_sam.name
			self.process_sorted_sam_with_umi_in_rname(umi_length, dup_obj=w)
			#sorted_umi_sam.seek(0)
			#w.mark_from_sorted_sam_with_umi_in_header(sorted_umi_sam)
		return w

	def process_sorted_sam_with_umi_in_rname(self, umi_length, dup_obj=None):
		""" Processes a SAM/BAM file with the a fixed length UMI sequence appended to read name """
		if umi_length is None:
			raise Exception('Requires a molecular tag sequence length')
		if dup_obj is None:
			w = self.DupMain(self._out_prefix)
			w.set_umi_length(umi_length, validator=True)
			self._set_dup_writer(w)
		else:
			w = dup_obj
		logger.info("Processing sorted SAM/BAM with molecular tag sequence in read name (assumes sorted)")
		if self.is_bam():
			tosam_cmd = 'samtools view -h {bam}'.format(bam=self._sam)  
			logger.debug(tosam_cmd)
			tosam_p = sp.Popen(shlex.split(tosam_cmd), stdout=sp.PIPE)
			try:
				w.mark_from_sorted_sam_with_umi_in_header(tosam_p.stdout)
			finally:
				tosam_p.stdout.close()
				tosam_p.wait()
			return w
		else:
			with open(self._sam, 'rb') as f:
				try:
					w.mark_from_sorted_sam_with_umi_in_header(f)
				except AssertionError as e:
					logger.error(e.message)
					sys.exit(1)
			return w

	def _log_output(self, dup_obj):
		try:
			p_rate = dup_obj.dup_count/float(dup_obj.aligned_count)
		except ZeroDivisionError:
			p_rate = 0
		try:
			d_rate = dup_obj.umi_dup_count/float(dup_obj.aligned_count)
		except ZeroDivisionError:
			d_rate = 0
		logger.info( '           Aligned count: {0: 13d}'.format(dup_obj.aligned_count))
		logger.info( '         Unaligned count: {0: 13d}'.format(dup_obj.unaligned_count))
		logger.debug('   Positional dups count: {0: 13d} ({1:0.4f} rate)'.format(dup_obj.dup_count,p_rate))
		logger.info( 'Molecular tag dups count: {0: 13d} ({1:0.4f} rate)'.format(dup_obj.umi_dup_count, d_rate))
		logger.info('Deduplication success.')
		with open(self._out_prefix+"_dup_log.txt", 'wb') as f:
			f.write('\t'.join(['aligned_count','unaligned_count','position_dup_count','frac_position_dup','moltag_dup_count','frac_moltag_dup'])+'\n')
			metrics = [dup_obj.aligned_count, dup_obj.unaligned_count, dup_obj.dup_count, p_rate, dup_obj.umi_dup_count, d_rate]
			f.write('\t'.join(['{0:d}','{1:d}','{2:d}','{3:0.4f}','{4:d}','{5:0.4f}']).format(*metrics)+'\n')
		
	
	def main(self,umi_start=None,umi_length=None):
		""" Checks/validates the inputs and identifies a pipeline to use for marking and removal of duplicates """
		logger.info('Deduplicating NuGEN {0} end reads...'.format(self.type_str))
		if self._index is None:
			w = self.process_sorted_sam_with_umi_in_rname(umi_length)
		else:
			w = self.process_unsynced_sam(umi_start, umi_length)	
		self._log_output(w)
		if not self._rmdup_only:
			logger.info('Created output file {0} with duplicates marked'.format(w.get_markdup_path()))	
		logger.info('Created output file {0} with duplicates removed'.format(w.get_rmdup_path()))	
		return w

class PrepDeDupPairedEnd(PrepDeDup):
	""" Treats Input SAM as containing Paired End data, and marks/removes potential PCR duplicate
	  templates
	"""
	def __init__(self, sam_file, tmp_prefix, fq_file=None, out_prefix='', old_samtools=False, rmdup_only=False):
		PrepDeDup.__init__(self, sam_file, tmp_prefix, fq_file, out_prefix, old_samtools, rmdup_only)

		self.DupMain = MarkRmDupsPairedEnd
		self.type_str = 'paired'

def from_sorted_bam(bam_fpath, out_prefix, umi_length=6):
	w = MarkRmDups(out_prefix)
	#w = MarkRmDupsUnsorted(out_prefix)
	w.set_umi_length(umi_length)

	tosam_cmd = 'samtools view -h {in_bam}'.format(in_bam=bam_fpath)  
	logger.debug(tosam_cmd)
	tosam_p = sp.Popen(shlex.split(tosam_cmd), stdout=sp.PIPE)
	try:
		w.mark_from_sorted_sam_with_umi_in_header(tosam_p.stdout)
	finally:
		tosam_p.stdout.close()
		tosam_p.wait()
	return tosam_p.returncode

def from_sorted_sam_stdin(out_prefix, umi_length=6):
	import sys
	#w = MarkRmDups(out_prefix)
	w = MarkRmDups(out_prefix)
	w.set_umi_length(umi_length, validator=True)

	w.mark_from_sorted_sam_with_umi_in_header(sys.stdin)

def file_check(parser, arg):
	if not os.path.exists(arg):
		parser.error("The file {0} does not exist!".format(arg))
	else:
		return str(arg)

def tmp_dir_check(parser, arg):
	temp_dir = os.path.abspath(arg)
	if not os.path.isdir(temp_dir):
		parser.error("The argument {0} is not a directory!".format(temp_dir))

	check_fpath = os.path.join(temp_dir, 'named_pipe')
	try:
		os.mkfifo(check_fpath)
	except OSError:
		parser.error("Must be able to make named pipes in temp directory. See `man mkfifo` for named pipe info.") 

	finally:
		if os.path.exists(check_fpath):
			os.unlink(check_fpath)

	return os.path.join(str(temp_dir), 'nudup_')
	


if __name__ == '__main__':
	import argparse, sys, logging
	default_tmp_dir = '/tmp'

	parser = argparse.ArgumentParser(description=__doc__.format(author=__author__, company=__company__, email=__email__), formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)
		
	pgroup = parser.add_argument_group("Input")
	pgroup.add_argument('sam', metavar='IN.sam|IN.bam', type=lambda x:file_check(parser, x), help='input sorted/unsorted SAM/BAM containing only unique alignments (sorted required for case 2 detailed above)')

	ogroup = parser.add_argument_group("Options")
	ogroup.add_argument('-2','--paired-end', dest='pe', action='count', default=0, help="use paired end deduping with template. SAM/BAM alignment must contain paired end reads. Degenerate read pairs (alignments for one read of pair) will be discarded.")
	ogroup.add_argument('-f', dest='fq', metavar='INDEX.fq|READ.fq', type=lambda x:file_check(parser, x), default=None, help='FASTQ file containing the molecular tag sequence for each read name in the corresponding SAM/BAM file (required only for CASE 1 detailed above)')
	ogroup.add_argument('-o','--out', dest='out_prefix', default='prefix', help='prefix of output file paths for sorted BAMs (default will create prefix.sorted.markdup.bam, prefix.sorted.dedup.bam, prefix_dup_log.txt)')
	ogroup.add_argument('-s','--start', dest='start', type=int, default=6, help="position in index read where molecular tag sequence begins. This should be a 1-based value that counts in from the 3' END of the read. (default = 6)")
	ogroup.add_argument('-l','--length', dest='length', type=int, default=6, help="length of molecular tag sequence (default = 6)")
	ogroup.add_argument('-T', dest='tmp_prefix', metavar='TEMP_DIR', type=lambda x:tmp_dir_check(parser, x), default=default_tmp_dir, help='directory for reading and writing to temporary files and named pipes (default: {0})'.format(default_tmp_dir))
	ogroup.add_argument('--old-samtools', dest='old_samtools', action='store_true', default=False, help="required for compatibility with samtools sort style in samtools versions <=0.1.19")
	ogroup.add_argument('--rmdup-only', dest='rmdup_only', action='store_true', default=False, help="required for only outputting duplicates removed file")
	ogroup.add_argument('--debug', dest='debug', action='store_true', default=False, help=argparse.SUPPRESS)
	#ogroup.add_argument('-l', help="log file to write statistics to (optional)")
	ogroup.add_argument('-v','--version', action='version', version='%(prog)s '+ __version__)
	ogroup.add_argument('-h','--help',action='help', help='show this help message and exit')

	args = parser.parse_args()

	if args.debug:
		logger.setLevel(logging.DEBUG)
	else:
		logger.setLevel(logging.INFO)	

	if args.length>args.start:
		parser.error("Invalid molecular tag subsequence: The start position {0} counting from the end is less than the length {1}".format(args.start, args.length))

	if args.pe > 0:
		w = PrepDeDupPairedEnd(args.sam, args.tmp_prefix, fq_file=args.fq, out_prefix=args.out_prefix, old_samtools=args.old_samtools, rmdup_only=args.rmdup_only)
	else:
		w = PrepDeDup(args.sam, args.tmp_prefix, fq_file=args.fq, out_prefix=args.out_prefix, old_samtools=args.old_samtools, rmdup_only=args.rmdup_only)

	if args.debug:
		w.main(umi_start=args.start, umi_length=args.length)
	else:
		try:
			w.main(umi_start=args.start, umi_length=args.length)
		except Exception as e:
			logger.error(e.message)
			sys.exit(1)
	#from_sorted_bam(sys.argv[1], sys.argv[2], int(sys.argv[3]))
	#from_sorted_sam_stdin(sys.argv[1], int(sys.argv[2]))
	
