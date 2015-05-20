#!/bin/python
"""Marks/removes duplicates using NuGEN's Unique Molecular Barcode
technology.

It removes reads as duplicates if they fulfill the following criteria: a) start
in the same location b) same orientation c) have the same unique molecular barcode. The read 
with the highest mapping quality is not considered a duplicate.  For paired end reads, 
the criteria for a duplicate is: a) start in the same location b) same template length 
c) have the same unique molecular barcode.

The runtime for marking and removal of duplicates depends on the format of the
inputs. Here are the three cases for running this script.
	Case 1 (fastest runtime): Pre-built SAM/BAM
		input_sam -- sorted SAM/BAM with fixed length UMI-sequence appended to the read name
		umi_length -- length of UMI sequence (defaults to 6)
	Case 2 (almost as fast as case 1): Pairwise matched SAM/BAM and Index fastq
	  SAM/BAM records do not have a UMI, but are in a one-to-one ordered
	  match to sequence records in index fastq. The UMI is a substring of an index
	  sequence record.
		input_sam -- SAM/BAM
		index_fastq -- Fastq file containing sequences 
		start -- starting position (1-based) of UMI in each sequence record in the index fastq (defaults to 6)
		umi_length -- length of UMI sequence (defaults to 6)
	Case 3 (slowest): Unmatched SAM/BAM and Index fastq
	  SAM/BAM records do not have a UMI, and are not a particular order
	  compared to index fastq.  Index fastq read names must be unique, and each
	  SAM/BAM record must have corresponding read name in the fastq file.
		input_sam -- SAM/BAM
		index_fastq -- Fastq file containing sequences 
		start -- starting position (1-based) of UMI in each sequence record in the index fastq (defaults to 6)
		umi_length -- length of UMI sequence (defaults to 6)

Note: Script autodetects between Case 2 and Case 3.
"""

import os
import shutil
import sys
import subprocess as sp
import shlex

import re
from operator import itemgetter
from itertools import imap, groupby, chain,islice
import tempfile
import contextlib
import logging
	
MAX_MEMORY = 4* (2**30) # Default to 4 GB of memory
DEFAULT_TMP = ""
logger = logging.getLogger('trimmark')
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(logging.Formatter('%(asctime)s [%(levelname) 9s] - %(message)s'))
logger.addHandler(ch)
	
logger.setLevel(20)
IUPAC = 'ATCGRYMKSWHBVDN'

##### Grabbed from trimmark.flow.fq
ALLOWED_FASTQ = ['.fq','.fastq.gz']

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

########## grab from trimmark.utils.io
@contextlib.contextmanager
def make_named_pipe(*args, **kwargs):
	temp_dir = tempfile.mkdtemp(*args, **kwargs)
	try:
		path = os.path.join(temp_dir, 'named_pipe')
		os.mkfifo(path)
		yield path
	finally:
		shutil.rmtree(temp_dir)

class SubprocessChainToStdout(object):
	""" Creates sequential Popens, where output for the ith Popen is redirected to
	 the i+1th Popen. 

	  First cmd reads from filepath, and the stdout must be read from last process before exit

	WARNING: By default, std error are redirected to devnull. 
	"""
	def __init__(self, cmds, suppress_stderr=False):
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

	def __enter__(self):
		for idx, cmd in enumerate(self._cmds):
			# First Popen has no stdin otherwise connect previous Popen
			stdin = None if idx==0 else self._popen[-1].stdout

			self._popen.append(sp.Popen(shlex.split(cmd),stdin=stdin, stdout=sp.PIPE, stderr=self._ferr))
	
		return self

	def _flush(self):		
		""" Flush out stdout of previous processes """
		for p in self._popen[:-1]:
			p.stdout.close()
	def get_last_proc(self):
		return self._popen[-1]

	def __exit__(self, type, value, tb):
		self._flush()		
		# Wait for final subprocess to finish
		self._popen[-1].wait()

		if not self._ferr is None:
			self._ferr.close()	

	def __str__(self):
		return ' | '.join(self._cmds)


class SubprocessChain(SubprocessChainToStdout):
	""" Creates sequential Popens, where output for the ith Popen is redirected to
	 the i+1th Popen. 

	  First cmd reads from filepath, and the final output must be to a file. If the last cmd does not have a '{out}', then 'tee {out}' is added as the final command.

	WARNING: By default, std error are redirected to devnull. 
	"""
	def __init__(self, cmds, final_out, suppress_stderr=False):
		""" Init
		
		Arguments:
		  cmds -- iterable of basic parsing commands run in shell (str)
		  final_out -- file path to write stdout of last process (str)
		"""
		SubprocessChainToStdout.__init__(self, cmds, suppress_stderr=suppress_stderr)
		self._out_path = final_out
		self._fout = open(os.devnull, 'w')
		

	def __enter__(self):
		if not '{out}' in self._cmds[-1]:
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
		SubprocessChainToStdout.__exit__(self,type,value,tb)
		self._fout.close()

class MarkRmDups(object):
	""" Makes two bam files, one with PCR duplicates marked, one with duplicates removed.

	The input must be a sorted SAM file where each read has the last UMI_LENGTH letters representing a unique molecular identifier. A PCR duplicate is a read, which maps to the same contig, positions, orientation, and has the same molecular identifier as another read with greater MAPQ score.  
	"""

	def __init__(self, out_prefix, Writer=None):
		self._out_mark_suffix = '.sorted.markdup.bam'
		self._out_rm_suffix = '.sorted.dedup.bam'
		self._out_dir = os.path.normpath(out_prefix)
		logger.debug("Outputting to %s", self._out_dir)
		self.unaligned_count = None
		self.umi_dup_count = None
		self.dup_count = None
			
		self._umi = None
		if Writer is None:
			self._writer = DefaultWriter(self)
		else:
			self._writer = Writer(self)

		self._tomark = None # process to write stdin to
		self._torm = None # process to write stdin to

	def set_umi_length(self, length=6, validator=False):
		if validator == False:
			self._umi = UMISeq(length)
		else:
			self._umi = ValidUMISeq(length)

	def iter_unique_position(self, sorted_sam_fh):
		""" Generator that yields on groups of sam rows, where contig and pos are identical """
		contig = None
		pos = None
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
			elif pos<n_pos or contig!=n_contig or contig is None:
				# handles initial call, if contigs and positions are different 
				if len(grp)>0:
					yield grp
				grp = [row,]
				contig = n_contig
				pos = n_pos
			else:
				#  Only reaches here if pos>n_pos and contig==l_contig
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
		return self._out_dir + self._out_mark_suffix
	def get_rmdup_path(self):
		if self._out_rm_suffix is None:
			return self.rmdup_path
		return self._out_dir + self._out_rm_suffix

	def mark_from_sorted_sam_with_umi_in_header(self, sorted_sam_fh):
		""" Marks UMI duplicates in one bam, and Removes UMI duplicates in another bam """

		self.umi_dup_count = 0
		self.dup_count = 0
		samtobam_fmt = 'samtools view -bS -o {out_bam} -'
	
		
		if not os.path.isdir(os.path.dirname(self._out_dir)):
			if self._out_dir == '.':
				self._out_dir = os.path.join(self._out_dir, '')
			elif not os.sep in self._out_dir:
				self._out_dir = os.path.join('.',self._out_dir)
			else:
				raise IOError('directory path not found %s'%self._out_dir)

		self._fnull = open(os.devnull, 'wb')
		#self._fnull = None
		# open pipes to send data to 
		tomark_cmd =samtobam_fmt.format(out_bam=self.get_markdup_path()) 
		logger.debug(tomark_cmd)
		self._tomark = sp.Popen(shlex.split(tomark_cmd), stdin=sp.PIPE, stderr=self._fnull)

		torm_cmd =samtobam_fmt.format(out_bam=self.get_rmdup_path()) 
		logger.debug(torm_cmd)
		self._torm = sp.Popen(shlex.split(torm_cmd), stdin=sp.PIPE, stderr=self._fnull)
		
		#self._tomark_raw = open(self.get_markdup_path().replace('.bam','.sam'), 'wb')
		#self._torm_raw = open(self.get_rmdup_path().replace('.bam','.sam'), 'wb')
		
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

		self._tomark.stdin.close()
		self._torm.stdin.close()
		self._tomark.wait()
		self._torm.wait()
		try:
			self._fnull.close()
		except:
			pass


class MarkRmDupsPairedEnd(MarkRmDups):
	""" Makes two bam files, one with PCR duplicates marked, one with duplicates removed from Paired End data

	The input must be a sorted SAM file where each read has the last UMI_LENGTH letters representing a unique molecular identifier. A PCR duplicate is a read, which maps to the same contig, positions, orientation, and has the same molecular identifier as another read with greater MAPQ score.  
	"""

	def __init__(self, out_prefix, Writer=None):
		MarkRmDups.__init__(self, out_prefix, Writer=Writer)

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
		  (SAM Flag 0x40 set) and not reverse complemented (SAM Flag 0x10 unset)
		"""
		contig = None
		pos = None
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
				n_tlen = int(row[8])
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
			elif n_flag & 0x1 == 0 or n_flag & 0x2 == 0 or n_flag & 0x10 > 0 or n_tlen<=0:
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
			elif n_flag & 0x10 == 0 and n_flag & 0x20 == 0x20 and n_tlen>0 and  pos==n_pos and contig==n_contig:
				# consider only Forward/reverse mapping reads, with this being the first read.
				grp.append(row)
			elif n_flag & 0x10 == 0 and n_flag & 0x20 == 0x20 and n_tlen>0 and  pos<n_pos or contig!=n_contig or contig is None:
				# handles initial call, if contigs and positions are different 
				for pgrp in self._iter_grp_and_post(grp,interleaved_with_grp):
						yield pgrp
				#interleaved_with_grp = []
				grp = [row,]
				contig = n_contig
				pos = n_pos
			else:
				#  Only reaches here if pos>n_pos and contig==l_contig
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

		self.umi_dup_count = 0
		self.dup_count = 0
		samtobam_fmt = 'samtools view -bS -o {out_bam} -'
		
		if not os.path.isdir(os.path.dirname(self._out_dir)):
			if self._out_dir == '.':
				self._out_dir = os.path.join(self._out_dir, '')
			elif not os.sep in self._out_dir:
				self._out_dir = os.path.join('.',self._out_dir)
			else:
				raise IOError('directory path not found %s'%self._out_dir)

		self._fnull = open(os.devnull, 'wb')
		tomark_cmd =samtobam_fmt.format(out_bam=self.get_markdup_path()) 
		logger.debug(tomark_cmd)
		self._tomark = sp.Popen(shlex.split(tomark_cmd), stdin=sp.PIPE, stderr=self._fnull)

		torm_cmd =samtobam_fmt.format(out_bam=self.get_rmdup_path()) 
		logger.debug(torm_cmd)
		self._torm = sp.Popen(shlex.split(torm_cmd), stdin=sp.PIPE, stderr=self._fnull)
	
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
						if is_best:
							# executed exactly once.
							self.set_flag_and_write(sam_row, False)
							is_reverse_read_a_dup[sam_row[0]] = False
							is_best = False
						else:
							self.set_flag_and_write(sam_row, True)
							is_reverse_read_a_dup[sam_row[0]] = True
							self.umi_dup_count += 1 
					#logger.debug('Key %s, Unique Dup Count, %s: Len(rows)=%s', key, umi_dup_count, len(sam_rows))

		if len(is_reverse_read_a_dup)>0:
			logger.error('There are %d reverse end reads, that were not caught and set properly.', len(is_reverse_read_a_dup))
		for k,v in is_reverse_read_a_dup.iteritems():
			logger.error('QNAME: %s IsDup: %d', k,v)

		self._tomark.stdin.close()
		self._torm.stdin.close()
		self._tomark.wait()
		self._torm.wait()
		try:
			self._fnull.close()
		except:
			pass

class DefaultWriter(object):
	def __init__(self, parent):
		self.parent = parent
	def write(self, sam_row, to_rm_flag=True):
		""" Writes toRM and toMark open processes """
		sam_str = '\t'.join(sam_row)
		if to_rm_flag:
			self.parent._torm.stdin.write(sam_str)
		self.parent._tomark.stdin.write(sam_str)

class UMIStripWriter(DefaultWriter):
	def __init__(self, parent):
		DefaultWriter.__init__(self,parent)
	def write(self, sam_row, to_rm_flag=True):	
		""" Writes toRM and toMark open processes, strips UMI from the header """
		if not sam_row[0].startswith('@'):
			# strip umi_length + 1 character from the end of a header
			sam_row[0] = sam_row[0][-(self.parent._umi.length+1):]
		sam_str = '\t'.join(sam_row)
		if to_rm_flag:
			self.parent._torm.stdin.write(sam_str)
		self.parent._tomark.stdin.write(sam_str)

class UMISeq(object):
	""" Extracts UMI sequence from SAM row."""
	def __init__(self, length):
		self.length = length
	def get_umi_seq(self, sam_row):	
		try:
			return sam_row[0][-self.length:]  
		except TypeError, IndexError:
			raise Exception('Need to set umi length, to get umi from read name suffix')

class ValidUMISeq(UMISeq):
	def __init__(self, length):
		UMISeq.__init__(self,length)
		self.umi_re = None
	def get_umi_seq(self, sam_row):
		""" Validates the UMI seq belongs to the IUPAC alphabet """
		try:
			umi_seq = sam_row[0][-self.length:]  
		except TypeError, IndexError:
			raise Exception('Need to set umi length, to get umi from read name suffix')

		# Validates the UMI code matches a string of fixed umi length
		try:
			if not self.umi_re.match(umi_seq):
				raise AssertionError('UMI sequence, {0} does not match IUPAC format. Is UMI in SAM header?'.format(umi_seq))
		except AttributeError:
			# If no umi pattern is compiled, then sets adn returns the new umi_seq
			logger.debug("Setting UMI pattern for validation")
			self.umi_re = re.compile('^[%s]{%d}$'%(IUPAC,self.length))
			return self.get_umi_seq(sam_row)
				
		return umi_seq

class PrepDeDup(object):
	""" Autodetects modular flow of file processing required for marking and removing deduplicated reads
	in sorted bam files.

	Depends on subprocesses ran in a BASH shell. Requires GNU-coreutils (cut,diff,grep,gzip,sed) and samtools. 
	Creates BASH scripts, and runs them for sanity checks and processing files on the stream.	
	"""

	def __init__(self, sam_file, index_fq_file=None, out_prefix=''):
		""" Init 

		Arguments:
		  sam_file -- absolute path to sam/bam sorted/unsorted file (str)
		  index_fq -- absolute path to fastq or gzipped fastq (str)
		  out_prefix -- absolute prefix of path to write sorted BAM files to (str)
		"""
		
		self._sam = sam_file
		self._index = None if index_fq_file is None else FQFilePath(index_fq_file)
		self._out_prefix = out_prefix

		# Sets how to grab sam regardless of bam or sam files	
		if self.is_bam():
			self._nohead_sam_cmd = "samtools view {sam}".format(sam=self._sam)
		else:
			self._nohead_sam_cmd = "grep -v '^@' {sam}".format(sam=self._sam)

		# Sets how to parse fastq files, must pass through
		# series of piped processes
		if self._index:
			sed_pair_reformat = r"$!N;s/^@\?\([^ \t]\+\).*\n/\1\t/"
			if self._index.is_compressed():
				self._fq_cmds = ("gzip -c -d {fq}".format(fq=self._index.fq), "sed -e '{re}'".format(re=sed_pair_reformat),  "sed -n 1~2p")
			else:
				self._fq_cmds = ("sed -e '{re}' {fq}".format(fq=self._index.fq,re=sed_pair_reformat), "sed -n 1~2p")

		# how to obtain UMI from fq sequence
		self._umi_from_fq_fmt = r"""awk '{{print $1 "\t" substr($2,length($2)-{s:d}+1,{l:d});}}'"""

		self.DupMain = MarkRmDups

	def is_bam(self):
		return self._sam[-4:]=='.bam'
	def is_sam_and_fq_synced(self):
		""" Check if sam and index are row wise synced with the same read name """
		if self._index is None:
			return False
		
		with make_named_pipe(prefix=DEFAULT_TMP) as fqcut_path, make_named_pipe(prefix=DEFAULT_TMP) as samcut_path:
			
			# Process reads from named pipes.
			diff_cmd = "diff -q {sam} {fq}".format(sam=samcut_path, fq=fqcut_path)
			diff_p = sp.Popen(shlex.split(diff_cmd), stdout=sp.PIPE, stderr=sp.PIPE)
		
			# Chain of commands to output to named pipe
			cut = 'cut -f 1'
			fqcut_chain = chain(self._fq_cmds, [cut,])	
			samcut_chain = [self._nohead_sam_cmd, cut]
			with SubprocessChain(samcut_chain, samcut_path, suppress_stderr=True) as samcut, SubprocessChain(fqcut_chain, fqcut_path, suppress_stderr=True) as fqcut:		
				logger.debug("Shell FQ parsing one-liner: %s", fqcut)
				logger.debug("Shell SAM parsing one-liner: %s", samcut)

				diff_out, diff_err = diff_p.communicate()
			logger.debug("Diff out %s, err %s", diff_out, diff_err)
			if len(diff_out.strip())==0:
				return True
			else:
				return False

	def process_unsynced_sam(self, umi_start, umi_length):
		""" Processing sam that is not synced. Also performs Sanity Checks on RNAME key"""
		if umi_length is None:
			raise Exception('Requires a UMI length')
		if umi_start is None:
			raise Exception('Requires a UMI start')

		logger.info("Processing Unsynced SAM/BAM with UMI found in Index Fastq")
		w = self.DupMain(self._out_prefix, Writer=UMIStripWriter)
		w.set_umi_length(umi_length)

		with tempfile.NamedTemporaryFile(prefix=DEFAULT_TMP, suffix='_sorted.fq') as fq_path, tempfile.NamedTemporaryFile(prefix=DEFAULT_TMP, suffix='_sorted.sam') as sam_path, tempfile.NamedTemporaryFile(prefix=DEFAULT_TMP, suffix='_umi_sorted.bam') as sorted_umi_sam:

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
				logger.error("SAM QNAMEs did not match QNAMEs in Index\n %s", sam_uniq_out[:1024])
				raise Exception('SAM and Index do not match check failed')

			# A little not obvious. Requires thinking of a pipe, backwards
			# Define the reader, before you start writing.

			# Declare all the pipes we will use first
			with make_named_pipe(prefix=DEFAULT_TMP) as umi_sam_path, make_named_pipe(prefix=DEFAULT_TMP) as sam_head_path, make_named_pipe(prefix=DEFAULT_TMP) as umi_join_path:
	
				# SAMwUMIwHead --> BAMwUMI --> Sorted BAMwUMI --> Sorted SAMwUMIwHead
				samtobam = 'samtools view -bhS {sam}'.format(sam=umi_sam_path)
				#sort_bam_cmd = 'samtools sort -o -m {mem} - -'.format(mem=MAX_MEMORY)
				sort_bam_cmd = 'samtools sort -m {mem} - {{out}}'.format(mem=MAX_MEMORY)
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
							logger.debug('UMI sam %s', umi_nohead_sam)
			#raw_input('continue? [enter] \n')
			self._sam = sorted_umi_sam.name
			w = self.process_sorted_sam_with_umi_in_rname(umi_length)
			#sorted_umi_sam.seek(0)
			#w.mark_from_sorted_sam_with_umi_in_header(sorted_umi_sam)
		return w
	def process_pairwise_sam(self, umi_start, umi_length):
		""" Processins pairwise sam"""
		if umi_length is None:
			raise Exception('Requires a UMI length')
		if umi_start is None:
			raise Exception('Requires a UMI start')

		logger.info("Processing Synced SAM/BAM with UMI found in Index Fastq")
		w = self.DupMain(self._out_prefix, Writer=UMIStripWriter)
		w.set_umi_length(umi_length)

		# A little not obvious. Requires thinking of a pipe, backwards
		# Define the reader, before you start writing.

		# Declare all the pipes we will use first
		with make_named_pipe(prefix=DEFAULT_TMP) as fq_path, make_named_pipe(prefix=DEFAULT_TMP) as sam_path, make_named_pipe(prefix=DEFAULT_TMP) as umi_sam_path, make_named_pipe(prefix=DEFAULT_TMP) as sam_head_path, make_named_pipe(prefix=DEFAULT_TMP) as umi_join_path, tempfile.NamedTemporaryFile(prefix=DEFAULT_TMP, suffix='_umi_sorted.bam') as sorted_umi_sam:
	
			# SAMwUMIwHead --> BAMwUMI --> Sorted BAMwUMI --> Sorted SAMwUMIwHead
			samtobam = 'samtools view -bhS {sam}'.format(sam=umi_sam_path)
			#sort_bam_cmd = 'samtools sort -o -m {mem} - -'.format(mem=MAX_MEMORY)
			#bamtosam = 'samtools view -h -'
			sort_bam_cmd = 'samtools sort -m {mem} - {{out}}'.format(mem=MAX_MEMORY)
			#with SubprocessChain([samtobam, sort_bam_cmd, bamtosam], sorted_umi_sam.name) as sort_sam:
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

					# Add SAM file with UMI in header
					join_cmd = "join -t '	' --nocheck-order -j 1 {fq} {sam}".format(fq=fq_path, sam=sam_path)
					add_umi_cmd = r"sed -e 's/^\([^\t]\+\)\t\([^\t]\+\)/\1:\2/'"
					with SubprocessChain([join_cmd, add_umi_cmd], umi_join_path) as umi_nohead_sam:
						logger.debug('UMI sam %s', umi_nohead_sam)

						# Parse FQ and SAM ad add UMI
						# Sort and Validate matching names"
						fq_chain = chain(self._fq_cmds, [self._umi_from_fq_fmt.format(s=umi_start, l=umi_length)])
						with SubprocessChain(fq_chain, fq_path) as fq_p, SubprocessChain([self._nohead_sam_cmd,], sam_path) as sam_p:
							logger.debug("Shell FQ: %s", fq_p)
							logger.debug("Shell SAM: %s", sam_p)


			self._sam = sorted_umi_sam.name
			w = self.process_sorted_sam_with_umi_in_rname(umi_length)
			#sorted_umi_sam.seek(0)
			#w.mark_from_sorted_sam_with_umi_in_header(sorted_umi_sam)
		return w
	def process_sorted_sam_with_umi_in_rname(self, umi_length):
		""" Processes a SAM/BAM file with the a fixed length UMI sequence appended to read name """
		if umi_length is None:
			raise Exception('Requires a UMI length')
		w = self.DupMain(self._out_prefix)
		w.set_umi_length(umi_length, validator=True)
		logger.info("Processing Sorted SAM/BAM with UMI in suffix of read name (assumes sorted sam)")
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
				w.mark_from_sorted_sam_with_umi_in_header(f)
	
			return w

	def _log_output(self, dup_obj):
		logger.info('      Unaligned count: {0:010d}'.format(dup_obj.unaligned_count))
		logger.info('Positional dups count: {0:010d}'.format(dup_obj.dup_count))
		logger.info('        N6 dups count: {0:010d}'.format(dup_obj.umi_dup_count))
		logger.info('Deduplication success.')
	
	def main(self,umi_start=None,umi_length=None):
		""" Checks/validates the inputs and identifies a pipeline to use for marking and removal of duplicates """
		logger.info('Deduplicating single end reads...')
		if self._index is None:
			w = self.process_sorted_sam_with_umi_in_rname(umi_length)
		elif self.is_sam_and_fq_synced():
			w = self.process_pairwise_sam(umi_start, umi_length)
		else:
			logger.info('SAM and FQ are not synced')
			w = self.process_unsynced_sam(umi_start, umi_length)	
		self._log_output(w)
		return w

class PrepDeDupPairedEnd(PrepDeDup):
	""" Treats Input SAM as containing Paired End data, and marks/removes potential PCR duplicate
	  templates
	"""
	def __init__(self, sam_file, index_fq_file=None, out_prefix=''):
		PrepDeDup.__init__(self, sam_file, index_fq_file, out_prefix)

		self.DupMain = MarkRmDupsPairedEnd
	def main(self,umi_start=None,umi_length=None):
		logger.info('Deduplicating paired end reads...')
		if self._index is None:
			w = self.process_sorted_sam_with_umi_in_rname(umi_length)
		elif self.is_sam_and_fq_synced():
			w = self.process_pairwise_sam(umi_start, umi_length)
		else:
			logger.info('SAM and FQ are not synced')
			w = self.process_unsynced_sam(umi_start, umi_length)	
		self._log_output(w)
		return w

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
		parser.error("The file %s does not exist!" %arg)
	else:
		return str(arg)

if __name__ == '__main__':
	import argparse, sys
	
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-2','--paired-end', dest='pe', action='count', default=0, help="SAM contains paired end reads, remove potential PCR template duplicates.")
	parser.add_argument('--index-fq', dest='index_fq', type=lambda x:file_check(parser, x), default=None, help='index file paired with input sam')
	parser.add_argument('-o','--out', dest='out_prefix', default='./prefix', help='prefix of output file paths for sorted BAMs (default will create ./prefix.sorted.markdup.bam, ./prefix.sorted.dedups.bam)')
	parser.add_argument('-s','--start', dest='start', type=int, default=6, help="where to start considering the UMI (default = last 6 bases are considered)")
	parser.add_argument('-l','--length', dest='length', type=int, default=6, help="length of UMI sequence (default = 6)")
	#parser.add_argument('-l', help="log file to write statistics to (optional)")
	parser.add_argument('-v','--version', action='version', version='%(prog)s 0.1')
	parser.add_argument('sam', type=lambda x:file_check(parser, x), help='input sorted/unsorted SAM/BAM')

	args = parser.parse_args()

	if args.pe > 0:
		w = PrepDeDupPairedEnd(args.sam, index_fq_file=args.index_fq, out_prefix=args.out_prefix)
	else:
		w = PrepDeDup(args.sam, index_fq_file=args.index_fq, out_prefix=args.out_prefix)


	w.main(umi_start=args.start, umi_length=args.length)
	#from_sorted_bam(sys.argv[1], sys.argv[2], int(sys.argv[3]))
	#from_sorted_sam_stdin(sys.argv[1], int(sys.argv[2]))
	
