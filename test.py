import unittest
import tempfile
import os

from nudup import MarkRmDups
from nudup import PrepDeDup
from nudup import RmdupMarkdupWriter

TEST_DATA_DIR="/mnt/rddata/apatel/test_data"

class TestMarkRmDups(unittest.TestCase):
	def setUp(self):
		self.sam_umi_sorted = os.path.join(TEST_DATA_DIR, 'dedup','input', 'multx_umi_sorted.sam')
		self.bam_umi_sorted = os.path.join(TEST_DATA_DIR, 'dedup','input', 'multx_umi_sorted.bam')
		self.sam_sorted = os.path.join(TEST_DATA_DIR, 'dedup','input', 'multx_sorted.sam')
		self.bam_sorted = os.path.join(TEST_DATA_DIR, 'dedup','input', 'multx_sorted.bam')

		self.sam_unsorted = os.path.join(TEST_DATA_DIR, 'dedup','input', 'multx.sam')

		self.index = os.path.join(TEST_DATA_DIR, 'dedup', 'input', 'bI1.fastq')
		self.read = os.path.join(TEST_DATA_DIR, 'dedup', 'input', 'R1_space_hash.fastq')
		self.index_gz = os.path.join(TEST_DATA_DIR, 'dedup', 'input', 'bI1.fastq.gz')

		self.tmp_prefix = '/tmp/nudup_'
		self.out_prefix = os.path.join(TEST_DATA_DIR, 'dedup','output', 'multx')
		self.preview = os.path.join(TEST_DATA_DIR, 'dedup','input', 'h200_multx_umi_sorted.sam')
		
	def tearDown(self):
		pass
	def testDupSmall(self):
		w = MarkRmDups(out_prefix=self.out_prefix)
		w.set_umi_length(6)
		w.set_writer(RmdupMarkdupWriter(w.get_rmdup_path(), w.get_markdup_path()))

		with open(self.preview, 'rb') as f:
			w.mark_from_sorted_sam_with_umi_in_header(f)	
	
		self.assertTrue(os.path.isfile(self.out_prefix + w._out_mark_suffix), msg='No mark bam file produced')
		self.assertTrue(os.path.isfile(self.out_prefix + w._out_rm_suffix), msg='No rm bam file produced')
		
		self.assertEquals(w.unaligned_count,0, msg='Unaligned count is off')	
		self.assertEquals(w.dup_count,63, msg='Dup count {0} is off'.format(w.dup_count))	
		self.assertEquals(w.umi_dup_count,63, msg='UMI dup count {0} is off'.format(w.umi_dup_count))	
	def _check_for_multx(self,w):
		self.assertTrue(os.path.isfile(self.out_prefix + w._out_mark_suffix), msg='No mark bam file produced')
		self.assertTrue(os.path.isfile(self.out_prefix + w._out_rm_suffix), msg='No rm bam file produced')
		
		self.assertEquals(w.unaligned_count,8180, msg='Unaligned count is off')	
		self.assertEquals(w.dup_count,727503, msg='Dup count {0} is off'.format(w.dup_count))	
		self.assertEquals(w.umi_dup_count,520211, msg='UMI dup count {0} is off'.format(w.umi_dup_count))	

	def testDup(self):
		w = MarkRmDups(out_prefix=self.out_prefix)
		w.set_umi_length(6)
		w.set_writer(RmdupMarkdupWriter(w.get_rmdup_path(), w.get_markdup_path()))


		try:
			os.remove(self.out_prefix + w._out_mark_suffix)
			os.remove(self.out_prefix + w._out_rm_suffix)
		except:
			pass

		with open(self.sam_umi_sorted, 'rb') as f:
			w.mark_from_sorted_sam_with_umi_in_header(f)	
		self._check_for_multx(w)

	
	def testDupSyncedSamFastq(self):
		#d = PrepDeDup(self.sam_unsorted, fq_file=self.index, out_prefix=self.out_prefix)
		pass
		#w = d.main(umi_start=6, umi_length=6)
		#self._check_for_multx(w)
		raise NotImplementedError()

	def testDupUnsortedSamFastqGz(self):
		d = PrepDeDup(self.sam_sorted, self.tmp_prefix, fq_file=self.index_gz, out_prefix=self.out_prefix)
		w = d.main(umi_start=6, umi_length=6)
		self._check_for_multx(w)


	def testDupUnsortedSam(self):
		d = PrepDeDup(self.sam_sorted, self.tmp_prefix, fq_file=self.index, out_prefix=self.out_prefix)

		w = d.main(umi_start=6, umi_length=6)
		self._check_for_multx(w)
	
	def testDupUnsortedSamIndexInReadTitle(self):
		d = PrepDeDup(self.sam_sorted, self.tmp_prefix, fq_file=self.read, out_prefix=self.out_prefix)

		w = d.main(umi_start=6, umi_length=6)
		self._check_for_multx(w)

	def testOldSamtoolsDupUnsortedSamIndexInReadTitle(self):
		d = PrepDeDup(self.sam_sorted, self.tmp_prefix, fq_file=self.read, out_prefix=self.out_prefix, old_samtools=True)

		w = d.main(umi_start=6, umi_length=6)
		self._check_for_multx(w)
	
def test_all():
	
	suite = unittest.TestLoader().loadTestsFromTestCase(TestMarkRmDups)
	
	unittest.TextTestRunner(verbosity=2).run(suite)

def test():
	suite = unittest.TestSuite()
	#suite.addTest(TestMarkRmDups('testOldSamtoolsDupUnsortedSamIndexInReadTitle'))
	#suite.addTest(TestMarkRmDups('testDupUnsortedSamIndexInReadTitle'))
	suite.addTest(TestMarkRmDups('testDupSmall'))
	#suite.addTest(TestMarkRmDups('testDupSyncedBamFastq'))
	unittest.TextTestRunner(verbosity=2).run(suite)


if __name__=='__main__':
	test_all() 
	#test()
