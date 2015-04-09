import sys
sys.path.append('/nfs/users/nfs_a/ag14/workspace/HIPSCI_limix/')
import h5py
import scipy as SP
import pylab as PL
import scipy.stats
import numpy as np
import pdb
import time
import os
from utils import *
import limix
import limix.stats.fdr as FDR
import limix.modules.varianceDecomposition as VAR
import limix.modules.qtl as QTL
import data as DATA

#data = '/lustre/scratch109/sanger/ag14/hipsci/methyl/gt-20150128_me-20150128.maf10.hdf5'
#out_dir = '/lustre/scratch109/sanger/ag14/hipsci/methyl/limix_out_030415_imputedgt/'

def getRealPos(pos,start,end,strand):
	""" relative position accounting for the direction """
	#if strand=='+':			rv = pos-start
	#elif strand=='-':		rv = end-pos
	#else:					rv = SP.nanA
	rv = pos-start
	return rv

center=False

if __name__ == "__main__":
	print 'Number of arguments:', len(sys.argv), 'arguments.'
	print 'Argument List:', str(sys.argv)
	
	data = sys.argv[3]
	out_dir = sys.argv[4]

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
		
	if 'debug' in sys.argv:
		nfolds = 1
		fold_j = 0 
		out_file = os.path.join(out_dir,'./debug.hdf5')
		#pdb.set_trace()
		
	else:
		nfolds = int(sys.argv[1])   # number of folds the genelist is split into (= number of jobs)
		fold_j = int(sys.argv[2])   # fold considered for this job
		out_file = os.path.join(out_dir,'%d_%d.hdf5'%(nfolds,fold_j))

	data = DATA.data(data)
	uKpop = data.getUKpop(center=center,normalize=False)

	gender = data.gender
	#tissue = data.tissue
	reprog = data.reprog
	media = data.media
	#uCov = SP.concatenate([SP.ones_like(gender),gender,reprog,media,tissue],1)
	uCov = SP.concatenate([SP.ones_like(gender),gender,reprog,media],1)

	n_genes = data.Y.shape[1]
	genes = SP.arange(n_genes)
	Icv = SP.floor(nfolds*SP.arange(n_genes)/n_genes)
	I = Icv==fold_j
	genes = genes[I]

	f = h5py.File(out_file,'w')

	for gene_i in genes:
		probeID = data.probeID[gene_i]
		geneID  = data.geneID[gene_i]
		print ""
		print "%s: %s" % (gene_i,probeID)
		probe_group = f.create_group(probeID)
		
		try:
			# save geneID
			probe_group.create_dataset('geneID',data=SP.array([data.geneID[gene_i]]))
			#4.1 load data for gene geneID
			Y  = SP.log10(data.Y[:,gene_i:gene_i+1])
			Y -= Y.mean(0)
			Y /= Y.std(0)
			#print "getGenotypes"
			X, info = data.getGenotypes(geneID,cis_window=1E6,center=False)
			#4.1b get geneID info
			strand = data.gene_strand[gene_i]
			start = int(data.gene_start[gene_i])
			end = int(data.gene_end[gene_i])
			# export probe cis info
			out_group = probe_group.create_group('info')
			#print "dumpDict"
			dumpDictHdf5(info,out_group)
			probe_group.create_dataset('strand',data=SP.array([strand]))
			probe_group.create_dataset('start',data=SP.array([start]))
			probe_group.create_dataset('end',data=SP.array([end]))
			# one line per donor
			Xu = np.array(X, dtype='float')
			Yu = np.array(Y, dtype='float')
			#Yu -= Yu.mean(0); Yu /= Yu.std(0)
			if center:
				Xu -= Xu.mean(0); Xu /= Xu.std(0)
			uKcis    = SP.dot(Xu,Xu.T)
			uKtrans  = uKpop-uKcis
			uKcis   /= uKcis.diagonal().mean()
			uKtrans /= uKtrans.diagonal().mean()
			#4.3 perform experiment and store results in out_gene
			out_gene = {}
			#print "cis scan"
			lm=QTL.test_lmm(snps=Xu,pheno=Yu,K=uKtrans,covs=uCov,verbose=True)
			pv=lm.getPv()
			RV = {}
			RV['pv'] = pv
			RV['qv'] = FDR.qvalues(pv)[0]
			RV['lambd']   = getLambda(pv)
			RV['beta'] = lm.getBetaSNP()
			RV['posLead'] = SP.array([getRealPos(info['pos'][pv[0,:].argmin()],start,end,strand)])
			RV['aDirPosLead'] = abs(SP.array([info['pos'][pv[0,:].argmin()]-0.5*(start+end)]))
			out_group = probe_group.create_group('lmm')
			dumpDictHdf5(RV,out_group)
			#print 'ok'
		except:
			continue

	f.close()

