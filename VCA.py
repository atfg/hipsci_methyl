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

#data = '/lustre/scratch109/sanger/ag14/hipsci/methyl/limix_out_030415_imputedgt/pheno_pgm/gt-20150128_me-20150128.maf10.hdf5'
#out_dir = '/lustre/scratch109/sanger/ag14/hipsci/methyl/limix_out_030415_imputedgt/pheno_pgm/'
#withcovs = 'yes'
#nfolds=800
#fold_j=0

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
            #Y  = SP.log10(data.Y[:,gene_i:gene_i+1])
            Y  = data.Y[:,gene_i:gene_i+1]
            #Y -= Y.mean(0)
            #Y /= Y.std(0)
            #print "getGenotypes"
            try:
                X, info = data.getGenotypes(geneID,cis_window=5E4,center=False)
            except: continue
            #if res is not None: 
            #    X, info = res    
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
            
            print "cis/trans/noise + covariates variance decomposition"
            vc = VAR.VarianceDecomposition(Y)
            vc.addFixedEffect()
            vc.addRandomEffect(data.kgender)
            vc.addRandomEffect(data.kreprog)
            vc.addRandomEffect(data.kmedia)
            vc.addRandomEffect(data.kuser)
            vc.addRandomEffect(data.ksentrix_id)
            vc.addRandomEffect(data.ksentrix_pos)
            vc.addRandomEffect(data.kplate)
            vc.addRandomEffect(data.kwell)
            vc.addRandomEffect(data.ktime)
            vc.addRandomEffect(data.kpassage)
            vc.addRandomEffect(uKcis)
            vc.addRandomEffect(uKtrans)
            vc.addRandomEffect(is_noise=True)
            conv = vc.optimize()
            weights = vc.getWeights()
            vc_ctn = vc.getVarianceComps()[0,:]
            #vc_grtctn = SP.concatenate([vc_grt,vc_ctn])
            # export
            RV = {}
            RV['vc'] = vc_ctn
            RV['conv'] = SP.array([vc.getLMLgrad()<1e-2])
            out_group = probe_group.create_group('varDecomp')
            dumpDictHdf5(RV,out_group)
                
            print 'ok'
        except:
            continue

    f.close()

# import sys
# sys.path.append('/nfs/users/nfs_a/ag14/workspace/HIPSCI_limix/')
# import h5py
# import scipy as SP
# import pylab as PL
# import scipy.stats
# import numpy as np
# import pdb
# import time
# import os
# from utils import *
# import limix
# import limix.stats.fdr as FDR
# import limix.modules.varianceDecomposition as VAR
# import limix.modules.qtl as QTL
# import data as DATA
# 
# #from CFG.settings import *
# #sys.path.insert(0,CFG['limix_path'])
# #import limix.utils.fdr as FDR
# 
# #data file
# data = '/lustre/scratch109/sanger/ag14/hipsci/methyl/limix_out_030415_imputedgt/pheno_gm/gt-20150128_me-20150128.maf10.hdf5'
# 
# def getRealPos(pos,start,end,strand):
#     rv = pos-start
#     return rv
# 
# center=False
# 
# if __name__ == "__main__":
#     print 'Number of arguments:', len(sys.argv), 'arguments.'
#     print 'Argument List:', str(sys.argv)
#     
#     data = sys.argv[3]
#     out_dir = sys.argv[4]
#     withcovs = sys.argv[5]
#  
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)
# 
#     nfolds = int(sys.argv[1])   # number of folds the genelist is split into (= number of jobs)
#     fold_j = int(sys.argv[2])   # fold considered for this job
#     out_file = os.path.join(out_dir,'%d_%d.hdf5'%(nfolds,fold_j))
#     
#     #1. import data
#     data = DATA.data(data)
#     
#     uKpop = data.getUKpop(center=center,normalize=False)
#     
#     #if withcovs == 'yes':
#     gender = data.gender
#     reprog = data.reprog
#     media = data.media
#     user = data.user
#     sentrix_id = data.sentrix_id
#     sentrix_pos = data.sentrix_pos
#     plate = data.plate
#     well = data.well
#     time = data.time
#     gender /= gender.diagonal().mean()
#     
#         #uCov = SP.concatenate([SP.ones_like(gender),gender,reprog,media,tissue],1)
#         #uCov = SP.concatenate([SP.ones_like(gender),gender,reprog,media,user,sentrix_id,sentrix_pos,plate,well,time],1)
#         # build Kdonor
#         #donor  = data.donorID
#         #donors = SP.unique(donor)
#         #Kdonor = SP.zeros_like(uKpop) 
#         #for _d in donors:
#         #    vec = (1.*(donor==_d))[:,SP.newaxis]
#         #    Kdonor += SP.dot(vec,vec.T)
#     # else:
#     #    uCov = None
#     
#     #3. split gene list in parts
#     n_genes = data.Y.shape[1]
#     genes = SP.arange(n_genes)
#     Icv = SP.floor(nfolds*SP.arange(n_genes)/n_genes)
#     I = Icv==fold_j
#     genes = genes[I]
#         
#     #4. loop across gene, perform experiment
#     f = h5py.File(out_file,'w')
#         
#     for gene_i in genes:
#         probeID = data.probeID[gene_i]
#         geneID  = data.geneID[gene_i]
#         print ""
#         print "%s: %s" % (gene_i,probeID)
#         probe_group = f.create_group(probeID)
# 
#         try:
#             # save geneID
#             probe_group.create_dataset('geneID',data=SP.array([data.geneID[gene_i]]))
#             #4.1 load data for gene geneID
#             Y  = data.Y[:,gene_i:gene_i+1]
#             try:
#                 X, info = data.getGenotypes(geneID,cis_window=5E4,center=False)
#             except: continue         
#             #4.1b get geneID info
#             strand = data.gene_strand[gene_i]
#             start = int(data.gene_start[gene_i])
#             end = int(data.gene_end[gene_i])
#                 
#             # export probe cis info
#             out_group = probe_group.create_group('info')
#             dumpDictHdf5(info,out_group)
#             probe_group.create_dataset('strand',data=SP.array([strand]))
#             probe_group.create_dataset('start',data=SP.array([start]))
#             probe_group.create_dataset('end',data=SP.array([end]))
# 
#             # one line per donor
#             Xu = np.array(X, dtype='float')
#             Yu = np.array(Y, dtype='float')
#             if center:
#                 Xu -= Xu.mean(0); Xu /= Xu.std(0)   
# 
#             uKcis    = SP.dot(Xu,Xu.T)
#             uKtrans  = uKpop-uKcis
#             uKcis   /= uKcis.diagonal().mean()
#             uKtrans /= uKtrans.diagonal().mean()
#                          
#             #4.3 perform experiment and store results in out_gene
#             out_gene = {}
#             
#             if 1:
#                 print "cis/trans/noise + covariates variance decomposition"
#                 vc = VAR.VarianceDecomposition(Y)
#                 vc.addFixedEffect()
#                 vc.addRandomEffect(gender)
#                 vc.addRandomEffect(reprog)
#                 vc.addRandomEffect(media)
#                 vc.addRandomEffect(user)
#                 vc.addRandomEffect(sentrix_id)
#                 vc.addRandomEffect(sentrix_pos)
#                 vc.addRandomEffect(plate)
#                 vc.addRandomEffect(well)
#                 vc.addRandomEffect(time)
#                 vc.addRandomEffect(uKcis)
#                 vc.addRandomEffect(uKtrans)
#                 vc.addRandomEffect(is_noise=True)
#                 conv = vc.optimize()
#                 weights = vc.getWeights()
#                 vc_grt = SP.array([SP.var(weights[1,0]*gender),SP.var(weights[2,0]*reprog),SP.var(weights[3,0]*media),
#                                    SP.var(weights[4,0]*user),SP.var(weights[5,0]*sentrix_id),SP.var(weights[6,0]*sentrix_pos),
#                                    SP.var(weights[7,0]*plate),SP.var(weights[8,0]*well),SP.var(weights[9,0]*time)])
#                 vc_ctn = vc.getVarianceComps()[0,:]
#                 vc_grtctn = SP.concatenate([vc_grt,vc_ctn])
#                 # export
#                 RV = {}
#                 RV['vc'] = vc_grtctn
#                 RV['conv'] = SP.array([vc.getLMLgrad()<1e-2])
#                 out_group = probe_group.create_group('varDecomp')
#                 dumpDictHdf5(RV,out_group)
#                         
#             if 0:
#                 print "gender/reprog/donor/noise variance decomposition"
#                 _Y = Y[data.Ifibr,:]
#                 _Y -= _Y.mean(0); 
#                 _Y /= _Y.std(0)
#                 _gender = gender[data.Ifibr,:]
#                 _reprog = reprog[data.Ifibr,:]
#                 _Kdonor = Kdonor[data.Ifibr,:][:,data.Ifibr]
#                 vc = VAR.VarianceDecomposition(_Y)
#                 vc.addFixedEffect()
#                 vc.addFixedEffect(_gender)
#                 vc.addFixedEffect(_reprog)
#                 vc.addRandomEffect(_Kdonor)
#                 vc.addRandomEffect(is_noise=True)
#                 conv = vc.optimize()
#                 weights = vc.getWeights()
#                 vc_gr = SP.array([SP.var(weights[1,0]*gender),SP.var(weights[2,0]*reprog)])
#                 vc_dn = vc.getVarianceComps()[0,:]
#                 vc_drgn = SP.concatenate([vc_gr,vc_dn])[[2,1,0,3]]
#                 # export
#                 RV = {}
#                 RV['vc_drgn'] = vc_drgn
#                 RV['conv'] = SP.array([vc.getLMLgrad()<1e-2])
#                 out_group = probe_group.create_group('drgn_varDecomp')
#                 dumpDictHdf5(RV,out_group)
# 
#             if 0:
#                 print "gender/repr/cis/trans/noise vd"
#                 uGender = gender[data.uI,:]
#                 uReprog = reprog[data.uI,:]
#                 vc = VAR.VarianceDecomposition(Yu)
#                 vc.addFixedEffect()
#                 vc.addFixedEffect(uGender)
#                 vc.addFixedEffect(uReprog)
#                 vc.addRandomEffect(uKcis)
#                 vc.addRandomEffect(uKtrans)
#                 vc.addRandomEffect(is_noise=True)
#                 conv = vc.optimize()
#                 weights = vc.getWeights()
#                 vc_gr = SP.array([SP.var(weights[1,0]*uGender),SP.var(weights[2,0]*uReprog)])
#                 vc_ctn = vc.getVarianceComps()[0,:]
#                 vc_ctrgn = SP.concatenate([vc_gr,vc_ctn])[[2,3,1,0,4]]
#                 # export
#                 RV = {}
#                 RV['vc_ctrgn'] = vc_ctrgn
#                 RV['conv'] = SP.array([vc.getLMLgrad()<1e-2])
#                 out_group = probe_group.create_group('ctrgn_varDecomp')
#                 dumpDictHdf5(RV,out_group)
#             
#             print 'ok'
#                 
#         except:
#             continue
# 
#     f.close()
