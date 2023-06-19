'''
CWAS DAWN analysis step

This step uses the DNM counts file and the output of the categorization step, 
and then conducts DAWN analysis for the categorized DNM.
'''

import argparse
from pathlib import Path
import os, glob
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from rpy2.robjects.packages import importr
from tqdm import tqdm
from scipy.stats import norm
import random
import sys

from cwas.core.dawn.clustering import kmeans_cluster
from cwas.core.dawn.supernodeWGS import supernodeWGS_func, data_collection
from cwas.runnable import Runnable
from cwas.utils.check import check_is_file
from cwas.utils.check import check_is_dir
from cwas.utils.check import check_num_proc
from cwas.utils.log import print_arg, print_progress

class Dawn(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._eig_vector = None
        self._eig_vector_file = None
        self._corr_mat = None
        self._corr_mat_file = None
        self._permut_test = None
        self._permut_test_file = None
        self._k_val = None
        self.kmeans_r = importr('stats').kmeans

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg(
            "No. worker processes for the DAWN",
            f"{args.num_proc: ,d}",
        )
        print_arg("Input files directory: ", args.input_dir_path)
        print_arg("Output directory: ", args.output_dir_path)
        print_arg("Using DNM counts file: ", args.category_count_file)
        print_arg("Using category sets file: ", args.category_set_file)
        print_arg("Thresholds of count / correlation / size: ", ", ".join(list(map(str, [args.count_threshold, args.corr_threshold, args.size_threshold]))))

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)
        check_is_dir(args.input_dir_path)
        check_is_dir(args.output_dir_path)
        check_is_file(args.category_count_file)
        check_is_file(args.category_set_file)
    
    @property
    def num_proc(self):
        return self.args.num_proc
    
    @property
    def input_dir_path(self) -> Path:
        return self.args.input_dir_path
    
    @property
    def eig_vector_file(self):
        if self._eig_vector_file is None:
            self._eig_vector_file = self._get_input_file("*eig_vecs*.txt.gz")
        return self._eig_vector_file
    
    @property
    def eig_vector(self):
        if self._eig_vector is None:
            self._eig_vector = pd.read_table(self.eig_vector_file, compression='gzip', index_col=0)
        return self._eig_vector
    
    @property
    def corr_mat_file(self):
        if self._corr_mat_file is None:
            self._corr_mat_file = self._get_input_file("*corr_mat*.pickle")
        return self._corr_mat_file

    @property
    def corr_mat(self):
        if self._corr_mat is None:
            self._corr_mat = pd.read_pickle(self.corr_mat_file)
        return self._corr_mat

    @property
    def permut_test_file(self):
        if self._permut_test_file is None:
            self._permut_test_file = self._get_input_file("*permutation_test*.txt.gz")           
        return self._permut_test_file

    @property
    def permut_test(self):
        if self._permut_test is None:
            self._permut_test = pd.read_table(self.permut_test_file, compression='gzip')
            
        return self._permut_test
    
    @property
    def output_dir_path(self):
        return self.args.output_dir_path
    
    @property
    def k_range(self):
        return self.args.k_range
    
    @property
    def k_val(self):
        if self._k_val is None: # initialization
            if self.args.k_val is not None:
                self._k_val = self.args.k_val
            else: # self.args.k_val is None
                km_cluster = kmeans_cluster(self._tsne_out, self.seed)
                ## k 입력이 없으면 k_range 가지고 optimal k를 찾음, k_range는 default 값이 있으므로 user input이 없어도 적용됨
                print_progress("K is not seleted. Find the optimal K between the range ({})".format(self.k_range))
                output_name = os.path.join(self.output_dir_path, "{}_choose_K_silhouette_score_plot.pdf".format(self.tag)) # silhouette score plot (dawn output 1)
                self._k_val = km_cluster.optimal_k(self.k_range, output_name)

                print_progress("Find the optimal K = {}".format(self._k_val))
                
            ## k 입력이 있는 경우 k_range를 무시하고 k를 사용함, k와 k_range 값이 둘 다 있으면 k가 우선됨
            print_arg("K for K-means clustering algorithm", self._k_val)

        return self._k_val      

    
    @property
    def seed(self):
        return self.args.seed

    @property
    def tag(self):
        return self.args.tag

    @property
    def category_set(self):
        category_set_ = pd.read_table(self.args.category_set_file)
        return category_set_

    @property
    def category_count(self):
        category_count_ = pd.read_table(self.args.category_count_file)
        return category_count_

    @property
    def count_threshold(self):
        return self.args.count_threshold

    @property
    def corr_threshold(self):
        return self.args.corr_threshold
    
    @property
    def size_threshold(self):
        return self.args.size_threshold

    def _get_input_file(self, infix: str) -> Path:
        input_file = glob.glob(os.path.join(self.input_dir_path, infix))

        if len(input_file) == 0:
            raise FileNotFoundError(
                "'{}' cannot be found. The input file named '{}' must exist in the input file directory. Check the directory of input files.".format(infix, infix)
            )
        if len(input_file) > 1:
            raise Exception(
                "Too many files with '{}'. There must be one file in the input directory with this name. Check the directory of input files.".format(infix)
            )

        return input_file[0]

    def run(self):
        self.tsne_projection()
        self.kmeans_clustering()
        self.dawn_analysis()
        print_progress("Done")

    def tsne_projection(self): ## t-SNE and k-means clustering by optimal K
        random.seed(self.seed)
        print_progress("t-SNE projection for {}".format(self.eig_vector_file))

        # self.eig_vector = pd.read_table(self._eig_vector, compression='gzip')

        if type(self.eig_vector.iloc[0,0]) == np.complex128:
            self.eig_vector = self.eig_vector.astype(np.float128)

        U_pick = self.eig_vector.iloc[:, 1:50]
        scaler = StandardScaler(with_mean=False, with_std=False)
        scale_factor = np.sqrt(np.sum(U_pick**2, axis=1))
        U_norm = scaler.fit_transform(U_pick.T / scale_factor).T

        tsne_res = TSNE(n_components=2, perplexity=30, n_iter=500, random_state=self.seed, init='pca', method='exact').fit_transform(U_norm)

        self._U_norm = U_norm
        self._tsne_out = pd.DataFrame(tsne_res, columns=['t-SNE1', 't-SNE2'])

    def kmeans_clustering(self):
        random.seed(self.seed)
        km_cluster = kmeans_cluster(self._tsne_out, self.seed)

        print_progress("K-means clustering with {} clusters.".format(self.k_val))
        i_init_pt = km_cluster.center_init(self.k_val)

        # initial centers of clusters are given
        fit = self.kmeans_r(self._U_norm, centers=self._U_norm[i_init_pt,], iter_max=100)
        fit_res = dict(zip(fit.names, list(fit)))
        fit_res['annotation'] = self.corr_mat.columns
        self._fit_res = fit_res

        category_tsne = pd.DataFrame({"category": self.corr_mat.columns,
                                      "cluster": fit_res['cluster'],
                                      "tsne1": self._tsne_out['t-SNE1'],
                                      "tsne2": self._tsne_out['t-SNE2']})        
        category_tsne['is_center'] = 0
        category_tsne.loc[category_tsne.index.isin(i_init_pt), 'is_center'] = 1
        
    def dawn_analysis(self):
        random.seed(self.seed)
        print_progress("DAWN analysis")
        # supernode_resDir = os.path.join(self.output_dir_path, "supernodeWGS_results", "blocks_"+self.tag) # correlation block (dawn output 3)
        clusters = list(self._fit_res['cluster'])

        supernodeWGS_process = supernodeWGS_func(self.corr_mat,
                                                 self._fit_res,
                                                 self.k_val,
                                                 self.num_proc,
                                                 self.output_dir_path,
                                                 self.tag,
                                                 self.seed,
                                                 verbose=True)
        
        ## create correlation matrix blocks
        max_val = int(self.k_val * (self.k_val-1)/2 + self.k_val)
        print_progress("[DAWN] Create correlation matrix blocks")
        for i in tqdm(range(1, max_val+1)):
            supernodeWGS_process.corr_mat_blocks(i)
        
        ## preparation for dawn
        # Load data, which contains the number of variants in each category.
        categories = self.category_set['Category'].tolist()
        category_count_sub = self.category_count.loc[self.category_count['Category'].isin(categories)]
        permut_test_sub = self.permut_test.loc[self.permut_test['Category'].isin(categories)]

        print_progress("[DAWN] Preprocess data for the DAWN analysis")
        cat_count_permut, cluster_idx, cluster_size = supernodeWGS_process.dawn_preprocess(category_count_sub,
                                                                                           permut_test_sub,
                                                                                           self.count_threshold,
                                                                                           self.size_threshold)

        print_progress("[DAWN] Compute graph on correlation matrix and form graph")
        form_data = data_collection(path=os.path.join(self.output_dir_path, "supernodeWGS_results", "blocks_"+self.tag),
                                    cores=30, max_cluster=self.k_val, verbose=True)
        cor_mat = form_data.form_correlation(k=cluster_idx)
        g = supernodeWGS_process.form_graph_from_correlation(cor_mat,
                                                             func=lambda x: x>self.corr_threshold,
                                                             k=cluster_idx)
        adj_mat = pd.DataFrame(np.array(g.get_adjacency().data), index=g.vs['name'], columns=g.vs['name'])
        adj_mat.to_csv(os.path.join(self.output_dir_path, "{}.ipvalue_fdr_igraph.csv".format(self.tag)), sep=",")
                
        print_progress("[DAWN] Compute the test statistics with z-scores")
        testvec_res = form_data.form_testvec(vec=cat_count_permut['P'],
                                             clustering=clusters,
                                             flag_vec=cat_count_permut['flag'],
                                             k=cluster_idx,
                                             sparse=True,
                                             sumabsv=5.25)
        zval_supernode = testvec_res['vec']
                
        print_progress("[DAWN] Compute the test statistics with relative risks")
        riskvec_res = form_data.form_testvec(vec=cat_count_permut['Relative_Risk'],
                                             clustering=clusters,
                                             flag_vec=cat_count_permut['flag2'],
                                             k=cluster_idx)
        tmp = riskvec_res['vec']
        tmp = tmp.copy()
        tmp.loc[np.isnan(tmp)] = 0
        risk_supernode = np.exp(tmp)
        
        ## apply dawn to clusters
        print_progress("[DAWN] Apply DAWN to clusters and calculate p-values and FDR")
        hmrf_res = supernodeWGS_process.hmrf(z=zval_supernode,
                                             adj=adj_mat,
                                             seedindex=[0 for x in cluster_idx],
                                             verbose=True)
        fdr = supernodeWGS_process.report_results(vec=cluster_idx,
                                                  posterior=1-hmrf_res['post'],
                                                  pvalue=1-norm.cdf(zval_supernode),
                                                  Iupdate=hmrf_res['Iupdate'])
        fdr.indicator = fdr.indicator.astype(int)
        fdr.to_csv(os.path.join(self.output_dir_path, "{}.ipvalue_fdr.txt".format(self.tag)), sep="\t", index=False)
        
        cluster_pval = np.empty(len(zval_supernode)) * np.nan
        cluster_pval[np.where(risk_supernode >= 1)[0]] = 1 - norm.cdf(zval_supernode[np.where(risk_supernode >= 1)[0]])
        cluster_pval[np.where(risk_supernode < 1)[0]] = norm.cdf(zval_supernode[np.where(risk_supernode < 1)[0]])

        mat = pd.DataFrame({"Cluster.idx":g.vs['name'], "Pval.cluster":cluster_pval, "Risk":risk_supernode})
        mat.to_csv(os.path.join(self.output_dir_path, "{}.ipvalue_fdr_ipvalue_risk.csv".format(self.tag)), sep=",", index=False)


        ## save dawn plots and tables
        print_progress("[DAWN] Draw DAWN clusters network graph and save graph layout")
        supernodeWGS_process.dawn_plot(g, zval_supernode)

        ## 
        print_progress("[DAWN] Save DAWN clusters with annotations")
        tmp_idx = list(range(g.vcount()))
        cluster_indices = np.array(g.vs['name'])[tmp_idx].astype(int)
        assert all(np.array(sorted(cluster_indices)) == np.array(sorted(cluster_idx[tmp_idx])))

        result = [[] for x in cluster_indices]

        for i in range(len(cluster_indices)):
            node_idx = testvec_res['index'][tmp_idx[i]]
            node_nam = self._fit_res['annotation'][node_idx]

            assert len(set(np.array(clusters)[node_idx])) == 1, "Stop if they don't belong in the same cluster"
            assert len(node_idx) <= cluster_size[cluster_indices[i]-1], "Stop if the number of categories in the cluster is bigger than the size of the cluster."

            result[i] = list(node_nam)

        max_len = max(list(map(len, result)))

        # Create a matrix where the result will be stored.
        mat = pd.DataFrame(index=range(max_len), columns=range(len(cluster_indices)))

        for i in range(len(cluster_indices)):
            mat.iloc[:len(result[i]), i] = result[i]
        
        mat.columns = cluster_indices
        mat.to_csv(os.path.join(self.output_dir_path, "{}.cluster_annotation.csv".format(self.tag)), sep=",", index=False)