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
import zarr

from cwas.core.dawn.clustering import kmeans_cluster
from cwas.core.dawn.supernodeWGS import supernodeWGS_func, data_collection
from cwas.runnable import Runnable
from cwas.utils.check import check_is_file, check_is_dir, check_num_proc
from cwas.utils.log import print_arg, print_progress

class Dawn(Runnable):
    def __init__(self, args: argparse.Namespace):
        super().__init__(args)
        self._eig_vector = None
        self._corr_mat = None
        self._permut_test = None
        self._category_set = None
        self._k_val = None
        self.kmeans_r = importr('stats').kmeans

    @staticmethod
    def _print_args(args: argparse.Namespace):
        print_arg(
            "No. worker processes for the DAWN",
            f"{args.num_proc: ,d}",
        )
        print_arg("Eigen vector file", args.eig_vector_file)
        print_arg("Correlation matrix file", args.corr_mat_file)
        print_arg("Permutation test file", args.permut_test_file)
        print_arg("Category variant (or sample) counts file (burden test result): ", args.category_count_file)
        print_arg("Output directory", args.output_dir_path)
        print_arg("Lambda value:", args.lambda_val)
        print_arg("Thresholds of count / correlation / size: ", ", ".join(list(map(str, [args.count_threshold, args.corr_threshold, args.size_threshold]))))

    @staticmethod
    def _check_args_validity(args: argparse.Namespace):
        check_num_proc(args.num_proc)
        check_is_dir(args.eig_vector_file)
        check_is_dir(args.corr_mat_file)
        check_is_file(args.permut_test_file)
        check_is_file(args.category_count_file)
        check_is_dir(args.output_dir_path)

    @property
    def num_proc(self):
        return self.args.num_proc

    @property
    def lambda_val(self) -> float:
        return self.args.lambda_val

    @property
    def input_dir_path(self) -> Path:
        return self.args.input_dir_path
    
    @property
    def eig_vector_file(self):
        return self.args.eig_vector_file
    
    @property
    def eig_vector(self):
        if self._eig_vector is None:
            root = zarr.open(self.eig_vector_file, mode='r')
            self._eig_vector = pd.DataFrame(data=root['data'],
                              index=root['metadata'].attrs['category'])
            #self._eig_vector = self._eig_vector.iloc[:,1:]
        return self._eig_vector
    
    @property
    def corr_mat_file(self):
        return self.args.corr_mat_file

    @property
    def corr_mat(self):
        if self._corr_mat is None:
            root = zarr.open(self.corr_mat_file, mode='r')
            self._corr_mat = root['data']
            column_indices = list(map(lambda x: root['metadata'].attrs['category'].index(x), self.category_set))
            self._corr_mat = self._corr_mat[column_indices][:, column_indices].astype(np.float64)
            #self._corr_mat = pd.DataFrame(self._corr_mat, index=self.category_set, columns=self.category_set)
            #self._corr_mat = self._corr_mat.loc[self.category_set, self.category_set].astype(np.float64)
        return self._corr_mat

    @property
    def permut_test_file(self):
        return self.args.permut_test_file

    @property
    def permut_test(self):
        if self._permut_test is None:
            self._permut_test = pd.read_table(self.permut_test_file, compression='gzip', sep="\t")
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
                print_progress("K is not selected. Find the optimal K between the range ({})".format(self.k_range))
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
    def tsne_method(self):
        return self.args.tsne_method

    @property
    def tag(self) -> str:
        return str(self.args.tag)

    @property
    def category_set(self):
        if self._category_set is None:
            self._category_set = self.eig_vector.index.tolist()
        return self._category_set

    @property
    def category_count(self):
        category_count_ = pd.read_table(self.args.category_count_file, sep="\t")
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

    def run(self):
        self.tsne_projection()
        self.kmeans_clustering()
        self.dawn_analysis()
        print_progress("Done")

    def tsne_projection(self): ## t-SNE and k-means clustering by optimal K
        random.seed(self.seed)
        print_progress("t-SNE projection for {}".format(self.eig_vector_file))

        # self.eig_vector = pd.read_table(self._eig_vector, compression='gzip')

        if type(self.eig_vector.iloc[0,0]) != np.float64:
            self.eig_vector = self.eig_vector.astype(np.float64)

        U_pick = self.eig_vector.iloc[:, 1:50]
        scaler = StandardScaler(with_mean=False, with_std=False)
        scale_factor = np.sqrt(np.sum(U_pick**2, axis=1))
        U_norm = scaler.fit_transform(U_pick.T / scale_factor).T

        tsne_res = TSNE(n_components=2, 
                        perplexity=30,
                        n_iter=500,
                        random_state=self.seed,
                        init='pca',
                        method=self.tsne_method,
                        n_jobs=self.num_proc).fit_transform(U_norm)

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
        fit_res['annotation'] = self.category_set
        self._fit_res = fit_res
        
    def dawn_analysis(self):
        random.seed(self.seed)
        print_progress("DAWN analysis")
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
        category_count_sub = self.category_count.loc[self.category_count['Category'].isin(self.category_set)]
        permut_test_sub = self.permut_test.loc[self.permut_test['Category'].isin(self.category_set)]

        print_progress("[DAWN] Preprocess data for the DAWN analysis")
        cat_count_permut, cluster_idx, cluster_size = supernodeWGS_process.dawn_preprocess(category_count_sub,
                                                                                           permut_test_sub,
                                                                                           self.count_threshold,
                                                                                           self.size_threshold)
        print_progress("[DAWN] Compute graph on correlation matrix and form graph")
        form_data = data_collection(path=os.path.join(self.output_dir_path, "supernodeWGS_results", "blocks_"+self.tag),
                                    cores=self.num_proc,
                                    max_cluster=self.k_val,
                                    verbose=True)
        cor_mat = form_data.form_correlation(k=cluster_idx)
        g = supernodeWGS_process.form_graph_from_correlation(cor_mat,
                                                             func=lambda x: x>self.corr_threshold,
                                                             k=cluster_idx)
        adj_mat = pd.DataFrame(np.array(g.get_adjacency().data), index=g.vs['name'], columns=g.vs['name'])
        adj_mat.to_csv(os.path.join(self.output_dir_path, "{}.adjacency_matrix.csv".format(self.tag)), sep=",")
        
        print_progress("[DAWN] Compute the test statistics with z-scores")
        testvec_res = form_data.form_testvec(vec=cat_count_permut['P'],
                                             clustering=clusters,
                                             flag_vec=cat_count_permut['flag'],
                                             k=cluster_idx,
                                             sparse=True,
                                             sumabsv=self.lambda_val)
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
        #fdr.to_csv(os.path.join(self.output_dir_path, "{}.ipvalue_fdr.txt".format(self.tag)), sep="\t", index=False)
        
        cluster_pval = np.empty(len(zval_supernode)) * np.nan
        cluster_pval[np.where(risk_supernode >= 1)[0]] = 1 - norm.cdf(zval_supernode[np.where(risk_supernode >= 1)[0]])
        cluster_pval[np.where(risk_supernode < 1)[0]] = norm.cdf(zval_supernode[np.where(risk_supernode < 1)[0]])

        mat = pd.DataFrame({"Cluster.idx":g.vs['name'], "Pval.cluster":cluster_pval, "Risk":risk_supernode})
        mat.to_csv(os.path.join(self.output_dir_path, "{}.ipvalue_fdr_ipvalue_risk.csv".format(self.tag)), sep=",", index=False)

        ## save dawn plots and tables
        ### annotation table
        print_progress("[DAWN] Save DAWN clusters with annotations")
        tmp_idx = list(range(g.vcount()))
        cluster_indices = np.array(g.vs['name'])[tmp_idx].astype(int)
        assert all(np.array(sorted(cluster_indices)) == np.array(sorted(cluster_idx[tmp_idx])))

        annot = [[] for x in cluster_indices]

        for i in range(len(cluster_indices)):
            node_idx = testvec_res['index'][tmp_idx[i]]
            node_nam = np.array(self._fit_res['annotation'])[node_idx]

            assert len(set(np.array(clusters)[node_idx])) == 1, "Stop if they don't belong in the same cluster"
            assert len(node_idx) <= cluster_size[cluster_indices[i]-1], "Stop if the number of categories in the cluster is bigger than the size of the cluster."

            annot[i] = list(node_nam)

        max_len = max(list(map(len, annot)))

        # Create a matrix where the result will be stored.
        annot_mat = pd.DataFrame(index=range(max_len), columns=range(len(cluster_indices)))

        for i in range(len(cluster_indices)):
            annot_mat.iloc[:len(annot[i]), i] = annot[i]
        
        annot_mat.columns = cluster_indices
        annot_mat.to_csv(os.path.join(self.output_dir_path, "{}.cluster_annotation.csv".format(self.tag)), sep=",", index=False)
        
        ## dawn plots
        print_progress("[DAWN] Draw DAWN clusters network graph and save graph layout")
        supernodeWGS_process.dawn_plot(g, zval_supernode, annot)
