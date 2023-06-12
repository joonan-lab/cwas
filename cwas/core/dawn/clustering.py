from sklearn.metrics import silhouette_score
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pandas2ri.activate()

class kmeans_cluster:
    def __init__(self, tsne_out: pd.DataFrame, seed: int) -> None:
                 #output_name: str, k_range: str, k=None) -> None:
        
        self._tsne_out = tsne_out
        self._seed = seed
        self.kmeans_r = importr('stats').kmeans
        self.set_seed = robjects.r('set.seed')

    @property
    def tsne_out(self) -> pd.DataFrame:
        return self._tsne_out
    
    @property
    def seed(self) -> str:
        return self._seed
    
    def optimal_k(self, k_range, output_name):
        start, end = list(map(int, k_range.replace(" ", "").split(",")))
        k_values = list(range(start, end+1))

        self.set_seed(self.seed)

        avg_sil_values = list(map(self._avg_sil, k_values))
        df = pd.DataFrame({"k_values": k_values, "avg_sil_values": avg_sil_values})
        self.avg_sil_df = df

        self._print_silhouett_plot(output_name)
        opt_k = df.loc[df.avg_sil_values==max(df.avg_sil_values), "k_values"].values[0]

        return opt_k
    
    def _avg_sil(self, k):
        pandas2ri.activate()
        km_res = self.kmeans_r(np.array(self.tsne_out), centers=k, nstart=300, iter_max=100)
        km_res_dict = dict(zip(km_res.names, list(km_res)))
        distance_matrix = squareform(pdist(self.tsne_out))
        ss = silhouette_score(distance_matrix, km_res_dict['cluster'], metric='precomputed')

        return ss
    
    def _print_silhouett_plot(self, output_name):
        k_start = min(self.avg_sil_df['k_values'])
        k_end = max(self.avg_sil_df['k_values'])
        max_ss = max(self.avg_sil_df.avg_sil_values)
        k = self.avg_sil_df.loc[self.avg_sil_df.avg_sil_values==max_ss, "k_values"].values[0]
        
        plt.figure(figsize=(20, 10))
        plt.title('Range of K: {} to {}'.format(k_start, k_end), fontsize=15, pad=10)
        plt.plot(self.avg_sil_df['k_values'], self.avg_sil_df['avg_sil_values'], '-o', c="black", linewidth=1, markersize=4)
        plt.vlines(k, 0, max_ss, color='red', linestyles='--')
        plt.text(k, max_ss+0.01, "Optimal K = {}".format(str(k)), fontdict={'color': 'red', 'weight': 'bold', 'fontsize': 15})
        plt.xticks(np.arange(0, int(np.ceil(k_end/10)*10)+1, 10), fontsize=10)
        plt.yticks(np.arange(0, max_ss+0.1, 0.1), fontsize=10)
        plt.xlim(k_start-5, k_end+5)
        plt.ylim(0, max_ss+0.1)
        plt.xlabel("Number of clusters K", fontsize=12, labelpad=10)
        plt.ylabel("Average Silhouettes", fontsize=12, labelpad=10)
        plt.tight_layout()
        plt.savefig(output_name) # final output
        #plt.savefig(self._output_name)


    def center_init(self, k):
        self.set_seed(self.seed)

        km_tsne = self.kmeans_r(self.tsne_out, centers=k, nstart=300, iter_max=100)
        km_tsne_dict = dict(zip(km_tsne.names, list(km_tsne)))

        km_tsne_cluster = [x - 1 for x in list(km_tsne_dict['cluster'])]
        km_tsne_centers = pd.DataFrame(km_tsne_dict['centers'][km_tsne_cluster], columns=["t-SNE1", "t-SNE2"])
        dist_to_center_tsne = np.sqrt(np.sum((self.tsne_out - km_tsne_centers)**2, axis=1))

        self.km_tsne_cluster_ = km_tsne_cluster
        self.dist_to_center_tsne_ = dist_to_center_tsne
        
        i_init_pt = list(map(self._init_pt_kmeans, list(range(k))))

        return i_init_pt

    def _init_pt_kmeans(self, k):
        i_check = list(np.where([x == k for x in self.km_tsne_cluster_])[0])
        i_check_min = self.dist_to_center_tsne_[i_check].idxmin()

        return i_check_min




