import logging
from .scoring import PatternScoring
from .seqio import gc_content
from .pfm import pfm, simpfm


class Cluster(object):

    def __init__(self, max_cluster=5, similarity=0.8,
                 max_patterns_per_cluster=5, simpfm_max_wsize='auto', reverse_complement=False):
        self._logger = logging.getLogger(self.__class__.__name__)
        self.max_cluster = max_cluster
        self.similarity = similarity
        self.max_patterns_per_cluster = max_patterns_per_cluster
        self.simpfm_max_wsize = simpfm_max_wsize
        self.reverse_complement = reverse_complement
        self.results = {}

    def run(self, pattern_scoring, gc=None, pset=None):
        isinstance(pattern_scoring, PatternScoring)
        self._logger.info('begin running pattern clustering')

        self.results = {}

        if not gc:
            gc = gc_content(pset)
        ranked_patterns = sorted(
            pattern_scoring.results, key=lambda x: pattern_scoring.results.get(x), reverse=True)
        clusters = [0] * len(pattern_scoring.results)
        n_clusters = 0
        n_patterns = {}

        for i in xrange(len(ranked_patterns)):
            if n_clusters == self.max_cluster:
                break
            if clusters[i] == 0:
                n_clusters += 1
                clusters[i] = n_clusters
                n_patterns.update({n_clusters: 0})
                self.results.update({n_clusters: [ranked_patterns[i]]})
                for j in xrange(i + 1, len(ranked_patterns)):
                    if n_patterns.get(n_clusters) == self.max_patterns_per_cluster - 1:
                        break
                    score = simpfm(
                        pfm(ranked_patterns[i]), pfm(ranked_patterns[j]),
                        gc, self.simpfm_max_wsize, self.reverse_complement)[2]
                    if clusters[j] == 0 and score >= self.similarity:
                        clusters[j] = n_clusters
                        self.results.get(n_clusters).append(ranked_patterns[j])
                        n_patterns[n_clusters] += 1

        return self
