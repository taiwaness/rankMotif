from ranking import PatternScoring
from matrix import position_frequency_matrix as pfm
from matrix import position_weight_matrix as pwm


class Cluster(object):

    def __init__(self, max_cluster=5, similarity=0.8):
        self.max_cluster = max_cluster
        self.similarity = similarity
        self.clustered_patterns = {}

    def run(self, pattern_scoring):
        isinstance(pattern_scoring, PatternScoring)

        self.clustered_patterns = {}
        ranked_patterns = sorted(pattern_scoring.results, key=lambda x: pattern_scoring.results.get(x), reverse=True)
        clusters = [0] * len(pattern_scoring.results)
        num_cluster = 0

        for i in len(pattern_scoring.results):
            if num_cluster == self.max_cluster:
                break
            if clusters[i] == 0:
                num_cluster += 1
                clusters[i] = num_cluster
                self.clustered_patterns.update({num_cluster: [ranked_patterns[i]]})
                for j in range(i + 1, len(self)):
                    if pwm(pfm(ranked_patterns[i]), pfm(ranked_patterns[j])) >= self.similarity:
                        clusters[j] = num_cluster
                        self.clustered_patterns.get(i).append(ranked_patterns[j])
