import pandas as pd

scores = pd.read_csv("sandbox/chi_scores", sep='\t', header=None, index_col=0)
best = pd.DataFrame([list(scores.sort_values(x+1).index) for x in range(scores.shape[1])])
best = best.applymap(lambda x: x.split('/')[-1])
best.transpose().to_csv('sandbox/sorted', header=None, sep='\t')
