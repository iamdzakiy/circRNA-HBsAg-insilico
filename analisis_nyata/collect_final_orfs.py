import json
M = {'V14': 'SA_11-11_0_19_c85', 'V15': 'SA_16-16_0_19_c85', 'V16': 'SA_3-2GCT_20_19_c85', 'V17': 'SA_3-2GCT_20_0_c85', 'V18': 'SPLIT_403_sa_c85'}
json.dump({v: json.load(open(f'final/{n}.json'))['orf'] for v, n in M.items()}, open('final_orfs.json', 'w'))
print('final_orfs ok')
