import json
d = json.load(open('../build-1/result_qb1.json'))
for row in d['routes']:
    print(row['x'][0][:10], row['chart'], row['source'][:30], [k for k in row if k.startswith('overlap')][:6], all(v for k, v in row.items() if k.startswith('overlap')), str(row.get('einstein_contains_0', row.get('ein', '')))[:20])
print(list(d['routes'][0].keys()))
