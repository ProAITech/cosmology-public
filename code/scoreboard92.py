"""Entry 93: the scoreboard after 92 entries, counted from the Gap Map's nodes."""
import re, json, os
from collections import Counter
os.chdir(os.path.dirname(os.path.abspath(__file__)))
g = open('../gapmap.html').read()
nodes = []
for m in re.finditer(r"\{id:'([^']+)',theme:'([^']+)',label:'((?:[^'\\]|\\.)*)',status:'([^']+)',aff:'([^']+)'", g):
    s = g[m.start():g.find('\n', m.start())]
    links = re.findall(r"\['Workbench · Entry (\d+)'", s)
    nodes.append(dict(id=m.group(1), label=m.group(3).replace("\\'", "'"), status=m.group(4), aff=m.group(5), entries=sorted(set(int(x) for x in links))))
order = ['settled', 'single', 'open', 'contested', 'notfound', 'excluded']
names = {'settled': 'settled', 'single': 'measured once', 'open': 'open', 'contested': 'contested', 'notfound': 'not found', 'excluded': 'excluded'}
groups = [('P', 'Plasma-cosmology claims'), ('S', 'Shared physics'), ('L', 'Standard-model pieces'), ('A', 'Changed gravity')]
rows = []
for a, lab in groups:
    c = Counter(n['status'] for n in nodes if n['aff'] == a)
    rows.append([lab + ' (%d)' % sum(c.values()), [[names[k], c[k], k] for k in order if c[k]]])
res = dict(n_nodes=len(nodes), rows=rows, by_side={a: dict(Counter(n['status'] for n in nodes if n['aff'] == a)) for a, _ in groups},
           nodes=nodes,
           entry50=dict(n_nodes=38, plasma=dict(excluded=6, contested=2, notfound=1)))
json.dump(res, open('scoreboard92.json', 'w'), indent=1, ensure_ascii=False)
print(res['n_nodes'], res['by_side'])
for n in nodes: print(n['aff'], n['status'].ljust(9), n['id'].ljust(9), n['label'][:45], n['entries'])
