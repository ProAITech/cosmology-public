"""Entry 143: the scoreboard after 142 entries, counted from the Gap Map's nodes (hub nodes left out, as in Entry 93), and what changed."""
import re, json, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
g = open(os.path.join(HERE, '..', 'gapmap.html')).read()
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
old = json.load(open(os.path.join(HERE, 'scoreboard126.json')))
oldn = {n['id']: n for n in old['nodes']}
new_nodes = [n for n in nodes if n['id'] not in oldn]
changed = [(n, oldn[n['id']]['status']) for n in nodes if n['id'] in oldn and oldn[n['id']]['status'] != n['status']]
gone = [n for i, n in oldn.items() if i not in {x['id'] for x in nodes}]
res = dict(n_nodes=len(nodes), rows=rows, by_side={a: dict(Counter(n['status'] for n in nodes if n['aff'] == a)) for a, _ in groups},
           nodes=nodes, new=[(n['id'], n['aff'], n['status'], n['label']) for n in new_nodes],
           changed=[(n['id'], n['aff'], o, n['status'], n['label']) for n, o in changed], gone=[(n['id'], n['aff'], n['status'], n['label']) for n in gone],
           entry127=old['by_side'])
json.dump(res, open(os.path.join(HERE, 'scoreboard142.json'), 'w'), indent=1, ensure_ascii=False)
print(res['n_nodes'], res['by_side']); print('entry 127:', old['n_nodes'], old['by_side'])
print('new:'); [print('  ', x) for x in res['new']]
print('changed:'); [print('  ', x) for x in res['changed']]
print('gone:'); [print('  ', x) for x in res['gone']]
