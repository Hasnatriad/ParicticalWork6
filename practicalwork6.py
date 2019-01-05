import networkx as nx
import random
from collections import OrderedDict
import matplotlib.pylab as plt

G = nx.read_edgelist('Flickr-test.txt')
G = nx.convert_node_labels_to_integers(G, first_label=0, ordering='default')

#Excercise 1
def degree(G):
    deg = [0 for x in range(len(G))]

    j = 0
    for i in G:
        deg[j] = len(G[i])
        j += 1
    return deg

def degdis(G, deg):
    i = (max(deg))
    dd = [0 for x in range(i)]
    for x in range(i):
        dd[x] = deg.count(x+1)
    return dd

def avdeg(deg):
    s = 0
    for i in range(len(deg)):
        s = s + deg[i]
    avg = s/len(deg)
    return avg

def tri(G, n, deg):
    tr = [0] * n
    u = 0
    while u < n:
        i = 0
        while i < deg[u]:
            v = list(G.neighbors(u))
            v = v[i]
            if u < v:
                U = list(G.neighbors(u))
                V = list(G.neighbors(v))
                iu = 0
                iv = 0
                while (iu < deg[u]) & (iv < deg[v])\
                        & (U[iu] < u) & (V[iv] < u):
                    if U[iu] < V[iv]:
                        iu += 1
                    else:
                        if U[iu] > V[iv]:
                            iv += 1
                        else:
                            tr[u-1] += 1
                            tr[v-1] += 1
                            tr[U[iu]] += 1
                            iu += 1
                            iv += 1
            i += 1
        u += 1
    return tr

def triplet(G, n):
    i = 0
    trip = [0] * n
    while i < n:
        neigh = list(G.neighbors(i))
        len_neigh = len(neigh)
        if len_neigh >= 2:
            trip[i] = len_neigh*(len_neigh-1)/2
        i += 1
    return trip

def ccoef(tr, trip, n):
    cc = 0
    total = 0
    i = 0
    while i < n:
        if trip[i] != 0:
            total += tr[i]/trip[i]
        i += 1
    cc = total/n
    return cc

def trratio(tr, trip, n):
    i = 0
    sum_tr = 0
    sum_trip = 0
    while i < n:
        sum_tr += tr[i]
        sum_trip += trip[i]
        i += 1
    tranratio = (3 * sum_tr) / sum_trip
    return tranratio

n = len(G.nodes())
print('No. of nodes ', n)

m = len(G.edges())
print('No. of edges ', m)

den = (2*m)/(n*(n-1))
print('Density ', den)

deg = degree(G)
print("Degree ", deg)

maxdeg = max(deg)
print("Maximum Degree ", maxdeg)

dd = degdis(G, deg)
print("Degree Distribution ", dd)
plt.loglog(dd)
plt.grid(True, which="both")
plt.title('Degree Distribution - Graph Flickr-test')
plt.savefig("Flickr-test DegDist.png")
plt.show()

avg = avdeg(deg)
print("Average Degree ", avg)

tr = tri(G, n, deg)
print("No. of Triangles ", tr)

trip = triplet(G, n)
print("No. of Triplets ", trip)

cc = ccoef(tr, trip,  n)
print("Clustering Coefficient ", cc)

tratio = trratio(tr, trip, n)
print("Transitive Ratio ", tratio)

#Excercise 2
def emptygraph(G):
    sample = G.copy()
    sample = nx.create_empty_copy(sample, with_data = True)
    return sample
sample = emptygraph(G)

#Excercise 3
def testlink(G, sample, u, v, no_tests, no_links_discovered, islink, fileres):
    no_tests += 1
    if G.has_edge(u, v):
        sample.add_edge(u, v)
        fileres.append((no_tests, u, v))   #Excercise 4
        no_links_discovered += 1
        islink = True
    
    return sample, no_tests, no_links_discovered, islink, fileres

#Excercise 5
def analyze(n, m, no_tests):    
    effbest = min(no_tests, m)
    effworst = max(0, no_tests - (n*(n-1)/2-m))    
    return effbest, effworst

#Excercise 6
def eff(n, m, no_tests, effbest, effworst, abeff, normeff, releff):
    maxeff, mineff = analyze(n, m, no_tests)
    if (maxeff - mineff) > 0:
        normeff = (abeff - mineff)/(maxeff - mineff)
    return normeff

no_tests_ran = 0
no_tests_comp = 0
no_tests_tbf = 0
no_tests_mix = 0
no_links_discovered_ran = 0
no_links_discovered_comp = 0
no_links_discovered_tbf = 0
no_links_discovered_mix = 0
no_links_discovered_before_ran = no_links_discovered_ran
no_links_discovered_before_comp = no_links_discovered_comp
no_links_discovered_before_tbf = no_links_discovered_tbf
no_links_discovered_before_mix = no_links_discovered_mix
fileres_ran = []
fileres_comp = []
fileres_tbf = []
fileres_mix = []
islink = False
pairstested_ran = []
pairstested_comp = []
pairstested_tbf = []
pairstested_mix = []
effbest_ran = []
effbest_comp = []
effbest_tbf = []
effbest_mix = []
effworst_ran = []
effworst_comp = []
effworst_tbf = []
effworst_mix = []
raneff_ran = 0
raneff_comp = 0
raneff_tbf = 0
raneff_mix = 0
abeff_ran = 0
abeff_comp = 0
abeff_tbf = 0
abeff_mix = 0
normeff_ran = 0
normeff_comp = 0
normeff_tbf = 0
normeff_mix = 0
releff_ran = 0
releff_comp = 0
releff_tbf = 0
releff_mix = 0
evolution_ran = []
evolution_comp = []
evolution_tbf = []
evolution_mix = []

#Excercise 7
def randomstrategy(G, sample, n, m, no_tests, no_links_discovered, no_links_discovered_before, fileres, islink, pairstested, effbest, effworst, raneff, abeff, normeff, releff, evolution):
    for x in range(50000):
        u = random.choice(list(G.nodes()))
        v = random.choice(list(G.nodes()))
        while u == v:
            v = random.choice(list(G.nodes()))
        islink = False
        if sample.has_edge(u, v) or (u, v) in pairstested:
            pass
        else:
            sample, no_tests, no_links_discovered, islink, fileres = testlink(G, sample, u, v, no_tests, no_links_discovered, islink, fileres)
            if no_links_discovered > no_links_discovered_before:
                abeff = abeff + no_links_discovered
                effbest.append(min(no_tests, no_links_discovered))
                effworst.append(max(0, no_tests - (n*(n-1)/2-no_links_discovered)))
            no_links_discovered_before = no_links_discovered
            pairstested.append((u, v))
            pairstested.append((v, u))
            evolution.append((no_tests, no_links_discovered))

        if islink:
            U = list(G.neighbors(u))
            V = list(G.neighbors(v))
            j = 0
            k = 0
            while j < len(U):
                if (sample.has_edge(v, U[j])) or ((v, U[j]) in pairstested) or (v == U[j]):
                    pass
                else:
                    sample, no_tests, no_links_discovered, islink, fileres = testlink(G, sample, v, U[j], no_tests, no_links_discovered, islink, fileres)
                    if no_links_discovered > no_links_discovered_before:
                        abeff = abeff + no_links_discovered
                        effbest.append(min(no_tests, no_links_discovered))
                        effworst.append(max(0, no_tests - (n*(n-1)/2-no_links_discovered))) 
                    no_links_discovered_before = no_links_discovered
                    pairstested.append((v, U[j]))
                    pairstested.append((U[j], v))
                    evolution.append((no_tests, no_links_discovered))
                j += 1
            while k < len(V):
                if (sample.has_edge(u, V[k])) or ((u, V[k]) in pairstested) or (u == V[k]):
                    pass
                else:
                    sample, no_tests, no_links_discovered, islink, fileres = testlink(G, sample, u, V[k], no_tests, no_links_discovered, islink, fileres)
                    if no_links_discovered >no_links_discovered_before:
                        abeff = abeff + no_links_discovered
                        effbest.append(min(no_tests, no_links_discovered))
                        effworst.append(max(0, no_tests - (n*(n-1)/2-no_links_discovered))) 
                    no_links_discovered_before = no_links_discovered
                    pairstested.append((u, V[k]))
                    pairstested.append((V[k], u))
                    evolution.append((no_tests, no_links_discovered))
                k += 1
    normeff = eff(n, no_tests, no_links_discovered, effbest, effworst, abeff, normeff, releff)
    raneff = normeff
    releff = normeff/raneff
    return sample, no_tests, no_links_discovered, no_links_discovered_before, raneff, abeff, normeff, releff, evolution

sample_ran, no_tests_ran, no_links_discovered_ran, no_links_discovered_before_ran, raneff_ran, abeff_ran, normeff_ran, releff_ran, evolution_ran = randomstrategy(G, sample, n, m, no_tests_ran, no_links_discovered_ran, no_links_discovered_before_ran, fileres_ran, islink, pairstested_ran, effbest_ran, effworst_ran, raneff_ran, abeff_ran, normeff_ran, releff_ran, evolution_ran)

print('No. of tests in Random Strategy ', no_tests_ran)
print('No. of links detected in Random Strategy ', no_links_discovered_ran)
print("Absolute Efficiency of Random Strategy ", abeff_ran)
print("Normalised Efficiency of Random Strategy ", normeff_ran)
print("Relative Efficiency of Random Strategy ", releff_ran)

m_ran = len(sample_ran.edges())
print('No. of links ', m_ran)

den = (2*m_ran)/(n*(n-1))
print('Density ', den)

deg = degree(sample_ran)
#print("Degree of the Graph is: ", deg)

maxdeg = max(deg)
print("Maximum Degree ", maxdeg)

dd = degdis(sample_ran, deg)

avg = avdeg(deg)
print("Average Degree is: ", avg)

tr = tri(sample_ran, n, deg)
#print("No. of Triangles ", tr)

trip = triplet(G, n)
#print("No. of Triplets ", trip)

cc = ccoef(tr, trip,  n)
print("CC ", cc)

tratio = trratio(tr, trip, n)
print("TR ", tratio)

sample_tbf = sample_ran.copy()
no_tests_tbf = no_tests_ran
no_links_discovered_tbf = no_links_discovered_ran
no_links_discovered_before_tbf = no_links_discovered_before_ran
fileres_tbf = fileres_ran.copy()
pairstested_tbf = pairstested_ran.copy()
effbest_tbf = effbest_ran.copy()
effworst_tbf = effworst_ran.copy()
raneff_tbf = raneff_ran
abeff_tbf = abeff_ran
normeff_tbf = normeff_ran
releff_tbf = releff_ran
evolution_tbf = evolution_ran.copy()

#Excercise 8
x, y = zip(*evolution_ran)
plt.plot(x,y)
plt.title('Random Strategy')
plt.xlabel('No. of tests')
plt.ylabel('No. of discovered links')
plt.grid()
plt.savefig("Random_Strategy.png")
plt.show()

#Excercise 10
def nodesordered(sample):
    nodes = sorted(sample.degree(), key=lambda x: x[1], reverse=True)
    return nodes
nodes = nodesordered(sample_ran)

#Excercise 11
def comp(G, sample, nodes, n, m, no_tests, no_links_discovered, no_links_discovered_before, fileres, islink, pairstested, effbest, effworst, raneff, abeff, normeff, releff, evolution):
    X = []
    alfa = 10
    for j in range(len(nodes)):
        if nodes[j][1] >= alfa:
            X.append(nodes[j][0])
    while len(X) > 0:
        u = X[0]
        X.remove(u)
        for i in range(len(nodes)):
            v = nodes[i][0]
            if (sample.has_edge(u,v)) or ((u, v) in pairstested) or (u == v):
                pass
            else:
                sample, no_tests, no_links_discovered, islink, fileres = testlink(G, sample, u, v, no_tests, no_links_discovered, islink, fileres)
                if no_links_discovered >no_links_discovered_before:
                    abeff = abeff + no_links_discovered
                    effbest.append(min(no_tests, no_links_discovered))
                    effworst.append(max(0, no_tests - (n*(n-1)/2-no_links_discovered))) 
                no_links_discovered_before = no_links_discovered
                pairstested.append((u, v))
                pairstested.append((v, u))
                evolution.append((no_tests, no_links_discovered))
                
            if (sample.has_edge(u,v)) & (len(sample.edges(v))== 1) & (sample.degree(v) > alfa):
                X.append(v)
    normeff = eff(n, no_tests, no_links_discovered, effbest, effworst, abeff, normeff, releff)
    raneff = normeff
    releff = normeff/raneff
    return sample, no_tests, no_links_discovered, no_links_discovered_before, abeff, normeff, releff, evolution

sample_comp, no_tests_comp, no_links_discovered_comp, no_links_discovered_before_comp, abeff_comp, normeff_comp, releff_comp, evolution_comp = comp(G, sample, nodes, n, m, no_tests_ran, no_links_discovered_ran, no_links_discovered_before_ran, fileres_ran, islink, pairstested_ran, effbest_ran, effworst_ran, raneff_ran, abeff_ran, normeff_ran, releff_ran, evolution_ran)

print('No. of tests in Complete Strategy ', no_tests_comp)
print('No. of links detected in Complete Strategy ', no_links_discovered_comp)
#print('Test No. for each edge found: ', fileres)
   
print("Absolute Efficiency for Complete Strategy ", abeff_comp)
print("Normalised Efficiency for Complete Strategy ", normeff_comp)
print("Relative Efficiency for Complete Strategy ", releff_comp)

m_comp = len(sample_comp.edges())
print('No. of edges ', m_comp)

den = (2*m_comp)/(n*(n-1))
print('Density of the Graph is: ', den)

deg = degree(sample_comp)
#print("Degree ", deg)

maxdeg = max(deg)
print("Maximum Degree of the graph is: ", maxdeg)

dd = degdis(sample_comp, deg)
#print("Degree Distribution ", dd)

avg = avdeg(deg)
print("Average Degree ", avg)

tr = tri(sample_comp, n, deg)
#print("No. of Triangles ", tr)

trip = triplet(sample_comp, n)
#print("No. of Triplets ", trip)

cc = ccoef(tr, trip,  n)
print("Clustering Coefficient ", cc)

tratio = trratio(tr, trip, n)
print("Transitive Ratio  ", tratio)


x, y = zip(*evolution_comp)
plt.plot(x,y)
plt.title('Complete Strategy')
plt.xlabel('No. of tests')
plt.ylabel('No. of discovered links')
plt.grid()
plt.savefig("Complete_Strategy.png")
plt.show()

#Excercise 12
def linksordered(sample):
    links = list(sample.edges())
    c = [0 for x in range(len(links))]
    keys = [0 for x in range(len(links))]
    values = [0 for x in range(len(links))]
    for i in range(len(links)):
        u = links[i][0]
        v = links[i][1]
        a = sample.degree(u)
        b = sample.degree(v)
        if (a >= 1) & (b >= 1):
            c[i] = a + b
            keys[i] = c[i]
            values[i] = links[i]
    
    links = dict(zip(keys, values))
    
    links = OrderedDict(sorted(links.items(), key=lambda t: t[0], reverse = True))
    
    return links

links = linksordered(sample_tbf)

#Excercise 13
def tbf(G, sample, links, n, m, no_tests, no_links_discovered, no_links_discovered_before, fileres, islink, pairstested, effbest, effworst, raneff, abeff, normeff, releff, evolution):
    keys = list(links.keys())
    for i in range(len(links)):
        u = links[keys[i]][0]
        for j in range(len(links)):
            v = links[keys[j]][1]
            if (sample.has_edge(u, v)) or ((u, v) in pairstested) or (u == v):
                pass
            else:
                sample, no_tests, no_links_discovered, islink, fileres = testlink(G, sample, u, v, no_tests, no_links_discovered, islink, fileres)
                if no_links_discovered >no_links_discovered_before:
                    abeff = abeff + no_links_discovered
                    effbest.append(min(no_tests, no_links_discovered))
                    effworst.append(max(0, no_tests - (n*(n-1)/2-no_links_discovered)))
                no_links_discovered_before = no_links_discovered
                pairstested.append((u, v))
                pairstested.append((v, u))
                evolution.append((no_tests, no_links_discovered))
    
    normeff = eff(n, no_tests, no_links_discovered, effbest, effworst, abeff, normeff, releff)
    releff = normeff/raneff
    return sample, no_tests, no_links_discovered, no_links_discovered_before, abeff, normeff, releff, evolution
    

sample_tbf, no_tests_tbf, no_links_discovered_tbf, no_links_discovered_before_tbf, abeff_tbf, normeff_tbf, releff_tbf, evolution_tbf = tbf(G, sample, links, n, m, no_tests_tbf, no_links_discovered_tbf, no_links_discovered_before_tbf, fileres_tbf, islink, pairstested_tbf, effbest_tbf, effworst_tbf, raneff_tbf, abeff_tbf, normeff_tbf, releff_tbf, evolution_tbf)

print('No. of tests in TBF Strategy ', no_tests_tbf)
print('No. of links detected in TBF Strategy ', no_links_discovered_tbf)
print("Absolute Efficiency for TBF Strategy is ", abeff_tbf)
print("Normalised Efficiency for TBF Strategy is ", normeff_tbf)
print("Relative Efficiency for TBF Strategy is ", releff_tbf)

m_tbf = len(sample_tbf.edges())
print('No. of edges ', m_tbf)

den = (2*m_tbf)/(n*(n-1))
print('Density ', den)

deg = degree(sample_tbf)
#print("Degree ", deg)

maxdeg = max(deg)
print("Maximum Degree ", maxdeg)

dd = degdis(sample_tbf, deg)

avg = avdeg(deg)
print("Average Degree ", avg)

tr = tri(sample_tbf, n, deg)
#print("No. of Triangles ", tr)

trip = triplet(sample_tbf, n)
#print("No. of Triplets ", trip)

cc = ccoef(tr, trip,  n)
print("Clustering Coefficient ", cc)

tratio = trratio(tr, trip, n)
print("Transitive Ratio ", tratio)


x, y = zip(*evolution_tbf)
plt.plot(x,y)
plt.title('TBF Strategy')
plt.xlabel('No. of tests')
plt.ylabel('No. of discovered links')
plt.grid()
plt.savefig("TBF_Strategy.png")
plt.show()

#Excercise 16
sample_mix = emptygraph(G)
sample_mix, no_tests_mix, no_links_discovered_mix, no_links_discovered_before_mix, raneff_mix, abeff_mix, normeff_mix, releff_mix, evolution_mix = randomstrategy(G, sample, n, m, no_tests_mix, no_links_discovered_mix, no_links_discovered_before_ran, fileres_mix, islink, pairstested_mix, effbest_mix, effworst_mix, raneff_mix, abeff_mix, normeff_mix, releff_mix, evolution_mix)
links = linksordered(sample_mix)
sample_mix, no_tests_mix, no_links_discovered_mix, no_links_discovered_before_mix, abeff_mix, normeff_mix, releff_mix, evolution_mix = tbf(G, sample_mix, links, n, m, no_tests_mix, no_links_discovered_mix, no_links_discovered_before_mix, fileres_mix, islink, pairstested_mix, effbest_mix, effworst_mix, raneff_mix, abeff_mix, normeff_mix, releff_mix, evolution_mix)
nodes = nodesordered(sample_mix)
sample_mix, no_tests_mix, no_links_discovered_mix, no_links_discovered_before_mix, abeff_mix, normeff_mix, releff_mix, evolution_mix = comp(G, sample_mix, nodes, n, m, no_tests_mix, no_links_discovered_mix, no_links_discovered_before_ran, fileres_mix, islink, pairstested_mix, effbest_mix, effworst_mix, raneff_mix, abeff_mix, normeff_mix, releff_mix, evolution_mix)

print('No. of tests in Mix Strategy ', no_tests_mix)
print('No. of links detected in Mix Strategy ', no_links_discovered_mix)
#print('Test No. for each edge found ', v)

print("Absolute Efficiency for Mix Strategy is ", abeff_mix)
print("Normalised Efficiency for Mix Strategy is ", normeff_mix)
print("Relative Efficiency for Mix Strategy is ", releff_mix)

m_mix = len(sample_mix.edges())
print('No. of edges ', m)

den = (2*m_mix)/(n*(n-1))
print('Density ', den)

deg = degree(sample_mix)
#print("Degree ", deg)

maxdeg = max(deg)
print("Maximum Degree ", maxdeg)

avg = avdeg(deg)
print("Average Degree ", avg)

tr = tri(sample_mix, n, deg)
#print("No. of Triangles ", tr)

trip = triplet(sample_mix, n)
#print("No. of Triplets ", trip)

cc = ccoef(tr, trip,  n)
print("Clustering Coefficient ", cc)

tratio = trratio(tr, trip, n)
print("Transitive Ratio  ", tratio)

x, y = zip(*evolution_mix)
plt.plot(x,y)
plt.title('Mixed Strategies')
plt.xlabel('No. of tests')
plt.ylabel('No. of discovered links')
plt.grid()
plt.savefig("Mixing_Strategies.png")
plt.show()