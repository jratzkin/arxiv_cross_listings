import time 
import arxiv
import networkx as nx
import matplotlib.pyplot as plt

#creation of classes

class arxiv_class:
    def __init__ (self, name, subclass_list, tot_cross_lists):
        self.name=name
        self.subclass_list = []
        self.tot_cross_lists = 0

class arxiv_subclass:
    kind='arxiv_class'
    def __init__(self, name, parent_class, cross_lists, tot_cross_lists):
        self.name=name
        self.parent_class = parent_class
        self.cross_lists = {}
        self.tot_cross_lists = 0
        
#initialize counts

tot_papers = 0
papers_w_crosslist = 0

#create a class for each arxiv subcategory
#is there a way to not hard-code this step????

math = arxiv_class('math',[], 0)
physics = arxiv_class('physics',[], 0)
cs = arxiv_class('cs',[], 0)
stat = arxiv_class('stat',[], 0)
q_bio = arxiv_class('q-bio',[], 0)
q_fin = arxiv_class('q-fin',[], 0)

math_ag = arxiv_subclass('math.AG',math,{}, 0) #algebraic geometry
math_at = arxiv_subclass('math.AT', math, {}, 0) #algebraic topology
math_ap = arxiv_subclass('math.AP', math, {}, 0) #analysis of PDE
math_ct = arxiv_subclass('math.CT', math, {}, 0) #category theory
math_ca = arxiv_subclass('math.CA', math, {}, 0) #classical analysis and ODE
math_co = arxiv_subclass('math.CO', math, {}, 0) #combinatorics
math_ac = arxiv_subclass('math.AC', math, {}, 0) #commutative algebra
math_cv = arxiv_subclass('math.CV', math, {}, 0) #complex variables
math_dg = arxiv_subclass('math.DG', math, {}, 0) #differential geometry
math_ds = arxiv_subclass('math.DS', math, {}, 0) #dynamical systems
math_fa = arxiv_subclass('math.FA', math, {}, 0) #functional analysis
math_gm = arxiv_subclass('math.GM', math, {}, 0) #general mathematics
math_gn = arxiv_subclass('math.GN', math, {}, 0) #general topology
math_gt = arxiv_subclass('math.GT', math, {}, 0) #geometric topology
math_gr = arxiv_subclass('math.GR', math, {}, 0) #group theory
math_ho = arxiv_subclass('math.HO', math, {}, 0) #history and overview
math_it = arxiv_subclass('math.IT', math, {}, 0) #information theory
math_kt = arxiv_subclass('math.KT', math, {}, 0) #K theory and homology
math_lo = arxiv_subclass('math.LO', math, {}, 0) #logic
math_mp = arxiv_subclass('math.MP', math, {}, 0) #mathematical physics
math_mg = arxiv_subclass('math.MG', math, {}, 0) #metric geometry
math_nt = arxiv_subclass('math.NT', math, {}, 0) #number theory
math_na = arxiv_subclass('math.NA', math, {}, 0) #numerical analysis
math_oa = arxiv_subclass('math.OA', math, {}, 0) #operator algebras
math_oc = arxiv_subclass('math.OC', math, {}, 0) #optimization and control
math_pr = arxiv_subclass('math.PR', math, {}, 0) #probability
math_qa = arxiv_subclass('math.QA', math, {}, 0) #quantum algebras
math_rt = arxiv_subclass('math.RT', math, {}, 0) #representation theory
math_ra = arxiv_subclass('math.RA', math, {}, 0) #rings and algebras
math_sp = arxiv_subclass('math.SP', math, {}, 0) #spectral theory
math_sg = arxiv_subclass('math.SG', math, {}, 0) #symplectic geometry

q_bio_bm = arxiv_subclass('q-bio.BM', q_bio, {}, 0) #biomolecules
q_bio_cb = arxiv_subclass('q-bio.CB', q_bio, {}, 0) #cell behavior
q_bio_gn = arxiv_subclass('q-bio.GN', q_bio, {}, 0) #genomics
q_bio_mn = arxiv_subclass('q-bio.MN', q_bio, {}, 0) #molecular networks
q_bio_nc = arxiv_subclass('q-bio.NC', q_bio, {}, 0) #neurons and cognition
q_bio_ot = arxiv_subclass('q-bio.OT', q_bio, {}, 0) #other quantitative biology
q_bio_pe = arxiv_subclass('q-bio.PE',q_bio, {}, 0) #populations and evolution
q_bio_qm = arxiv_subclass('q-bio.QM',q_bio, {}, 0) #quantitative methods
q_bio_sc = arxiv_subclass('q-bio.SC',q_bio, {}, 0) #subcellular processes
q_bio_to = arxiv_subclass('q-bio.TO',q_bio, {}, 0) #tissues and organs

q_fin_cp = arxiv_subclass('q-fin.CP', q_fin, {}, 0) #computational finance
q_fin_ec = arxiv_subclass('q-fin.EC', q_fin, {}, 0) #economics
q_fin_gn = arxiv_subclass('q-fin.GN', q_fin, {}, 0) #general finance
q_fin_mf = arxiv_subclass('q-fin.MF', q_fin, {}, 0) #mathematical finance
q_fin_pm = arxiv_subclass('q-fin.PM', q_fin, {}, 0) #portfolio management
q_fin_pr = arxiv_subclass('q-fin.PR', q_fin, {}, 0) #pricing of securities
q_fin_rm = arxiv_subclass('q-fin.RM', q_fin, {}, 0) #risk management
q_fin_st = arxiv_subclass('q-fin.ST', q_fin, {}, 0) #statistical finance
q_fin_tr = arxiv_subclass('q-fin.TR', q_fin, {}, 0) #trading and market microstructure

stat_ap = arxiv_subclass('stat.AP', stat, {}, 0) #applications
stat_co = arxiv_subclass('stat.CO', stat, {}, 0) #computation
stat_ml = arxiv_subclass('stat.ML', stat, {}, 0) #machine learning
stat_me = arxiv_subclass('stat.ME', stat, {}, 0) #methodology
stat_ot = arxiv_subclass('stat.OT', stat, {}, 0) #other statistics
stat_th = arxiv_subclass('stat.TH', stat, {}, 0) #statistics theory
math_st = stat_th

astro_ph_ga = arxiv_subclass('astro-ph.GA', physics, {}, 0) #astrophysics of galaxies
astro_ph_co = arxiv_subclass('astro-ph.CO', physics, {}, 0) #cosmology and nongalactic astrophysics
astro_ph_ep = arxiv_subclass('astro-ph.EP', physics, {}, 0) #earth and planetary astrophysics
astro_ph_he = arxiv_subclass('astro-ph.HE', physics, {}, 0) #high energy astrophysics
astro_ph_im = arxiv_subclass('astro-ph.IM', physics, {}, 0) #instrumentation and methods
astro_ph_sr = arxiv_subclass('astro-ph.SR', physics, {}, 0) #solar and stellar astrophysics

cond_mat_dis_nn = arxiv_subclass('cond-mat.dis.nn', physics, {}, 0) #disordered systems and neural networks
cond_mat_mtrl_sci = arxiv_subclass('cond-mat.mtrl-sci', physics, {}, 0) #material sciences
cond_mat_mes_hall = arxiv_subclass('cond-mat.mes-hall', physics, {}, 0) #mesoscale and nanoscale physics
cond_mat_other = arxiv_subclass('cond-math.other', physics, {}, 0) #other condensed matter physics
cond_mat_quant_gas = arxiv_subclass('cond-math.quant-gas', physics, {}, 0) #quantum gases
cond_mat_soft = arxiv_subclass('cond-mat.soft', physics, {}, 0) #soft condensed matter
cond_mat_stat_mech = arxiv_subclass('cond-mat.stat-mech', physics, {}, 0) #statistical mechanics
cond_mat_str_el = arxiv_subclass('cond-mat.str-el', physics, {}, 0) #strongly correlated electrons
cond_mat_supr_con = arxiv_subclass('cond-mat.supr-con', physics, {}, 0) #superconductivity

gr_qc = arxiv_subclass('gr-qc', physics, {}, 0) #general relativity and quantum cosmology
hep_ex = arxiv_subclass('hep-ex', physics, {}, 0) #high energy physics - experiments
hep_lat = arxiv_subclass('hep-lat', physics, {}, 0) #high energy physics - lattices
hep_th = arxiv_subclass('hep-th', physics, {}, 0) #high energy physics - theory
nlin_ao = arxiv_subclass('nlin.AO', physics, {}, 0) #adaptation and self-organizing systems
nlin_cg = arxiv_subclass('nlin.CG', physics, {}, 0) #celluar automata and lattice gases
nlin_cd = arxiv_subclass('nlin.CD', physics, {}, 0) #chaotic dynamics
nlin_si = arxiv_subclass('nlin.SI', physics, {}, 0) #exactly solvable and integrable systems
nlin_ps = arxiv_subclass('nlin.PS', physics, {}, 0) #pattern formation and solitons
nucl_ex = arxiv_subclass('nucl-ex', physics, {}, 0) #nuclear experiments
nucl_th = arxiv_subclass('nucl-th', physics, {}, 0) #nuclear theory
quant_ph = arxiv_subclass('quant-ph', physics, {}, 0) #quantum physics

physics_acc_ph = arxiv_subclass('physics.acc-ph', physics, {}, 0) #accelerator physics
physics_ao_ph = arxiv_subclass('physics.ao-ph', physics, {}, 0) # atmospheric and oceanic physics
physics_atom_ph = arxiv_subclass('physics.atom-ph', physics, {}, 0) #atomic physics
physics_atm_clus = arxiv_subclass('physics.atm-clus', physics, {}, 0) #atomic and molecular clusters
physics_bio_ph = arxiv_subclass('physics.bio-ph', physics,{}, 0) #biological physics
physics_chem_ph = arxiv_subclass('physics.chem-ph', physics, {}, 0) #chemical physics
physics_class_ph = arxiv_subclass('physics.class-ph', physics, {}, 0) #classical physics
physics_comp_ph = arxiv_subclass('physics.comp-ph', physics, {}, 0) #computational physics
physics_data_an = arxiv_subclass('physics.data-an', physics, {}, 0) #data analysis, statistics, and probability
physics_flu_dyn = arxiv_subclass('physics.flu-dyn',physics, {}, 0) #fluid dynamics
physics_gen_ph = arxiv_subclass('physics.gen-ph', physics, {}, 0) #general physics
physics_geo_ph = arxiv_subclass('physics.geo-ph', physics, {}, 0) #geophysics
physics_hist_ph = arxiv_subclass('physics.hist-ph', physics, {}, 0) #history and philosophy of physics
physics_ins_det = arxiv_subclass('physics.ins-det', physics, {}, 0) #instrumentation and detectors
physics_med_ph = arxiv_subclass('physics.med-ph', physics, {}, 0) #medical physics
physics_optics = arxiv_subclass('physics.optics', physics, {}, 0) #optics
physics_ed_ph = arxiv_subclass('physics.ed-ph', physics, {}, 0) #physics education
physics_soc_ph = arxiv_subclass('physics.soc-ph', physics, {}, 0) #physics and society
physics_plasm_ph = arxiv_subclass('physics.plasm-ph', physics, {}, 0) #plasma physics
physics_pop_ph = arxiv_subclass('physics.pop-ph', physics, {}, 0) #popular physics
physics_space_ph = arxiv_subclass('physics.space-ph', physics, {}, 0) #space physics

cs_ai = arxiv_subclass('cs.AI', cs, {}, 0) #artificial intelligence
cs_cl = arxiv_subclass('cs.CL', cs, {}, 0) #computation and language
cs_cc = arxiv_subclass('cs.CC',cs, {}, 0) #computational complexity
cs_ce = arxiv_subclass('cs.CE',cs, {}, 0) #computational engineering
cs_cg = arxiv_subclass('cs.CG', cs,{}, 0) #computational geometry
cs_gt = arxiv_subclass('cs.GT',cs, {}, 0) #computer science and game theory
cs_cv = arxiv_subclass('cs.CV',cs,{}, 0) #computer vision and pattern recognition
cs_cy = arxiv_subclass('cs.CY',cs, {}, 0) #computers and society
cs_cr = arxiv_subclass('cs.CR', cs,{}, 0) #cryptography and security
cs_ds = arxiv_subclass('cs.DS', cs, {}, 0) #data structures and algorithms
cs_db = arxiv_subclass('cs.DB', cs, {}, 0) #databases
cs_dl = arxiv_subclass('cs.DL', cs, {}, 0) #digital libraries
cs_dm = arxiv_subclass('cs.DM', cs,{}, 0) #discrete mathematics
cs_dc = arxiv_subclass('cs.DC', cs, {}, 0) #distributed, parallel, and cluster computing
cs_et = arxiv_subclass('cs.ET', cs, {}, 0) #emerging technologies
cs_fl = arxiv_subclass('cs.FL', cs, {}, 0) #formal languages and automata
cs_gl = arxiv_subclass('cs.GL', cs, {}, 0) #general literature
cs_gr = arxiv_subclass('cs.GR', cs, {}, 0) #graphics
cs_ar = arxiv_subclass('cs.AR', cs, {}, 0) #hardware architecture
cs_hc = arxiv_subclass('cs.HC', cs, {}, 0) #human-computer interaction
cs_ir = arxiv_subclass('cs.IR', cs, {}, 0) #information retrieval
cs_it = arxiv_subclass('cs.IT', cs, {}, 0) #information theory
cs_lg = arxiv_subclass('cs.LG', cs, {}, 0) #machine learning
cs_lo = arxiv_subclass('cs.LO', cs, {}, 0) #logic in computer science
cs_ms = arxiv_subclass('cs.MS', cs,{}, 0) #mathematical software
cs_ma = arxiv_subclass('cs.MA', cs,{}, 0) #multiagent systems
cs_mm = arxiv_subclass('cs.MM', cs, {}, 0) #multimedia
cs_ni = arxiv_subclass('cs.NI', cs, {}, 0) #networking and internet architecture
cs_ne = arxiv_subclass('cs.NE', cs, {}, 0) #neural and evolutionary computing
cs_na = arxiv_subclass('cs.NA', cs, {}, 0) #numerical analysis
cs_oa = arxiv_subclass('cs.OA', cs, {}, 0) #operating systems
cs_oh = arxiv_subclass('cs.OH', cs, {}, 0) #other computer science
cs_pf = arxiv_subclass('cs.PF', cs, {}, 0) #performace
cs_pl = arxiv_subclass('cs.PL', cs, {}, 0) #programming languages
cs_ro = arxiv_subclass('cs.RO', cs, {}, 0) #robotics
cs_si = arxiv_subclass('cs.SI', cs, {}, 0) #social and information networks
cs_se = arxiv_subclass('cs.SE', cs, {}, 0) #software engineering
cs_sd = arxiv_subclass('cs.SD', cs, {}, 0) #sound
cs_sc = arxiv_subclass('cs.SC', cs, {}, 0) #symbolic computation
cs_sy = arxiv_subclass('cs.SC', cs, {}, 0) #systems and control

#initialize various lists

list_of_classes = [math, physics, stat, q_bio, q_fin, cs]
list_of_math_subclasses = [math_ag, math_at, math_ap, math_ct, math_ca, math_co, math_ac, math_cv, math_dg, math_ds, math_fa, math_gm, math_gn, math_gt, math_gr, math_ho, math_it, math_kt, math_lo, math_mp, math_mg, math_nt, math_na, math_oa, math_oc, math_pr, math_qa, math_rt, math_ra, math_sp, math_st, math_sg]
list_of_qbio_subclasses = [q_bio_bm, q_bio_cb, q_bio_gn, q_bio_mn, q_bio_nc, q_bio_ot, q_bio_pe, q_bio_qm, q_bio_sc, q_bio_to]
list_of_qfin_subclasses = [q_fin_cp, q_fin_ec, q_fin_gn, q_fin_mf, q_fin_pm, q_fin_pr, q_fin_rm, q_fin_st, q_fin_tr]
list_of_stat_subclasses = [stat_ap, stat_co, stat_ml, stat_me, stat_ot, stat_th]
list_of_cs_subclasses = [cs_ai, cs_cl, cs_cc, cs_ce, cs_cg, cs_gt, cs_cv, cs_cy, cs_cr, cs_ds, cs_db, cs_dl, cs_dm, cs_dc, cs_et, cs_fl, cs_gl, cs_gr, cs_ar, cs_hc, cs_ir, cs_it, cs_lg, cs_lo, cs_ms, cs_ma, cs_mm, cs_ni, cs_ne, cs_na, cs_oa, cs_oh, cs_pf, cs_pl, cs_ro, cs_si, cs_se, cs_sd, cs_sc, cs_sy]
list_of_physics_subclasses = [astro_ph_ga, astro_ph_co, astro_ph_ep, astro_ph_he, astro_ph_im, astro_ph_sr, cond_mat_dis_nn, cond_mat_mtrl_sci, cond_mat_mes_hall, cond_mat_other, cond_mat_quant_gas, cond_mat_soft, cond_mat_stat_mech, cond_mat_str_el, cond_mat_supr_con, gr_qc, hep_lat, hep_ex, hep_th, nlin_ao, nlin_cg, nlin_cd, nlin_si, nlin_ps, nucl_ex, nucl_th, quant_ph, physics_acc_ph, physics_ao_ph, physics_atom_ph, physics_atm_clus, physics_bio_ph, physics_chem_ph, physics_comp_ph, physics_data_an, physics_flu_dyn, physics_gen_ph, physics_geo_ph, physics_hist_ph, physics_ins_det, physics_med_ph, physics_optics, physics_ed_ph, physics_soc_ph, physics_plasm_ph, physics_pop_ph, physics_space_ph]

list_of_math_names = []
for cat in list_of_math_subclasses:
    list_of_math_names.append(cat.name)
list_of_qbio_names = []
for cat in list_of_qbio_subclasses:
    list_of_qbio_names.append(cat.name)
list_of_qfin_names = []
for cat in list_of_qfin_subclasses:
    list_of_qfin_names.append(cat.name)
list_of_stat_names = []
for cat in list_of_stat_subclasses:
    list_of_stat_names.append(cat.name)
list_of_cs_names = []
for cat in list_of_cs_subclasses:
    list_of_cs_names.append(cat.name)
list_of_physics_names = []
for cat in list_of_physics_subclasses:
    list_of_physics_names.append(cat.name)

list_of_subclasses = list_of_math_subclasses + list_of_qbio_subclasses + list_of_qfin_subclasses + list_of_stat_subclasses + list_of_cs_subclasses + list_of_physics_subclasses
list_of_subclass_names = []
for subclass in list_of_subclasses:
    list_of_subclass_names.append(subclass.name)

math.subclass_list = list_of_math_subclasses
stat.subclass_list = list_of_stat_subclasses
cs.subclass_list = list_of_cs_subclasses
q_bio.subclass_list = list_of_qbio_subclasses
q_fin.subclass_list = list_of_qfin_subclasses
physics.subclass_list = list_of_physics_subclasses

#initialize cross lists

for subclass in list_of_subclasses:
    for possible_cross in list_of_subclasses:
        if not(subclass.name == possible_cross.name):
            subclass.cross_lists[possible_cross.name] = 0

#generate dictionaries of number of cross-lists between various subclasses
#time period is January 2017
#NB: this is a directed graph

date_url='+AND+submittedDate:[2017010100+TO+201701311159]'

for subclass in list_of_subclasses:
    time.sleep(5)
    print ('counting cross-lists in category %s now' % subclass.name)
    results = arxiv.query(subclass.name+date_url, max_results=1000)
    tot_papers = tot_papers + len(results)
    for result in results:
        if len(result.tags) > 1: 
            terms = []
            for i in range (1, len(result.tags)):
                terms.append(str(result.tags[i].term))
            crossed = False
            for term in terms: 
                if term in list(subclass.cross_lists.keys()):
                    crossed = True
                    subclass.cross_lists[term] += 1
                if crossed: 
                    papers_w_crosslist += 1
    print ('the number of papers with a cross-list is now %d' % papers_w_crosslist)            

#sum number of cross listings in various categories

for subclass in list_of_subclasses:
    for cat in list_of_subclasses:
        if cat.name in list(subclass.cross_lists.keys()):
            subclass.tot_cross_lists += subclass.cross_lists[cat.name]

for subclass in list_of_math_subclasses:
    math.tot_cross_lists += subclass.tot_cross_lists
for subclass in list_of_cs_subclasses:
    cs.tot_cross_lists += subclass.tot_cross_lists
for subclass in list_of_physics_subclasses: 
    physics.tot_cross_lists += subclass.tot_cross_lists
for subclass in list_of_qbio_subclasses:
    q_bio.tot_cross_lists += subclass.tot_cross_lists
for subclass in list_of_qfin_subclasses:
    q_fin.tot_cross_lists += subclass.tot_cross_lists

#print out list of cross-listing results 
   
print ('-'*20)
      
for category in list_of_subclasses:
    print('these are the cross-list numbers for %s' % category.name)
    for cat in list(category.cross_lists.keys()):
        if category.cross_lists[cat]>0:
            print (cat, category.cross_lists[cat])
    print ('-'*20)
          
print('%d articles have a cross-listing out of %d total articles' % (papers_w_crosslist, tot_papers))


#create directed graph representing all the cross-lists found
#edges of graph are weighted by the number of cross-lists

cross_list_graph = nx.DiGraph()
cross_list_graph.add_nodes_from(list_of_subclass_names) 
for category in list_of_subclasses:
    for cat in list (category.cross_lists.keys()):
        if category.cross_lists[cat] >= 10:
            cross_list_graph.add_edge(category, cat, weight=category.cross_lists[cat])

math_cross_list_graph = nx.DiGraph()
math_cross_list_graph.add_nodes_from(list_of_math_names)
for category in list_of_math_subclasses:
    for cat in list(category.cross_lists.keys()):
        if (cat in list_of_math_names) and (category.cross_lists[cat] >= 10):
            math_cross_list_graph.add_edge(category, cat, weight=category.cross_lists[cat])

physics_cross_list_graph = nx.DiGraph()
physics_cross_list_graph.add_nodes_from (list_of_physics_names)
for category in list_of_physics_subclasses:
    for cat in list(category.cross_lists.keys()):
        if (cat in list_of_physics_names) and (category.cross_lists[cat] >= 10):
            physics_cross_list_graph.add_edge(category, cat, weight=category.cross_lists[cat])

stat_cross_list_graph = nx.DiGraph()
stat_cross_list_graph.add_nodes_from (list_of_stat_names)
for category in list_of_stat_subclasses:
    for cat in list(category.cross_lists.keys()):
        if (cat in list_of_stat_names) and (category.cross_lists[cat] >= 10):
            stat_cross_list_graph.add_edge(category, cat, weight=category.cross_lists[cat])

cs_cross_list_graph = nx.DiGraph()
cs_cross_list_graph.add_nodes_from (list_of_cs_names)
for category in list_of_cs_subclasses:
    for cat in list(category.cross_lists.keys()):
        if (cat in list_of_cs_names) and (category.cross_lists[cat] >= 10):
            cs_cross_list_graph.add_edge(category, cat, weight=category.cross_lists[cat])

qbio_cross_list_graph = nx.DiGraph()
qbio_cross_list_graph.add_nodes_from (list_of_qbio_names)
for category in list_of_qbio_subclasses:
    for cat in list(category.cross_lists.keys()):
        if (cat in list_of_qbio_names) and (category.cross_lists[cat] >= 10):
            qbio_cross_list_graph.add_edge(category, cat, weight=category.cross_lists[cat])

qfin_cross_list_graph = nx.DiGraph()
qfin_cross_list_graph.add_nodes_from (list_of_physics_names)
for category in list_of_qfin_subclasses:
    for cat in list(category.cross_lists.keys()):
        if (cat in list_of_qfin_names) and (category.cross_lists[cat] >= 10):
            qfin_cross_list_graph.add_edge(category, cat, weight=category.cross_lists[cat])


#display the graph using matplotlib

fig=plt.figure(figsize=(9,9))

pos = nx.circular_layout(cross_list_graph, scale=5.0)        

nx.draw_networkx_nodes(cross_list_graph, pos, nodelist=list_of_math_names, color='r', node_size=25)
nx.draw_networkx_nodes(cross_list_graph, pos, nodelist=list_of_physics_names, color='b', node_size=25)
nx.draw_networkx_nodes(cross_list_graph, pos, nodelist=list_of_cs_names, color='g', node_size=25)
nx.draw_networkx_nodes(cross_list_graph, pos, nodelist=list_of_stat_names, color='y', node_size=25)
nx.draw_networkx_nodes(cross_list_graph, pos, nodelist=list_of_qbio_names, color='k', node_size=25)
nx.draw_networkx_nodes(cross_list_graph, pos, nodelist=list_of_qfin_names, color='c', node_size=25)

nx.draw_networkx_edges(cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv cross-listings for January 2017')
plt.savefig("graph_of_arxiv_crosslists.png")
plt.show()

pos = nx.circular_layout(math_cross_list_graph, scale=5.0)
nx.draw_networkx_nodes(math_cross_list_graph, pos, color='r', node_size=25)
nx.draw_networkx_edges(math_cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv mathematics cross-listsings')
plt.savefig("graph_of_math_arxiv_crosslists.png")
plt.show()

pos = nx.circular_layout(physics_cross_list_graph, scale=5.0)
nx.draw_networkx_nodes(physics_cross_list_graph, pos, color='b', node_size=25)
nx.draw_networkx_edges(physics_cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv physics cross-listsings')
plt.savefig("graph_of_physics_arxiv_crosslists.png")
plt.show()

pos = nx.circular_layout(cs_cross_list_graph, scale=5.0)
nx.draw_networkx_nodes(cs_cross_list_graph, pos, color='g', node_size=25)
nx.draw_networkx_edges(cs_cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv cs cross-listsings')
plt.savefig("graph_of_cs_arxiv_crosslists.png")
plt.show()

pos = nx.circular_layout(stat_cross_list_graph, scale=5.0)
nx.draw_networkx_nodes(stat_cross_list_graph, pos, color='y', node_size=25)
nx.draw_networkx_edges(stat_cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv statistics cross-listsings')
plt.savefig("graph_of_stat_arxiv_crosslists.png")
plt.show()

pos = nx.circular_layout(qbio_cross_list_graph, scale=5.0)
nx.draw_networkx_nodes(qbio_cross_list_graph, pos, color='k', node_size=25)
nx.draw_networkx_edges(qbio_cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv quantitative biology cross-listsings')
plt.savefig("graph_of_qbio_arxiv_crosslists.png")
plt.show()

pos = nx.circular_layout(qfin_cross_list_graph, scale=5.0)
nx.draw_networkx_nodes(qfin_cross_list_graph, pos, color='c', node_size=25)
nx.draw_networkx_edges(qfin_cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv quantitative finance cross-listsings')
plt.savefig("graph_of_qfin_arxiv_crosslists.png")
plt.show()

pos = nx.circular_layout(math_cross_list_graph, scale=5.0)
nx.draw_networkx_nodes(math_cross_list_graph, pos, color='r', node_size=25)
nx.draw_networkx_edges(math_cross_list_graph, pos, width=0.25)

plt.axis('off')
plt.title('plot of arxiv mathematics cross-listsings')
plt.savefig("graph_of_math_arxiv_crosslists.png")
plt.show()

fig = plt.figure(figsize=(6,6))
X = [0,1,2,3,4,5]
Y = [math.tot_cross_lists, physics.tot_cross_lists, cs.tot_cross_lists, stat.tot_cross_lists, q_bio.tot_cross_lists, q_fin.tot_cross_lists]
names = ['math', 'physics', 'cs', 'stat', 'q-bio', 'q-fin']

plt.bar(X, Y)
plt.xticks(X, names) 
plt.ylabel('total number of cross lists')
plt.title('bar plot of total cross lists on arXiv, Jan. 2017')
plt.savefig("bar_plot_of_tot_cross_lists.png")
plt.show()
