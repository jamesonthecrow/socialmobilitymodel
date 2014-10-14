import numpy as np 
import networkx as nx 
import matplotlib.pyplot as plt 
import matplotlib as mpl
#mpl.rc('text', usetex=True)

DIR = '/Users/jltoole/Documents/Projects/geo_social/'
class Model():
	def __init__(self, u, l, c, s, a):
		self.NUSERS = u
		self.NLOCS = l
		self.NCONTACTS = c
		self.NSTEPS = s
		self.G = nx.Graph()

		# CONSTANTS
		self.rho = 0.6
		self.gamma = 0.21
		self.alpha = a
		self.beta = 0.8
		self.tau = 17.0; # exponential cuttoff on time between calls
		self.xmin = 1.0;

	def reset(self):
		print 'Initializing Graph...'
		
		self.G = nx.Graph()
		nodes = np.arange(self.NUSERS)
		degs = np.random.lognormal(np.log(self.NCONTACTS),0.3, self.NUSERS)
		# add nodes

		# assign them random locations to start
		a = np.arange(self.NLOCS,dtype=float)
		p = a/np.sum(a)
		p = np.cumsum(p)
		for i in nodes:
			self.G.add_node(i, lvec=np.zeros(self.NLOCS,dtype=int), a=np.random.exponential(scale=self.alpha), S=3)
			for l in xrange(3):
				r = np.random.rand(3)
				l = np.digitize(r, p)
				#self.G.node[i]['lvec'][l] = 1
				self.G.node[i]['lvec'][np.random.randint(self.NLOCS)] = 1

		# connect the network
		for i in nodes:
			stubs_left = degs[i] - self.G.degree(i)
			if stubs_left > 0:
				nbrs = []
				while len(nbrs) < stubs_left:
					tries = 0
					j = nodes[np.random.randint(self.NUSERS)]
					while (((degs[j] - self.G.degree(j) <= 0) or (i==j)) and (tries < 1000)) :
						j = nodes[np.random.randint(self.NUSERS)]
						tries += 1
					nbrs.append(j)
				edges = [ (i,j, {'sim':None}) for j in nbrs ]
				self.G.add_edges_from(edges)
			if i%(self.NUSERS/10) == 0:
				print i, self.NUSERS
		
		'''
		self.G =  nx.newman_watts_strogatz_graph(self.NUSERS,self.NCONTACTS,0.9)
		for i in self.G.nodes_iter():
			self.G.add_node(i, lvec=np.zeros(self.NLOCS,dtype=int), S=2)
			for l in xrange(2):
				self.G.node[i]['lvec'][np.random.randint(self.NLOCS)] = 1
		for e in self.G.edges_iter():
			self.G[e[0]][e[1]]['sim'] = None
		'''


	def get_return_location(self, u):
		''' choose a location to return to based on a preferential attachement model
		'''
		lvec = self.G.node[u]['lvec']
		p = np.cumsum(lvec)/float(np.sum(lvec))
		r = np.random.rand()
		return np.digitize( [r], p )[0]

	def get_randomuser_location(self, u):
		I = np.where( self.G.node[u]['lvec'] > 0)[0]
		return I[np.random.randint(len(I))]

	def get_friend_location(self, u):
		''' pick a friend and choose one of their locations both based on pref. attachment
		'''
		p = np.arange(1,self.G.degree(u))
		p = np.cumsum(p)/float(np.sum(p))
		r = np.random.rand()
		f = np.digitize( [r], p )[0]
		return self.get_return_location(self.G.neighbors(u)[f])
		#return self.get_randomuser_location(self.G.neighbors(u)[f])

	def get_citywide_location(self,p):
		r = np.random.rand()
		return np.digitize( [r], p )[0]

	def get_random_location(self):
		return np.random.randint(self.NLOCS)

	def get_citywide_visits(self):
		lvecs = np.array(nx.get_node_attributes(self.G,'lvec').values())
		p = np.sum(lvecs,0).astype(float)/np.sum(lvecs)
		return p

	def run(self):
		print 'Running Model...'
		for t in xrange(self.NSTEPS):
			nextlocs = np.zeros(self.NUSERS)
			p = np.cumsum(self.get_citywide_visits())
			for u in xrange(self.NUSERS):
				r = np.random.rand()
				tries = 0
				l = None
				if r > self.rho*self.G.node[u]['S']**(-self.gamma):
					# preferential return
					r = np.random.rand()
					if r > self.G.node[u]['a']:
						l = self.get_return_location(u)
					else:
						l = self.get_friend_location(u)
						while (self.G.node[u]['lvec'][l] == 0) and (tries < 100):
							l = self.get_friend_location(u)
							tries += 1
					
				else:
					# explore
					if r > self.G.node[u]['a']:
						#l = self.get_citywide_location(p)
						l = self.get_random_location()
						while (self.G.node[u]['lvec'][l] > 0) and (tries < 100):
							#l = self.get_citywide_location(p)
							l = self.get_random_location()
					else:
						l = self.get_friend_location(u)
						while (self.G.node[u]['lvec'][l] > 0) and (tries < 100):
							#l = self.get_citywide_location(p)
							l = self.get_friend_location(u)
							tries += 1
					self.G.node[u]['S'] += 1
				nextlocs[u] = l
					
			# update users
			for u in xrange(self.NUSERS):
				self.G.node[u]['lvec'][nextlocs[u]] += 1
			print "Step: ", t

	def calculate_similarity(self):
		print 'Calculating Similarity...'
		for i in self.G.nodes_iter():
			l1 = self.G.node[i]['lvec']
			for j in self.G.neighbors(i):
				l2 = self.G.node[j]['lvec']
				self.G.edge[i][j]['sim'] = cosine_similarity( l1,l2 )

def cosine_similarity( u, v ):
		u = u.astype(float)
		v = v.astype(float)
		return np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) 

def figure1(nets, leg=['Data','Randomized'], outfile=None):

	plt.close()
	f1 = plt.figure(1)
	for G in nets:
		data = np.array(nx.get_edge_attributes(G,'sim').values())
		rand_data = []
		for u in G.nodes_iter():
			for k in xrange(2):
				j = np.random.randint(G.number_of_nodes())
				rand_data.append( cosine_similarity(G.node[u]['lvec'], G.node[j]['lvec'] ))
			#rand_data.append(np.mean(temp))

		data = np.array(data)
		rand_data = np.array(rand_data)

		# sim
		xbins = np.logspace(-2,0,num=30)
		#xbins = np.linspace(0,1,num=30)
		x = xbins[1:] - (xbins[1:]-xbins[:-1])/2.
		y, edges = np.histogram(data,xbins, density=True)
		plt.subplot(2,1,1)
		plt.loglog( x, y, '.-',color=np.random.rand(3))
		#plt.legend(['Data','Randomized'], loc='best', fontsize=6)
		
		plt.subplot(2,1,2)
		data = np.array(nx.get_node_attributes(G,'lvec').values()).astype('float')
		data = np.sort(data,1)[:,::-1]
		data = data/np.tile(np.sum(data,1),(data.shape[1],1)).T
		y = np.mean(data,0)
		x = np.arange(len(y))+1
		plt.loglog(x, y, 'k.-')

	plt.subplot(2,1,1)
	y, edges = np.histogram(rand_data,xbins, density=True)
	x = edges[1:] - (edges[1:]-edges[:-1])/2.
	plt.loglog( x, y, 'k.-')

	plt.xlim([0,1])
	plt.ylim([10**-2, 10**3])
	plt.xlabel('$\cos\phi$')
	plt.ylabel('$P(\cos\phi)$')
	plt.title('Cosine Similarity')
	leg.append('Randomized')
	plt.legend(leg, loc='best', fontsize=8)


	plt.subplot(2,1,2)
	plt.xlabel('$k$')
	plt.ylabel('$f_k$')

	f1.set_size_inches(4,6)
	f1.tight_layout()
	if outfile != None:
		plt.savefig( DIR+'figures/'+outfile+'_model.png',fmt='png')
	plt.close()

'''
import chaomodel_graph as cm
m.figure1(m.G, outfile='figure1')
m = cm.Model(10000,250,10,55,0.1)

reload(cm)
alpha = [0.0,0.25,0.5,0.75,1.0]
nets = []
for a in alpha:
	m = cm.Model(10000, 250, 20, 50, a); 
	m.reset(); 
	m.run(); 
	m.calculate_similarity();
	nets.append(m.G)

cm.figure1(nets, leg=[str(a) for a in alpha], outfile='chao_alpha')

reload(cm)
m = cm.Model(10000,250,20,100,0.10); m.reset(); m.run(); m.calculate_similarity()
cm.figure1([m.G], leg=['Data'], outfile='chao')
'''

