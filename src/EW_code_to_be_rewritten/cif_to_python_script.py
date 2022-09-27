import numpy as np
import re

atom_num=108
cage='cc1'

coords=[]
lat_pat=[]
coords_1=np.zeros((atom_num,3))
coord_type=[]
n=8
temp=101
filename='out_%d_8.cif'%(temp)
with open(filename, 'r') as f:
	lines=f.readlines() 
	new=[i.split() for i in lines]
for i in range(20*n*n):
	coords.append(new[i+21][1:4])
	coord_type.append(new[i+21][0])
trunc_coord_type=coord_type	
for i in range(20*n*n):
	trunc_coord_type[i]=coord_type[i][:2]
lat_pat.append(new[8][1])
lat_pat.append(new[9][1])
lat_pat.append(new[10][1])
print(coords)

#filename='out_1_8.xyz'
#with open(filename, 'r') as f:
#	lines=f.readlines() 
#	new=[i.split() for i in lines]
#print(len(new))	
#for i in range(len(new)-2):
#	coords.append(new[i+2][1:4])
#trunc_coord_type=coord_type	
	
	#print(coords[i])
#	
#print(trunc_coord_type)
coords_ar=np.asarray(coords, dtype=np.float64, order='C' )
#print(coords_ar)
newarr = coords_ar.reshape(20*n*n, 3)
#print(coords_ar.shape)
for i in range(20*n*n):
	coords_ar[i][0]=coords_ar[i][0]*float(lat_pat[0])
	coords_ar[i][1]=coords_ar[i][1]*float(lat_pat[1])
	coords_ar[i][2]=coords_ar[i][2]*float(lat_pat[2])
print(coords_ar.shape)

python_script=open("stk_from_cg_%d.py"%(temp),"w+")

python_script.write("""import stk
import stko
import numpy as np
from functools import partial
from operator import getitem
from collections import Counter

n=8
temp=%d

class new_non_linear_vertex(stk.molecular.topology_graphs.cof.vertices.NonLinearVertex):
	def place_building_block(self, building_block, edges):
		print(self.get_id(), building_block, edges)
		assert (
			building_block.get_num_functional_groups() > 2
		), (
			f'{building_block} needs to have more than 2 functional '
			'groups but has '
			f'{building_block.get_num_functional_groups()}.'
		)
		edges = sorted(edges, key=lambda edge: edge.get_parent_id())	
		building_block = building_block.with_centroid(
			position=self._position,
			atom_ids=building_block.get_placer_ids(),
		)
		core_centroid = building_block.get_centroid(
			atom_ids=building_block.get_core_atom_ids(),
		)
		normal = building_block.get_plane_normal(
			atom_ids=building_block.get_placer_ids(),
		)
		normal = get_acute_vector(
			reference=core_centroid - self._position,
			vector=normal,
		)
		building_block = building_block.with_rotation_between_vectors(
			start=normal,
			target=[0, 0, 1],
			origin=self._position,
		)
		fg, = building_block.get_functional_groups(0)
		fg_centroid = building_block.get_centroid(fg.get_placer_ids())
		edge_position = edges[self._aligner_edge].get_position()
		return building_block.with_rotation_to_minimize_angle(
			start=fg_centroid - self._position,
			target=edge_position - self._position,
			axis=np.array([0, 0, 1], dtype=np.float64),
			origin=self._position,
		).get_position_matrix()
		
		

class CUSTOM_PERIODIC_TOPOLOGY(stk.molecular.topology_graphs.cof.cof.Cof):


#	def __init_subclass__(cls, **kwargs):
		




	def __init__(
	    self,
	    building_blocks,
	    lattice_size,
	    periodic=True,
	    vertex_alignments=None,
	    reaction_factory=stk.GenericReactionFactory(),
	    num_processes=1,
	    optimizer=stk.NullOptimizer(),
	):
	
	
		self._vertex_alignments = (
		    dict(vertex_alignments)
		    if vertex_alignments is not None
		    else {}
		)
		self._lattice_size = lattice_size
		self._periodic = periodic
		
		lattice = self._get_lattice(self._vertex_alignments)
		edges = self._get_edges(lattice)
		vertices = self._get_vertices(lattice)
		
		#for v in vertices:
		#	print(v.get_id(),v,v.__class__.__name__)
		#	if (v.get_id()>4):
		#		sys.exit()
		vertex_degrees = Counter(
			
			vertex_id
			
			for edge in self._edge_prototypes
			
			for vertex_id in edge.get_vertex_ids()
			
		)
		print(vertex_degrees)
		self._allowed_degrees = set([2,3])
#		self._allowed_degrees = set(vertex_degrees.values())
		
		if isinstance(building_blocks, dict):
			for building_block in building_blocks:
				print('hi')
				assert (
		            building_block.get_num_functional_groups()
		            in self._allowed_degrees
		        ), (
		            'The number of functional groups in '
		            f'{building_block} needs to be one of '
		            f'{tuple(self._allowed_degrees)}, but is '
		            'currently '
		            f'{building_block.get_num_functional_groups()}.'
		        )
				get_vertex = partial(getitem, vertices)
				building_block_vertices = {
				    building_block: tuple(map(
				        get_vertex,
				        # Account for the fact that a building block can
				        # be mapped to a single int.
				        (ids, ) if isinstance(ids, int) else ids
				    ))
				    for building_block, ids in building_blocks.items()
				}
		else:
		    building_block_vertices = (
		        self._get_building_block_vertices(
		            building_blocks=building_blocks,
		            vertices=vertices,
		            edges=edges,
		        )
		    )
		
		building_block_vertices = self._with_unaligning_vertices(
		    building_block_vertices=building_block_vertices,
		)
		
		self._check_building_block_vertices(
		    num_vertices=(
		        np.product(lattice_size)*len(self._vertex_prototypes)
		    ),
		    building_block_vertices=building_block_vertices,
		)
		super(stk.molecular.topology_graphs.cof.cof.Cof, self).__init__(
		    building_block_vertices=building_block_vertices,
		    edges=edges,
		    reaction_factory=reaction_factory,
		    construction_stages=(),
		    num_processes=num_processes,
		    optimizer=optimizer,
		    edge_groups=self._get_edge_groups(edges),
		)

	def _get_vertices(self, lattice):

		xdim, ydim, zdim = self._lattice_size
		num_vertices = xdim*ydim*zdim*len(self._vertex_prototypes)
		vertices = [None for i in range(num_vertices)]
		for vertex in self._vertex_prototypes:
			vertices[vertex.get_id()] = vertex

		return tuple(vertices)	

	
	def _get_scale(self, building_block_vertices):
	    return 1
	
	def get_periodic_info(self):
		warnings.warn(
		    'You called get_periodic_info() on a topology graph '
		    'instance. This method will be removed in any version '
		    'of stk released on, or after, 21/10/21. Please call '
		    'the construct() method instead. This will return a '
		    'PeriodicConstructionResult which provides the new '
		    'get_periodic_info() method.'
		)
		
		lattice_constants = self._get_lattice_constants()
		
		return PeriodicInfo(
		    vector_1=(
		        lattice_constants[0]*self._lattice_size[0]*self._scale
		    ),
		    vector_2=(
		        lattice_constants[1]*self._lattice_size[1]*self._scale
		    ),
		    vector_3=(
		        lattice_constants[2]*self._lattice_size[2]*self._scale
		    ),
		)

	_lattice_constants = _a, _b, _c = (
		np.array([%9.5f, 0., 0.]),
		np.array([%9.5f, %9.5f, 0.]),	
		np.array([0., 0., 5])
	)\n\n
"""%(temp,float(lat_pat[0]),-float(lat_pat[0])*0.5,-float(lat_pat[1])*np.sqrt(3)/2))


input3=open("output_%d.txt"%(temp), "w+")

input3.write('	_lattice_constants = _a, _b, _c = (\n')
input3.write('		np.array([%9.5f, 0., 0.]),\n'%(float(lat_pat[0])))
input3.write('		np.array([%9.5f, %9.5f, 0.]),\n'%(-float(lat_pat[0])*0.5,-float(lat_pat[1])*np.sqrt(3)/2))
input3.write('		np.array([0., 0., %s])\n	)\n\n'%(lat_pat[2]))

python_script.write('	_vertex_prototypes = ( \n')

iter1=-1
nonlin=[]
for i in range(20*n*n):
	if (coord_type[i]=='Ti'):
		iter1=iter1+1
		python_script.write('		stk.molecular.topology_graphs.cof.vertices.NonLinearVertex(%d,(%s,%s,%s)),\n'%(iter1,coords_ar[i][0]-coords_ar[i][1]*0.5,-coords_ar[i][1]*np.sqrt(3)/2,coords_ar[i][2]))
		nonlin.append(i)

python_script.write(')\n\n')


python_script.write('	_vertex_prototypes = ( \n')
python_script.write('		*_vertex_prototypes, \n')
lin=[]

for i in range(20*n*n):
	if (coord_type[i]!='Ti'):
		iter1=iter1+1
		python_script.write('		stk.molecular.topology_graphs.cof.vertices.LinearVertex(\n			id=%d,\n			position=(%s,%s,%s),	\n		),\n'%(iter1,coords_ar[i][0]-coords_ar[i][1]*0.5,-coords_ar[i][1]*np.sqrt(3)/2,coords_ar[i][2]))
		lin.append(i)
python_script.write('	)\n')

#I think this is essentially saying which guys are between the other guys
#Think you should print out the nn shell and read it in and its like if it looped (which i think you can check from having to take away 1... maybe)

#for i in range(20*n*n):
#	print(coords_ar[i])

iter=-1
vec1=np.zeros((2))
vec2=np.zeros((2))
vec3=np.zeros((2))
vec4=np.zeros((2))
vec5=np.zeros((2))
vec6=np.zeros((2))
vec7=np.zeros((2))
vec8=np.zeros((2))
vec9=np.zeros((2))
period=np.zeros((2))
python_script.write('	_edge_prototypes = ( \n')
print(lat_pat[0])
for i in range(len(nonlin)):
	for j in range(len(lin)):
		#print(coords_ar[nonlin[i]],i)
		#print(coords_ar[lin[j]],j+len(nonlin))
		
		
		#x=coords_ar[i][0]-coords_ar[i][1]*0.5
		#y=coords_ar[i][1]*np.sqrt(3)/2
		
		
		#vec1[0]=(coords_ar[nonlin[i]][0]-coords_ar[lin[j]][0])
		#vec1[1]=-coords_ar[nonlin[i]][0]/2+(coords_ar[lin[j]][0])/2-(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2
		#
		#vec2[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0]))-(coords_ar[lin[j]][0])	#+x
		#vec2[1]=-(coords_ar[nonlin[i]][0]+float(lat_pat[0]))/2+(coords_ar[lin[j]][0])/2-(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2										#+x
		#
		#vec3[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0]))-(coords_ar[lin[j]][0])	#+x+y
		#vec3[1]=-(coords_ar[nonlin[i]][0]+float(lat_pat[0]))/2+(coords_ar[lin[j]][0])/2-(coords_ar[nonlin[i]][1]-float(lat_pat[1])-coords_ar[lin[j]][1])*np.sqrt(3)/2 	#+x+y
		#
		#vec4[0]=(coords_ar[nonlin[i]][0])-(coords_ar[lin[j]][0])	#+y
		#vec4[1]=-coords_ar[nonlin[i]][0]/2+(coords_ar[lin[j]][0])/2-(coords_ar[nonlin[i]][1]-float(lat_pat[1])-coords_ar[lin[j]][1])*np.sqrt(3)/2		#+y
		#
		#vec5[0]=(coords_ar[nonlin[i]][0])-(coords_ar[lin[j]][0]+float(lat_pat[0]))	#-x
		#vec5[1]=-coords_ar[nonlin[i]][0]/2+(coords_ar[lin[j]][0]+float(lat_pat[0]))/2-(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2
		#
		#vec6[0]=(coords_ar[nonlin[i]][0])-(coords_ar[lin[j]][0]+float(lat_pat[0]))	#-x-y
		#vec6[1]=-coords_ar[nonlin[i]][0]/2+(coords_ar[lin[j]][0]+float(lat_pat[0]))/2-(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2
		#
		#vec7[0]=(coords_ar[nonlin[i]][0]-coords_ar[lin[j]][0])
		#vec7[1]=-coords_ar[nonlin[i]][0]/2+(coords_ar[lin[j]][0])/2-(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2					
		#
		#vec8[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0]))-(coords_ar[lin[j]][0])	#x-y
		#vec8[1]=-(coords_ar[nonlin[i]][0]+float(lat_pat[0]))/2+(coords_ar[lin[j]][0])/2-(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2	
		#
		#vec9[0]=(coords_ar[nonlin[i]][0])-(coords_ar[lin[j]][0]+float(lat_pat[0]))	#-x+y
		#vec9[1]=-coords_ar[nonlin[i]][0]/2+(coords_ar[lin[j]][0]+float(lat_pat[0]))/2-(coords_ar[nonlin[i]][1]-float(lat_pat[1])-(coords_ar[lin[j]][1]))*np.sqrt(3)/2	

	
	
		vec1[0]=(coords_ar[nonlin[i]][0]-coords_ar[nonlin[i]][1]*0.5)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]*0.5)
		vec1[1]=-(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2
		
		vec2[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0])-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]/2)	#+x
		vec2[1]=-(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2										#+x
		
		vec3[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0])-(coords_ar[nonlin[i]][1]-float(lat_pat[1]))/2)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]/2)	#+x+y
		vec3[1]=-(coords_ar[nonlin[i]][1]-float(lat_pat[1])-coords_ar[lin[j]][1])*np.sqrt(3)/2 	#+x+y
		
		vec4[0]=(coords_ar[nonlin[i]][0]-(coords_ar[nonlin[i]][1]-float(lat_pat[1]))/2)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]/2)	#+y
		vec4[1]=-(coords_ar[nonlin[i]][1]-float(lat_pat[1])-coords_ar[lin[j]][1])*np.sqrt(3)/2		#+y
		vec5[0]=(coords_ar[nonlin[i]][0]-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]+float(lat_pat[0])-coords_ar[lin[j]][1]/2)	#-x
		vec5[1]=-(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2
		vec6[0]=(coords_ar[nonlin[i]][0]-(coords_ar[nonlin[i]][1])/2)-(coords_ar[lin[j]][0]+float(lat_pat[0])-(coords_ar[lin[j]][1]-float(lat_pat[1]))/2)	#-x-y
		vec6[1]=-(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2
		vec7[0]=(coords_ar[nonlin[i]][0]-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]-(coords_ar[lin[j]][1]-float(lat_pat[1]))/2)	#-y
		vec7[1]=-(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2					
		vec8[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0])-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]-(coords_ar[lin[j]][1]-float(lat_pat[1]))/2)	#x-y
		vec8[1]=-(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2	
		vec9[0]=(coords_ar[nonlin[i]][0]-(coords_ar[nonlin[i]][1]-float(lat_pat[1]))/2)-(coords_ar[lin[j]][0]+float(lat_pat[0])-(coords_ar[lin[j]][1])/2)	#-x+y
		vec9[1]=-(coords_ar[nonlin[i]][1]-float(lat_pat[1])-(coords_ar[lin[j]][1]))*np.sqrt(3)/2	
		
		
		
		
		
		#vec1[0]=(coords_ar[nonlin[i]][0]-coords_ar[nonlin[i]][1]*0.5)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]*0.5)
		#vec1[1]=(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2
		#
		#vec2[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0])-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]/2)	#+x
		#vec2[1]=(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2										#+x
		#
		#vec3[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0])-(coords_ar[nonlin[i]][1]-float(lat_pat[1]))/2)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]/2)	#+x+y
		#vec3[1]=(coords_ar[nonlin[i]][1]-float(lat_pat[1])-coords_ar[lin[j]][1])*np.sqrt(3)/2 	#+x+y
		#
		#vec4[0]=(coords_ar[nonlin[i]][0]-(coords_ar[nonlin[i]][1]-float(lat_pat[1]))/2)-(coords_ar[lin[j]][0]-coords_ar[lin[j]][1]/2)	#+y
		#vec4[1]=(coords_ar[nonlin[i]][1]-float(lat_pat[1])-coords_ar[lin[j]][1])*np.sqrt(3)/2		#+y
		#vec5[0]=(coords_ar[nonlin[i]][0]-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]+float(lat_pat[0])-coords_ar[lin[j]][1]/2)	#-x
		#vec5[1]=(coords_ar[nonlin[i]][1]-coords_ar[lin[j]][1])*np.sqrt(3)/2
		#vec6[0]=(coords_ar[nonlin[i]][0]-(coords_ar[nonlin[i]][1])/2)-(coords_ar[lin[j]][0]+float(lat_pat[0])-(coords_ar[lin[j]][1]-float(lat_pat[1]))/2)	#-x-y
		#vec6[1]=(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2
		#vec7[0]=(coords_ar[nonlin[i]][0]-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]-(coords_ar[lin[j]][1]-float(lat_pat[1]))/2)	#-y
		#vec7[1]=(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2					
		#vec8[0]=(coords_ar[nonlin[i]][0]+float(lat_pat[0])-coords_ar[nonlin[i]][1]/2)-(coords_ar[lin[j]][0]-(coords_ar[lin[j]][1]-float(lat_pat[1]))/2)	#x-y
		#vec8[1]=(coords_ar[nonlin[i]][1]-(coords_ar[lin[j]][1]-float(lat_pat[1])))*np.sqrt(3)/2	
		#vec9[0]=(coords_ar[nonlin[i]][0]-(coords_ar[nonlin[i]][1]-float(lat_pat[1]))/2)-(coords_ar[lin[j]][0]+float(lat_pat[0])-(coords_ar[lin[j]][1])/2)	#-x+y
		#vec9[1]=(coords_ar[nonlin[i]][1]-float(lat_pat[1])-(coords_ar[lin[j]][1]))*np.sqrt(3)/2	


																																							
		vec_mag_1=np.sqrt(vec1[1]**2+vec1[0]**2)	
		vec_mag_2=np.sqrt(vec2[1]**2+vec2[0]**2)	
		vec_mag_3=np.sqrt(vec3[1]**2+vec3[0]**2)	
		vec_mag_4=np.sqrt(vec4[1]**2+vec4[0]**2)	
		vec_mag_5=np.sqrt(vec5[1]**2+vec5[0]**2)	
		vec_mag_6=np.sqrt(vec6[1]**2+vec6[0]**2)	
		vec_mag_7=np.sqrt(vec7[1]**2+vec7[0]**2)	
		vec_mag_8=np.sqrt(vec8[1]**2+vec8[0]**2)	
		vec_mag_9=np.sqrt(vec9[1]**2+vec9[0]**2)	
		if (vec_mag_1<20):
			iter=iter+1
			python_script.write('		stk.Edge(%d, _vertex_prototypes[%d], _vertex_prototypes[%d]),\n'%(iter,i,j+len(nonlin)))
		elif(vec_mag_2<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(-1,0,0),\n		),\n'%(iter,i,j+len(nonlin)))
		elif(vec_mag_3<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(-1,1,0),\n		),\n'%(iter,i,j+len(nonlin)))
		elif(vec_mag_4<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(0,1,0),\n		),\n'%(iter,i,j+len(nonlin)))
		elif(vec_mag_5<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(1,0,0),\n		),\n'%(iter,i,j+len(nonlin)))
		elif(vec_mag_6<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(1,-1,0),\n		),\n'%(iter,i,j+len(nonlin)))
		elif(vec_mag_7<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(0,-1,0),\n		),\n'%(iter,i,j+len(nonlin)))	
		elif(vec_mag_8<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(-1,-1,0),\n		),\n'%(iter,i,j+len(nonlin)))	
		elif(vec_mag_9<20):
			iter=iter+1
			python_script.write('		stk.Edge(\n			id=%d,\n	 		vertex1=_vertex_prototypes[%d],\n			vertex2=_vertex_prototypes[%d],\n			periodicity=(1,1,0),\n		),\n'%(iter,i,j+len(nonlin)))																														
		#else:
		#	if(vec1[0]>0.5*float(lat_pat[0])):
		#		vec1[0]=vec1[0]-float(lat_pat[0])
		#		period[0]=-1
		#	elif(vec1[0]<-0.5*float(lat_pat[0])):
		#		vec1[0]=vec1[0]+float(lat_pat[0])
		#		period[0]=1	
		#	if(vec1[1]>0.5*float(lat_pat[1])):
		#		vec1[1]=vec1[1]-float(lat_pat[1])
		#		period[1]=-1
		#	elif(vec1[1]<-0.5*float(lat_pat[1])):
		#		vec1[1]=vec1[1]+float(lat_pat[1])
		#		period[1]=1											
		#	vec_mag=np.sqrt(vec1[1]**2+vec1[0]**2)	
		#	if (vec_mag<30):
		#		iter=iter+1
		#		input3.write('	Edge(\n		id=%d,\n	 	vertex1=_vertex_prototypes[%d],\n		vertex2=_vertex_prototypes[%d],\n		periodicity=(%d,%d,0),\n	),\n'%(iter,i,j,period[0],period[1]))
		
python_script.write('	)')


#		if (loop==1):
#			input3.write('	Edge(\n		id=i,\n		vertex1=_vertex_prototypes[%s],\n		vertex2=_vertex_prototypes[%s],\n		periodicity=(%s),\n	),\n'%(coords_ar[i][0],coords_ar[i][1],coords_ar[i][2]))
#		else:
#			input3.write('	Edge(i,_vertex_prototypes[%s],_vertex_prototypes[%s]),\n'%(coords_ar[i][0],coords_ar[i][1],coords_ar[i][2]))


python_script.write("""\nbb1 = stk.BuildingBlock('Br/N=C/c2ccc(c1ccc(/C=N\Br)cc1)cc2', [stk.BromoFactory()])
bb2 = stk.BuildingBlock('Br/N=C/c1ccc(/C=N/Br)cc1', [stk.BromoFactory()])
bb3 = stk.BuildingBlock('Brc4ccc(c3cc(c1ccc(Br)cc1)cc(c2ccc(Br)cc2)c3)cc4', [stk.BromoFactory()])
bb4 = stk.BuildingBlock('Br/N=C/c3ccc(c2ccc(c1ccc(/C=N/Br)cc1)cc2)cc3', [stk.BromoFactory()])

topology_graph = CUSTOM_PERIODIC_TOPOLOGY(
	building_blocks={\n""")
		



input4=open("bb_%d.txt"%(temp), "w+")
python_script.write('		bb1:(')
for j in range(len(lin)):	
	if (coord_type[lin[j]]=='Mn'):
		python_script.write('%d, '%(j+len(nonlin)))
python_script.write('),\n		bb2:(')
for j in range(len(lin)):	
	if (coord_type[lin[j]]=='Pb'):
		python_script.write('%d, '%(j+len(nonlin)))		
python_script.write('),\n		bb3:(')
for i in range(len(nonlin)):
	python_script.write('%d, '%(i))
python_script.write(')')

python_script.write("""\n	}, 
	lattice_size=(1, 1, 1),
)
construction_result = topology_graph.construct()
cof = stk.ConstructedMolecule.init_from_construction_result(
    construction_result=construction_result,
)
periodic_info = construction_result.get_periodic_info()
cell_matrix = periodic_info.get_cell_matrix()
# Can access all unit-cell parameters.
a = periodic_info.get_a()
b = periodic_info.get_b()
c = periodic_info.get_c()
alpha = periodic_info.get_alpha()
beta = periodic_info.get_beta()
gamma = periodic_info.get_gamma()
unit_cell = stko.UnitCell(
    vector_1=periodic_info.get_vector_1(),
    vector_2=periodic_info.get_vector_2(),
    vector_3=periodic_info.get_vector_3(),
) 

# Write to .pdb file.
writer = stk.PdbWriter()
writer.write(
    molecule=cof,
    path='cof_n_%d_T_%d.pdb',
    periodic_info=periodic_info,
)

writer = stk.XyzWriter()
writer.write(molecule=cof, path='cof.xyz')

writer = stk.MolWriter()
writer.write(molecule=cof, path='cof.mol')

# Define stko UnitCell.
           
gulp_opt = stko.GulpUFFOptimizer(
    gulp_path=(
        'CHANGE MEEE'
    ),
    maxcyc=100,  ## CHANGE
    metal_FF=None,
    conjugate_gradient=False,
)
gulp_opt.assign_FF(cof)
cof, unit_cell = gulp_opt.p_optimize(
    mol=cof,
    unit_cell=unit_cell,
)"""%(n,temp))


