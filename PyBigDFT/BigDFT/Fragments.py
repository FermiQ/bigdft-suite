"""
This module is related to the usage of BigDFT with Fragment-related Quantities.
Input as well as Logfiles might be processed with the classes and methods provided by it.

"""

#: Conversion between Atomic Units and Bohr
AU_to_A=0.52917721092
#: Conversion between Debye and Atomic units
Debye_to_AU = 0.393430307

class XYZfile():
    """
    .. |filename_docs| replace::
         The file which will be created. If None, the file will be eventually dumped in :class:~`sys.stdout`.

    A class associated to a xyz input file as processed by BigDFT

    :param filename: |filename_docs|
    :type filename: string
    :param units: The units of measure of the positions. Allowed avlues are 'atomic' or 'angstroem'
    :type units: string
    """
    def __init__(self,filename=None,units='atomic'):
        self.filename=filename
        self.lines=[]
        self.units=units
        self.fac=1.0
        if units == 'angstroem': self.fac=AU_to_A
    def append(self,array,basename='',names=None,attributes=None):
        """
        Add lines to the file position list

        :param array: list of the atomic positions
        :type array: list of  float triples
        :param basename: base for the name of the atoms
        :type basename: string
        :param names: list of atom names. Will be appended to `basename` if the latter is present
        :type names: list of strings
        :param attributes: list of further attributes to be associated to each of the atoms.
            Will be serialized close to each of the atomic positions
        :type attributes: list of dictionaries
        """
        nm=basename
        for i,r in enumerate(array):
            if names is not None: nm=basename+names[i]
            line=str(nm)
            for t in r:
                line+=' '+str(self.fac*t)
            if attributes is not None: line+=' '+str(attributes[i])
            self.lines.append(line+'\n')
    def dump(self,position='w'):
        """
        Dump the file on the file system if filename has been provided,
        otherwise dump on sys.stdout.

        :param position: filename position statement. Only menaingful for a file dumping.
        :type position: char
        """
        import sys
        f=sys.stdout
        if self.filename is not None: f=open(self.filename,position)
        f.write(str(len(self.lines))+' '+str(self.units)+'\n')
        f.write('# xyz dump \n')
        #then the positions
        for l in self.lines:
            f.write(l)
        if self.filename is not None: f.close()

def open_xyz(filename,nat,unit,comment,position='a'):
    import sys
    f=sys.stdout
    if filename is not None: f=open(filename,position)
    if (position != 'a'):
        f.write(str(nat)+' '+str(unit)+'\n')
        f.write(comment+'\n')
    return f

def close_xyz(f,filename):
    if filename is not None: f.close()

def dump_xyz_positions(f,array,basename='',names=None):
    nm=basename
    for i,r in enumerate(array):
        if names is not None: nm=basename+names[i]
        f.write(str(nm)+' '+str(r[0])+' '+str(r[1])+' '+str(r[2])+'\n')

def xyz_bc_spec(cell):
    """
    Defines the specification for expressing the Boundary Conditions starting from a cell vector.

    :param cell: array of the (orthorhombic) cell. Should be 0.0 on directions with free BC.
       If None is given, the BC are assumed to be Free.
    :type cell: triple of floats or None
    :returns: comment Line of the xyz file specifying the bc
    :rtype: string
    """
    if cell is None:
        return ""
    elif cell[1] == 0.0 and cell[2] != 0.0:
        return "surface "+str(cell[0])+" 0.0 "+str(cell[2])+" "
    elif cell[1] == 0.0 and cell[2] == 0.0:
        return "wire 0.0 0.0 "+cell[2]+" "
    else:
        return "periodic "+str(cell[0])+" "+str(cell[1])+" "+str(cell[2])+" "

def dump_xyz(array,basename='',units='atomic',names=None,filename=None,position='a',comment=None,cell=None):
    """
    Create a BigDFT xyz filename. Duplicates the meaning of the :class:`XYZfile` class.

    :param filename: |filename_docs|
    :type filename: string

    .. todo::
       Remove the duplication of operations in favour or the class.
       Move the related operation into a lower-level module.
    """
    cmt=xyz_bc_spec(cell)
    cmt+=comment if comment is not None else '# xyz dump with basename "'+basename+'"'
    f=open_xyz(filename,len(array),units,cmt,position)
    dump_xyz_positions(f,array,basename=basename,names=names)
    close_xyz(f,filename)

class Lattice():
    """
    Defines the fundamental objects to deal with periodic systems
    """
    def __init__(self,vectors):
        self.vectors=vectors
    def grid(self,origin=[0.0,0.0,0.0],extremes=None,radius=None):
        "produces a set of translation vectors from a given origin"
        import numpy as np
        transl=[]
        g=[[],[],[]] # the grid of discrete translations
        if extremes is not None:
            #print extremes
            for i,d in enumerate(extremes):
                for k in range(d[0],d[1]+1):
                   g[i].append(k)
            #print g
            for i in g[0]:
                arri=np.array(self.vectors[0])*i
                for j in g[1]:
                    arrj=np.array(self.vectors[1])*j+arri
                    for k in g[2]:
                        arrk=np.array(self.vectors[2])*k+arrj
                        vect=np.array(origin)+arrk
                        app=True
                        if radius is not None: app=np.linalg.norm(arrk) < radius
                        if app: transl.append(vect)
        return transl


class RotoTranslation():
    "Define a transformation which can be applied to a group of atoms"
    def __init__(self,pos1,pos2):
        try:
            import wahba
            self.R,self.t,self.J=wahba.rigid_transform_3D(pos1,pos2)
        except Exception,e:
            print 'Error',e
            self.R,self.t,self.J=(None,None,1.0e10)
    def dot(self,pos):
        "Apply the rototranslations on the set of positions provided by pos"
        import wahba as w,numpy as np
        if self.t is None:
            res=w.apply_R(self.R,pos)
        elif self.R is None:
            res=w.apply_t(self.t,pos)
        else:
            res=w.apply_Rt(self.R,self.t,pos)
        return res
    def invert(self):
        self.t=-self.t
        if self.R is not None: self.R=self.R.T

class Translation(RotoTranslation):
    def __init__(self,t):
        import numpy
        self.R=None
        self.t=numpy.mat(t).reshape(3,1)
        self.J=0.0
class Rotation(RotoTranslation):
    def __init__(self,R):
        self.t=None
        self.R=R
        self.J=0.0

class Fragment():
    """
    Introduce the concept of fragment. This is a subportion of the system
    (it may also coincide with the system itself) that is made of atoms.
    Such fragment might have quantities associated to it, like its
    electrostatic multipoles (charge, dipole, etc.) and also geometrical information
    (center of mass, principla axis etc.). A Fragment might also be rototranslated
    and combined with other moieteies to form a :class:`System`.
    """
    protected_keys=['q0','q1','q2','sigma']
    def __init__(self,atomlist=None,id='Unknown',units='AU'):
        self.atoms=[]
        self.id=self.set_id(id)
        self.to_AU=1.0
        if units == 'A': self.to_AU=1.0/AU_to_A
        self.allset=False
        if atomlist is not None:
            for i in atomlist:
                self.append(i)
        self.positions=self.__positions()  #update positions
        self.allset=True
    def __len__(self):
        return len(self.atoms)
    #def __str__(self):
    #    import yaml
    #    return yaml.dump({'Positions': self.atoms,'Properties': {'name': self.id}})
    def set_id(self,id):
	self.id=id
    def xyz(self,filename=None,units='atomic'):
        "Write the fragment positions in a xyz file"
        import numpy as np
        f=XYZfile(filename,units=units)
        names=[ self.element(at) for at in self.atoms]
        posarr=[ np.ravel(r) for r in self.positions]
        f.append(posarr,names=names)
        f.dump()
    def dict(self):
        "Transform the fragment information into a dictionary ready to be put as external potential"
        lat=[]
        for at in self.atoms:
            dat=at.copy()
            dat['r']=list(at[self.__torxyz(at)])
            dat['sym']=self.element(at)
            #assume that the provided charge is alway the net charge
            if 'nzion' in dat: dat.pop('nzion') #for the modification of the conventions
            for k in self.protected_keys:
                if k in at: dat[k]=list(at[k]) #.tolist()
            lat.append(dat)
        return lat
    def append(self,atom=None,sym=None,positions=None):
        if atom is not None:
            self.atoms.append(atom)
        elif sym is not None:
            self.atoms.append({sym: positions})
        if self.allset: self.positions=self.__positions()  #update positions
    def element(self,atom):
        "Provides the name of the element"
        el=self.__torxyz(atom)
        if el == 'r': el=atom['sym']
        return el
    def __torxyz(self,atom):
        "provide the key which contains the positions"
        ks=atom.keys()
        for k in ks:
            if k not in self.protected_keys and type(atom[k])==type([]):
                if len(atom[k])==3: return k
        raise ValueError
    def rxyz(self,atom):
        import numpy as np
        k=self.__torxyz(atom)
        return self.to_AU*np.array(atom[k])
    def __positions(self):
        import numpy
        return numpy.mat([ self.rxyz(at) for at in self.atoms ])
    def centroid(self):
        import numpy
        return numpy.ravel(numpy.mean(self.positions, axis=0))
    def center_of_charge(self,zion):
        cc=0.0
        qtot=0.0
        for at in self.atoms:
            netcharge=at.get('q0')[0]
            sym=self.element(at)
            zcharge=zion[sym]
            elcharge=zcharge-netcharge
            cc+=elcharge*self.rxyz(at)
            qtot+=elcharge
        return cc/qtot
    def transform(self,Rt): #R=None,t=None):
        "Apply a rototranslation of the fragment positions"
        import numpy as np
        self.positions=Rt.dot(self.positions)
        #import wahba as w,numpy as np
        #if t is None:
        #    self.positions=w.apply_R(R,self.positions)
        #elif R is None:
        #    self.positions=w.apply_t(t,self.positions)
        #else:
        #    self.positions=w.apply_Rt(R,t,self.positions)
        #then replace the correct positions at the atoms
        for at,r in zip(self.atoms,self.positions):
            k=self.__torxyz(at)
            at[k]=np.ravel(r).tolist()
        #further treatments have to be added for the atomic multipoles
        #they should be transfomed accordingly, up the the dipoles at least
    def line_up(self):
        "Align the principal axis of inertia of the fragments along the coordinate axis. Also shift the fragment such as its centroid is zero."
        import numpy
        Shift=Translation(self.centroid())
        Shift.invert()
        self.transform(Shift)
        #now the centroid is zero
        I=self.ellipsoid()
        w,v=numpy.linalg.eig(I)
        Redress=Rotation(v.T)
        self.transform(Redress)
        #now the principal axis of inertia are on the coordinate axis
    def ellipsoid(self,center=0.0):
        import numpy as np
        I=np.mat(np.zeros(9).reshape(3,3))
        for at in self.atoms:
            rxyz=self.rxyz(at)-center
            I[0,0]+=rxyz[0]**2 #rxyz[1]**2+rxyz[2]**2
            I[1,1]+=rxyz[1]**2 #rxyz[0]**2+rxyz[2]**2
            I[2,2]+=rxyz[2]**2 #rxyz[1]**2+rxyz[0]**2
            I[0,1]+=rxyz[1]*rxyz[0]
            I[1,0]+=rxyz[1]*rxyz[0]
            I[0,2]+=rxyz[2]*rxyz[0]
            I[2,0]+=rxyz[2]*rxyz[0]
            I[1,2]+=rxyz[2]*rxyz[1]
            I[2,1]+=rxyz[2]*rxyz[1]
        return I
    def q0(self,atom):
        "Provides the charge of the atom"
        charge=atom.get('q0')
        if charge is not None: charge=charge[0]
        return charge
    def q1(self,atom):
        "Provides the dipole of the atom"
        import numpy as np
        dipole=atom.get('q1') #they are (so far) always given in AU
        if dipole is not None: dipole=np.array([dipole[2],dipole[0],dipole[1]])
        return dipole
    def Q(self,atom=None):
        "Fragment Monopole, when the information is available. Select the element if possible"
        charge=0
        found=False
        for at in self.atoms:
            q0=self.q0(at)
            if q0 is not None and (atom is None or self.element(at)==atom):
                found=True
                charge+=q0
        if found:
            return charge
        else:
            return None
    def d0(self,center=None):
        "Fragment dipole, calculated only from the atomic charges"
        #one might added a treatment for non-neutral fragments
        #but if the center of charge is used the d0 value is zero
        import numpy as np
        if center is not None:
            cxyz=center
        else:
            cxyz=np.ravel(self.centroid())
        d0=np.zeros(3)
        found=False
        for at in self.atoms:
            if self.q0(at) is not None:
                found=True
                d0+=at['q0'][0]*(self.rxyz(at)-cxyz)
        if found:
            return d0
        else:
            return None
    def d1(self,center=None):
        "Fragment dipole including the atomic dipoles"
        import numpy as np
        d1=np.zeros(3)
        dtot=self.d0(center)
        if dtot is None: return dtot
        found=False
        for at in self.atoms:
            q1=self.q1(at)
            if q1 is not None:
                found=True
                d1+=q1
        if found:
            return d1+dtot
        else:
            return None


class System():
    "A system is defined by a collection of Fragments. It might be given by one single fragment"
    def __init__(self,mp_dict=None,xyz=None,nat_reference=None,units='AU',transformations=None,reference_fragments=None):
        self.fragments=[]
        self.CMs=[]
        self._get_units(units)
        if xyz is not None: self.fill_from_xyz(xyz,nat_reference)
        if mp_dict is not None: self.fill_from_mp_dict(mp_dict,nat_reference)
        if transformations is not None: self.recompose(transformations,reference_fragments)
    def __len__(self):
        return sum([len(frag) for frag in self.fragments])
    def _get_units(self,unt):
        self.units=unt
        if unt=='angstroem':
            self.units='A'
        elif unt=='atomic' or unt=='bohr':
            self.units='AU'
    def _bigdft_units(self):
        return 'angstroem' if self.units=='A' else 'atomic'
    def fill_from_xyz(self,file,nat_reference):
        "Import the fragment information from a xyz file"
        fil=open(file,'r')
        nat=0
        iat=0
        frag=None
        for l in fil:
            try:
                pos=l.split()
                if len(pos) <=2: #these are the number of atoms
                    nt=int(pos[0])
                    nat-=nt
                    if len(pos)==2:
                        unt=pos[1]
                        self._get_units(unt)
                    if frag is not None: self.append(frag)
                    frag=Fragment(units=self.units)
                    iat=0
                elif len(pos)>0:
                    if nat_reference is not None and iat == nat_reference: #we should break the fragment, alternative strategy
                        if frag is not None: self.append(frag)
                        frag=Fragment(units=self.units)
                        iat=0
                    frag.append({pos[0]: map(float,pos[1:])})
                    nat+=1
                    iat+=1
            except Exception,e:
                print 'Warning, line not parsed: "',l,e,'"'
        if iat != 0: self.append(frag) #append the remaining fragment
    def fill_from_mp_dict(self,mpd,nat_reference=None):
        "Fill the System from a dictionary of multipole coefficients"
        frag=Fragment(units=self.units)
        iat=0
        for at in mpd:
            #frag.append(sym=at['sym'],positions=at['r'])
            frag.append(at)
            iat+=1
            if nat_reference is not None and iat == nat_reference:
                if len(frag) !=0: self.append(frag)
                frag=Fragment(units=self.units)
                iat=0
        if nat_reference is None: self.append(frag) #case of one single fragment
    def xyz(self,filename=None,units='atomic'):
        import numpy as np
        f=XYZfile(filename,units)
        for i,frag in enumerate(self.fragments):
            names=[ frag.element(at) for at in frag.atoms]
            posarr=[ np.ravel(r) for r in frag.positions]
            f.append(posarr,names=names,attributes=[{'frag': [frag.id, i+1]}]*len(frag))
        f.dump()
    def dict(self,filename=None):
        atoms=[]
        for f in self.fragments:
            atoms+=f.dict()
        #if self.units != 'A':
        #    print 'Dictionary version not available if the system is given in AU'
        #    raise Exception
        dc={'units': self._bigdft_units(),'global monopole': float(self.Q()), 'values': atoms}
        return dc
    def append(self,frag):
        "Append a fragment to the System class"
        assert isinstance(frag,Fragment)
        self.fragments.append(frag)
        self.CMs.append(frag.centroid())  #update center of mass
    def pop(self,ifrag):
        "Pop the fragment ifrag from the list of fragments"
        self.CMs.pop(ifrag)
        return self.fragments.pop(ifrag)
    def centroid(self):
        "Center of mass of the system"
        import numpy
        return numpy.mean(self.CMs, axis=0)
    def central_fragment(self):
        "Returns the fragment whose center of mass is closest to the centroid"
        import numpy as np
        return np.argmin([ np.dot(dd,dd.T) for dd in (self.CMs - self.centroid())])
        #return self.fragments[imin]
    def fragment_transformation(self,frag1,frag2):
        "returns the transformation among fragments if exists"
        return RotoTranslation(frag1.positions,frag2.positions)
        #try:
        #    import wahba
        #    roto,translation,J=wahba.rigid_transform_3D(frag1.positions,frag2.positions)
        #except Exception,e:
        #    print 'Error',e
        #    roto,translation,J=(None,None,1.0e10)
        #return roto,translation,J
    def decompose(self,reference_fragments):
        "Decompose the system into reference fragments"
        assert type(reference_fragments) == type([])
        self.decomposition=[]
        for frag in self.fragments:
            transf=[]
            Js=[]
            for ref in reference_fragments:
                #r,t,j=self.fragment_transformation(ref,frag)
                RT=self.fragment_transformation(ref,frag)
                transf.append({'RT': RT})#{'R':r,'t':t})
                #Js.append(j)
            #choose the minimal one
            import numpy
            Jchosen=numpy.argmin([rt['RT'].J for rt in transf])
            ref=transf[Jchosen]
            ref['ref']=reference_fragments[Jchosen]
            #ref['J']=Js[Jchosen]
            ref['id']=Jchosen
            self.decomposition.append(ref)
    def recompose(self,transformations=None,reference_fragments=None):
        "Rebuild the system from a set of transformations"
        import copy,numpy as np
        if transformations is not None:
            RT=transformations
            self.decomposition=RT
        else:
            RT=self.decomposition
        self.fragments=[]
        self.CMs=[]
	self.templates=[]
        for item in RT:
            if reference_fragments:
                idf=item['id']
		template=reference_fragments[idf]
            else:
		template=item['ref']
	    frag=copy.deepcopy(template)
	    self.templates.append(template)
            #frag.transform(item['R'],item['t'])
            frag.transform(item['RT'])
            self.append(frag)
    def Q(self):
            "Provides the global monopole of the system given as a sum of the monopoles of the atoms"
            return sum([ f.Q() for f in self.fragments])
    def fragdict(self):
	""" Provides the value of the dictionary fragment to be used for the inputfile in a fragment calculation """
	refs=[]
	for t in self.templates:
	    if t not in refs: refs.append(t)
	#generate the fragments id that have to be put into the input posinp
	allfrags=find_reference_fragment(refs,self.templates)
	fragdict={}
	for t,r in zip(refs,allfrags):
	     fragdict[t.id]=r
	return fragdict

# create the directory of the template file
def prepare_fragment_inputs(name,directory='.',flavour='Isolated',system=None,template=None,template_dir=None,template_name=None):
    import shutil,os,yaml
    from futile.Utils import ensure_dir
    dirct=directory
    ensure_dir(dirct)
    posinp=name+'.xyz'
    if template is not None: template.xyz(filename=posinp)
    input_dict={'posinp': posinp,'import': 'linear_laura'}
    if system is not None:
	input_dict['import']=['linear_laura','linear_fragments'] if flavour!='Embedded' else 'linear_laura'
	input_dict['frag']=system.fragdict()
	system.xyz(filename=posinp)
	datadir=os.path.join(dirct,'data-'+name)
	tempdatadir='data-'+template_name
	ensure_dir(datadir)
	datatemplate=os.path.join(datadir,tempdatadir)
	if flavour=='Embedded':
	   ensure_dir(datatemplate)
	elif not os.path.exists(datatemplate):
	   os.symlink(os.path.abspath(os.path.join(template_dir,tempdatadir)),datatemplate)
	for ext in ['.xyz','.yaml']:
	   f=os.path.join(template_dir,template_name+ext)
	   if os.path.exists(f): shutil.copyfile(src=f,dst=os.path.join(datadir,template_name+ext))
    if dirct != '.': shutil.copyfile(src=posinp,dst=os.path.join(dirct,posinp))
    f=open(os.path.join(dirct,name+'.yaml'),'w')
    f.write(yaml.dump(input_dict))
    f.close()
def find_reference_fragment(refs,transformed):
    """generate the list of fragments that are associates to a given list of templates"""
    import numpy as np
    where_for_frags=[]
    for r in refs:
        where_for_frags.append([ i+1 for i,t in enumerate(transformed) if t == r])
    return where_for_frags

def frag_average(ref,flist,clean_monopole=True):
    "Perform the average in attributes of the fragments provided by the list, with the position of the first fragment"
    import copy,numpy as np
    keys=['q0','q1','q2']
    navg=len(flist)
    if navg==0: return ref
    favg=copy.deepcopy(ref)
    qtot=0.0
    for i,at in enumerate(favg.atoms):
        #form a fragment which has the positions of the references and
        #neutral total monopole if asked for.
        for k in keys:
            population=[ f.atoms[i][k] for f in flist ]
            vals=np.mean(population,axis=0)
            st=np.std(population,axis=0)
            at[k]=vals
            at['sigma'+k]=st
            #print 'test',k,np.max(abs(at['sigma'+k]/at[k]))
        qtot+=at['q0'][0]
    qtot/=float(len(favg.atoms))
    for at in favg.atoms:
        at['q0'][0]-=qtot
        #print 'retest',i,at
    return favg

def distance(i,j):
    "Distance between fragments, defined as distance between center of mass"
    import numpy
    vec=i.centroid()-j.centroid()
    return numpy.sqrt(numpy.dot(vec,vec.T))

def build_transformations(RTlist,ref):
    "Build the transformation list to be passed to the systems' creation via the transformation keyword"
    return [ {'RT': RT,'ref':ref} for RT in RTlist]

def wahba_fragment(frag1,frag2):
    "Solve the wahba's problem among fragments, deprecated routine"
    import wahba #should be cleaned
    #For each of the fragment build the list of the coordinates
    roto,translation,J=wahba.rigid_transform_3D(frag1.positions,frag2.positions)
    return roto,translation,J

def rotot_collection(ref_frag,lookup,fragments):
    "Deprecated routine, only for backward compatibility"
    W=[]
    for f in lookup:
        refF=lookup[ref_frag]
        roto,translation,J=wahba_fragment(fragments[f],fragments[refF])
        if (J > 1.e-12):
            print 'Error',f,J,refF
            #try with the second ref
            refF2=lookup[ref_frag+1]
            roto2,translation2,J2=wahba_fragment(fragments[f],fragments[refF2])
            #print 'Error Now',f,J2,refF2
            if (J2< J):
                roto=roto2
                translation=translation2
                J=J2
                refF=refF2
        W.append({'R':roto,'t':translation,'J':J,'ref':refF})
    return W

if __name__ == '__main__':
    #extract fragments
    import sys,numpy
    one1=System(xyz='one-1.xyz')
    print 'Parsed',len(one1.fragments)
    print one1.xyz()
    two=System(xyz='two.xyz',nat_reference=len(one1),units='A')
    two.decompose(one1.fragments)
    trans=two.decomposition
    PC1=one1.fragments[0]
    print 'one',PC1.centroid()
    for frag,t in zip(two.fragments,trans):
        print 'ff',frag.centroid()
        print 't',t['t']
        print 'R',t['R']
    print two.CMs
    print trans
    two2=System(transformations=trans)
    two2.xyz('two-2.xyz',units='angstroem')
    #now rigidify the big case scenario
    #read lattice coordinates
    fil=open('lattice.txt','r')
    acell=[eval(l.strip('\r\n')) for l in fil]
    latt=Lattice(acell)
    print latt.vectors
    print 0.5*numpy.array(latt.vectors)
    print two.CMs[1]-two.CMs[0]
    #find the positions of the center of mass of the big system
    bigold=System(xyz='BigCase.xyz',nat_reference=36)
    icen=bigold.central_fragment()
    print icen,bigold.CMs[icen],bigold.CMs[icen-1]
    samples=[bigold.CMs[icen],bigold.CMs[icen-1]]
    extremes=[[-5,5],[-5,5],[-1,1]]
    grid=[]
    for oxyz in samples: #two.CMs:
      grid += latt.grid(extremes=extremes,origin=oxyz)
      #grid centroid
    cent=numpy.ravel(numpy.mean(numpy.mat(grid), axis=0))
    trans=[]
    limit=len(grid)/len(two.CMs)
    for i,d in enumerate(grid):
        ref=two2.fragments[i/limit]
        if numpy.linalg.norm(d-cent) < 25.0:
            trans.append({'t':numpy.mat(d).T/AU_to_A,'ref':ref,'R':None})

    big=System(transformations=trans)
    big.xyz('Bigcase-2.xyz',units='angstroem')
    cents=XYZfile('Bigcase2-centroids.xyz',units='angstroem')
    cents.append(big.CMs,basename='Cen')
    cents.dump()
    icen = big.central_fragment()
    print 'the central fragment is',icen
    #find the atoms
    iat=0
    for i,f in enumerate(big.fragments):
        if i==icen: print 'from',iat+1
        iat+=len(f)
        if i==icen:
            print 'to',iat+1
            break
    exit(0)
    filename=sys.argv[1]
    limit=36 #maximum value of each fragment
    fragments=[]
    #try to initialize with the class

    fil=open(filename,'r')
    count=0
    nat=0
    iat=0
    frag=None
    for l in fil:
        count+=1
        try:
          pos=l.split()
          if len(pos) ==1: #these are the number of atoms
            nt=int(pos[0])
            nat-=nt
            if frag is not None: fragments.append(frag)
            frag=Fragment()
            iat=0
          elif len(pos)>0:
            if iat == limit: #we should break the fragment, alternative strategy
                if frag is not None: fragments.append(frag)
                frag=Fragment()
                iat=0
            frag.append({pos[0]: map(float,pos[1:])})
            nat+=1
            iat+=1
        except Exception,e:
            print 'error',l,e
            break

    print 'calculation finished',len(fragments),'balance',nat

    #find the F4TCNQ
    F4TCNQs=[]
    PCs=[]
    CMs=[] #ordered center of mass
    for f,frag in enumerate(fragments):
        CMs.append(frag.centroid())
        if len(frag) < 36:
            F4TCNQs.append(f)
        else:
            PCs.append(f)
    #find the central molecule
    centroid= numpy.mean(CMs, axis=0)

    print 'species identified:',len(F4TCNQs),'F4TCNQ and',len(PCs),' pentacenes, tot',len(fragments),len(CMs)

    #now append the fragments to the System class
    stm=System(xyz=filename,nat_reference=36)
    print len(stm.fragments),'before'
    for frag in fragments:
        stm.append(frag)
    print len(stm.fragments),'after'
    print centroid,[i for i in CMs[0]],[i for i in centroid],stm.centroid(),stm.central_fragment()

    refF4=0
    refPC=0
    #try:
    #    icen=F4TCNQs.index(imin)
    #    refF4=icen
    #except:
    #    icen=PCs.index(imin)
    #    refPC=icen

    #check if now all the atoms are the rototranslation of the same fragment and find the transformation
    W_F4=rotot_collection(refF4,F4TCNQs,fragments)
    W_PEN=rotot_collection(refPC,PCs,fragments)


    #print CMs
    #search the NN of each of the F4TCNQs
    DFP=[]
    DFF=[]
    for f in F4TCNQs:
        DFP.append(numpy.array([distance(f,p) for p in PCs]))
        DFF.append(numpy.array([distance(f,p) for p in F4TCNQs]))


    #Then we can classify the attributes of each pentacene accordingly to the limiting distance

    threshold=10.0

    for i,f in enumerate(F4TCNQs):
        import yaml
        if (i==0):
            print yaml.dump(fragments[f])
            print 'test wahba'
            roto,translation,J=wahba_fragment(fragments[f],fragments[f])
            print roto
            print translation
            print 'Rototranslation Error',J
        iPC=0
        for dist in DFP[i]:
            if dist < threshold:
                iPC+=1
        iFF=0
        for dist in DFF[i]:
            if dist < threshold and dist !=0:
                iFF+=1
        print i,iFF,iPC
