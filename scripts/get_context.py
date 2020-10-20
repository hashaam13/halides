from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from biopandas.pdb import PandasPdb
import re
import os
import argparse
import numpy as np
import pandas as pd
import copy
import itertools
# import freesasa
import subprocess
import math


def ang(Coordinate): # get angles
             
  Angles=[]

  for  k in range(2,len(Coordinate)): 

       a=np.array(Coordinate[0])
       b=np.array(Coordinate[1])
       c=np.array(Coordinate[k])
       ab=b-a
       ac=c-a
       cosine_angle = np.dot(ab,ac)/(np.linalg.norm(ab)*np.linalg.norm(ac))
       cosine_angle=round(cosine_angle,15) # remove the error
       angle=np.arccos(cosine_angle) 
       if  np.degrees(angle) !=0:
           Angles+=[np.degrees(angle)] # fill the list of angles
       else:
           Angles+=[-1]

  return(Angles)
 

def site_filter_dist(N,modern_subset1,sites):
  atom_name=0
  residue_name=0
  if N>=1:
   for i in range(N):
        A=sites.loc[sites[f'dist_{N}']<5, [f'atom_name_{N}']]
        A=list(A[f'atom_name_{N}'])
        B=sites.loc[sites[f'dist_{N}']<5, [f'atom_name_{i}']]
        B=list(B[f'atom_name_{i}'])
        result=list(set(A) ^ set(B))
        
        if len(result)==0:
           atom_name=1
   for i in range(N):
        C=sites.loc[sites[f'dist_{N}']<5, [f'resi_name_{N}']]
        C=list(C[f'resi_name_{N}'])
        D=sites.loc[sites[f'dist_{N}']<5, [f'resi_name_{i}']]
        D=list(D[f'resi_name_{i}'])
        result=list(set(C) ^ set(D))
        
        if len(result)==0:
          residue_name=1
   for k in range(0,3*N,3):
    same_degree=0
    stop=sites[f'dist_{N}']<5
    
    for g in range(len(modern_subset1)):
     if sites.iloc[g][3*N]<5:
      if abs(sites.iloc[g][3*N]-sites.iloc[g][k])<0.5:
        same_degree+=1
      if (same_degree==stop.sum()) and (atom_name==1) and (residue_name==1) :
        return True

def  Type_ligands(modern_subset):
   Type_ligands=[]
   Type_atoms_DNA=['DA','DG','DT','DC']
   Type_atoms_RNA=['A','G','C','U']
   for i in modern_subset.values:
     #print(i[5])
     if  i[0]=='HETATM' and i[5] not in ['HOH','AHOH','BHOH']:
         Type_ligands += ['MOLECULE']
     elif i[0]=='HETATM' and i[5]  in ['HOH','AHOH','BHOH']:
         Type_ligands += ['WATER']
     else:
         if i[5] in Type_atoms_DNA:
            Type_ligands += ['DNA']
         elif i[5] in Type_atoms_RNA:
             Type_ligands += ['RNA']
         else:
            Type_ligands += ['PROTEIN']
   #print(Type_ligands)
   return(Type_ligands)

def Make_dssp():
      ref = {
         	'A' : 'ALA', 
		'R' : 'ARG',
		'N' : 'ASN',
		'D' : 'ASP',
		'B' : 'ASX',
		'C' : 'CYS',
		'E' : 'GLU',
	 	'Q' : 'GLN',
		'Z' : 'GLX', 
	 	'G' : 'GLY', 
		'H' : 'HIS',
		'I' : 'ILE',
	 	'L' : 'LEU', 
 		'K' : 'LYS', 
		'M' : 'MET',
	 	'F' : 'PHE', 
	 	'P' : 'PRO',
	 	'S' : 'SER', 
	 	'T' : 'THR',
 		'W' : 'TRP',
	 	'Y' : 'TYR', 
 	 	'V' : 'VAL',
                'X' : '---'} 
      dssp_dict ={}
      dssp_dict['X'] = np.NaN
      p = PDBParser()
      structure = p.get_structure('model', './current_pdb.txt')
      mod = structure[0]
      dssp = DSSP(mod, './current_pdb.txt')
      for i in range(len(dssp)):
             a_key = list(dssp.keys())[i]
             dssp_dict[str(a_key[0])+str(a_key[1][1])+str(ref[dssp[a_key][1]])] = dssp[a_key][2]
      return(dssp_dict)

def water_del(f):
    for k in range(len(f)):
             if (f[k][0:6]=='HETATM') and (f[k][16:20]==' HOH' or f[k][16:20]=='AHOH' or f[k][16:20]=='BHOH'):
                       f[k]=''
    return(f)
def main(args):
  F_r=4*(1.19+1.4)**2*math.pi
  CL_r=4*(1.67+1.4)**2*math.pi
  BR_r=4*(1.82+1.4)**2*math.pi
  I_r=4*(2.06+1.4)**2*math.pi
  # try:
    # os.makedirs(args.output_dir)
  # except:
  #   pass
  all_sites=0
  all_site_deleted=0
  all_number_PROTEIN_sites=0
  all_number_DNA_sites=0
  all_number_RNA_sites=0
  all_number_OTHER_sites=0

  with open(args.input_filter, 'r') as f:
   text=f.readlines()
   base_list=[]
   base_list+=text
   base_list=[line.rstrip() for line in base_list]
   base_list=[i for i in base_list if i !='']
   # print (base_list)
   high_resolution=0
   x_ray=base_list.count('X-RAY DIFFRACTION\n')
   NMR=base_list.count('NMR\n')
   for i in range(2,len(base_list),4): 
    if (float(base_list[i-1])<=2):
     high_resolution+=1
     if(base_list[i+1]!='NMR'):
      if args.input_type == 'pdb_id':
       try:
          struct = PandasPdb().fetch_pdb(f'{base_list[i]}')
       except:
          continue
       model_name = args.input
       with open('current_pdb.txt', 'w') as w:
                       w.write(f'{struct.pdb_text}')        
       model_name = base_list[i]
       # print(model_name)
       with open('current_pdb.txt', 'r') as w:
                       f=w.readlines()
       if args.water == 'n':
          water_del(f)
       with open('current_pdb.txt', 'w') as w:
                   line = ''.join(f)
                   w.write(line)
       struct = struct.read_pdb('current_pdb.txt')    

      elif args.input_type == 'structure':
        struct = PandasPdb()

        struct = struct.read_pdb(f'{args.input_struct}/pdb{base_list[i].lower()}.ent')
        # print(f'{args.input}/pdb{base_list[i].lower()}.ent')
        model_name = re.search('[\d\w]+$', struct.header).group()
        with open (f'{args.input_struct}/pdb{base_list[i].lower()}.ent', 'r') as pdb1:
                      f=pdb1.readlines()
        if args.water == 'n':              
           water_del(f)
        with open('current_pdb.txt', 'w') as w:
                  line = ''.join(f)
                  w.write(line)
        struct = struct.read_pdb('current_pdb.txt')
      try:
                  resolution = float(re.search("REMARK\s+2\s+RESOLUTION\.\s+(\d+\.\d+)", struct.pdb_text).group(1))
      except:
                  resolution = 100
        # print(resolution)

       

      halide_type = args.input_filter.split('/')[-1].split('_')[-1].split('.')[0]
      halide_atoms = struct.df['HETATM'][struct.df['HETATM']['atom_name'] == halide_type]

      if args.ligands=='y':
          ligand_atoms=struct.df['HETATM']
          protein_atoms=struct.df['ATOM'] # make the subset
          modern_df=pd.concat([protein_atoms,ligand_atoms], axis=0)
          modern_df.index = np.arange(len(modern_df))
  
      elif args.ligands=='n':
          modern_df=struct.df['ATOM'] 
   
      
      if args.sec_struct == 'y':
        try:
          dssp_dict = Make_dssp() # making dssp dict
        except:
          dssp_dict = {'X' : np.NaN}
      else:
          dssp_dict = {'X' : np.NaN}

      dict_of_subsets = {}
      dict_of_subsets1 = {}
      S=0 # halide counter for using freesasa
      sites=pd.DataFrame()
      N=-1 #sites counter for sites filter
      site_deleted=0
      number_DNA_sites=0
      number_RNA_sites=0
      number_PROTEIN_sites=0
      number_OTHER_sites=0
      if args.water == 'y':
         water_del(f)
      with open('current_pdb.txt', 'w') as w:
                   line = ''.join(f)
                   w.write(line)
      out=subprocess.check_output(["freesasa","--format=pdb", "./current_pdb.txt", "-H"],stderr=subprocess.STDOUT, encoding='utf-8')
      with open ('current_pdb.txt', 'w') as w:
           w.write(out)
      with open('current_pdb.txt','r') as w:
           g=w.readlines()
      for i in halide_atoms.values:
                    sites_dist=pd.DataFrame()
                    halide_atoms.index=np.arange(len(halide_atoms))
                    Halide_humber= halide_atoms[halide_atoms.index==S].values[0][1]
                    S+=1
                    N+=1
                    for k in range(len(g)):
                                if ((g[k][0:6]=='HETATM') and (int(g[k][6:11])==int(Halide_humber))):
                                   asa = float(g[k][60:67])
                                   if halide_type  == 'F':
                                      asa_r=asa/F_r
                                   if halide_type == 'CL':
                                      asa_r=asa/CL_r
                                   if halide_type == 'BR':
                                      asa_r=asa/BR_r
                                   if halide_type == 'I':
                                      asa_r=asa/I_r
                   
                    Coordinate=[] # list of coordinates М[0]= halide coordinates М[1] the nearest atom's coordinates
                    if args.ligands=='y':
                     
                        dist = struct.distance(xyz=tuple(i[11:14]))
                        dist.index=np.arange(len(dist))
                    elif args.ligands=='n':
                        dist = struct.distance(xyz=tuple(i[11:14]), records=('ATOM'))
                    modern_df['dist']=dist # add distanse to subset
                    
                    if args.C == 1:
                           

                           modern_subset = modern_df[(modern_df.dist < args.angstrem_radius+0.5) & (modern_df.dist!=0)] # halide neighbors


                    elif args.C == 2:

                           modern_df1=modern_df[(modern_df.dist < args.angstrem_radius+0.5) & (modern_df.dist!=0)]
                           modern_subset=modern_df1.loc[~modern_df1['element_symbol'].isin(['C','H'])]

  
                    elif args.C == 3:
                           
                           modern_df1 =modern_df[(modern_df.dist < args.angstrem_radius+0.5) & (modern_df.dist!=0)]
                           modern_subset=modern_df1.loc[~modern_df1['element_symbol'].isin(['C'])]

                    elif args.C == 4:
                           modern_df1 =modern_df[(modern_df.dist < args.angstrem_radius+0.5) & (modern_df.dist!=0)]
                           atoms=list(map("".join,itertools.permutations('BDEGHZ',1)))
                           atoms1=list(map("".join,itertools.permutations('BDEGHZ123',2)))
                           atoms2=atoms+atoms1
                           atom2=['C'+ atom for atom in atoms2] +['O','C']
                           modern_subset=modern_df1.loc[~modern_df1['atom_name'].isin(atom2)]
                    
                    if len(modern_subset[modern_subset.dist<args.angstrem_radius])==0:
                                N-=1
                                site_deleted+=1
                                continue

                    modern_subset2=modern_subset.copy(deep=True)
                    modern_subset2.loc[:,['x_coord','y_coord','z_coord']]=modern_subset[['x_coord','y_coord','z_coord']]-[i[11],i[12],i[13]]      
                    
                    xyz=i[11:14] # halide coordinate
                    for k in range(len(xyz)):
                               xyz[k]=0
                 
                    Coordinate+=[xyz] # add coordinates
                    try:
                              nearest=modern_subset2.loc[modern_subset['dist']==min(modern_subset[modern_subset.record_name == 'ATOM']['dist'])].values[0][[3,5,11,12,13,20,21]] #define the nearest atom
      
                              Coordinate+=[nearest[2:5]] # add the nearest atom's coordinates
                    except:
                         continue
                    for n in modern_subset2.values:
                           xyz2= n[11:14]
                           Coordinate+=[xyz2] # add coordinates
  
                    modern_subset1=modern_subset2.copy(deep=True)
                    modern_subset1['angles']=ang(Coordinate) # add angles to subset
                    modern_subset1['Type_ligands']=Type_ligands(modern_subset)
                    modern_subset1.index = np.arange(len(modern_subset1))

                    if args.prot== 'n':
                       if 'DNA' in modern_subset1.values:
                         number_DNA_sites+=1
                         
                       elif 'RNA' in modern_subset1.values:
                         number_RNA_sites+=1
                         
                       elif 'MOLECULE' in modern_subset1.values:
                         number_OTHER_sites+=1
                       else:
                         number_PROTEIN_sites+=1
                    elif args.prot== 'y':
                       if 'DNA' in modern_subset1.values:
                         number_DNA_sites+=1
                         site_deleted+=1
                         N-=1
                         continue
                       elif 'RNA' in modern_subset1.values:
                         number_RNA_sites+=1
                         site_deleted+=1
                         N-=1
                         continue
                       elif 'MOLECULE' in modern_subset1.values:
                         number_OTHER_sites+=1
                         site_deleted+=1
                         N-=1
                         continue
                       else:
                         number_PROTEIN_sites+=1

                    modern_subset3=modern_subset1.loc[~modern_subset1['Type_ligands'].isin(['MOLECULE','RNA','DNA','WATER'])]
                    modern_subset3.index = np.arange(len(modern_subset3))
                    if len(modern_subset3[modern_subset3.dist<args.angstrem_radius])==0:
                                N-=1
                                site_deleted+=1
                                continue
                    sites_dist[f'dist_{N}']=modern_subset3['dist']
                    sites_dist[f'atom_name_{N}']=modern_subset3['atom_name']
                    sites_dist[f'resi_name_{N}']=modern_subset3['residue_name']
                    sites_dist_sort= sites_dist.sort_values(f'dist_{N}')
                    sites_dist_sort.index = np.arange(len(sites_dist_sort))
                    sites=pd.concat([sites,sites_dist_sort], axis=1)
                    try:
                      if  site_filter_dist(N,modern_subset3['dist'],sites):
                         
                         site_deleted+=1
                         # print(f'{halide_type} atom is skipped, similar haligen site around 5A have already been got')
                         
                         continue
                    except:
                         N=-1

                   
                    dssp_cod = list(map(lambda x : dssp_dict[x] , [str(modern_subset3.at[x,'chain_id'])+str(modern_subset3.at[x,'residue_number'])+str(modern_subset3.at[x,'residue_name'])  if str(modern_subset3.at[x,'chain_id'])+str(modern_subset3.at[x,'residue_number'])+str(modern_subset3.at[x,'residue_name']) in dssp_dict.keys() else 'X'  for x in range (len(modern_subset3))  ]))
                    dssp_cod_df = pd.DataFrame({'dssp_cod' : dssp_cod})
                    modern_subset_exit = pd.concat([ modern_subset1,dssp_cod_df], ignore_index=True, axis=1)
                    modern_subset_exit = modern_subset_exit[modern_subset_exit[21] < args.angstrem_radius]
                    site_count= len(modern_subset_exit.groupby(7).size())
                    print(modern_subset_exit)
                    dict_of_subsets[f'{model_name}:{asa_r}:{resolution}:{i[3]}:{i[1]}:{site_count}'] =\
                    [(f'{j[3]}:{j[5]}:{"%.3f"% j[21]}:{"%.3f"% j[22]}:{j[23]}:{j[7]}:{j[24]}') for j in modern_subset_exit.values]
                    dict_of_subsets1[f'{model_name}:{i[3]}:{i[1]}:{site_count}'] =\
                    [(f'{j[3]}:{j[5]}') for j in modern_subset_exit.values]
      #print(f'get {len(halide_atoms)} {halide_type} sites, {site_deleted} sites were deleted. In selected sites - protein sites: {number_PROTEIN_sites}, DNA sites: {number_DNA_sites}, RNA sites: {number_RNA_sites}, other sites: {number_OTHER_sites}')
      all_sites+=len(halide_atoms)
      all_site_deleted+=site_deleted
      all_number_PROTEIN_sites+=number_PROTEIN_sites
      all_number_DNA_sites+=number_DNA_sites
      all_number_RNA_sites+=number_RNA_sites
      all_number_OTHER_sites+=number_OTHER_sites

      with open(f'{args.output_full}', 'a') as w:
            for k,v in dict_of_subsets.items():
              w.write(f'{k}\t')
              for i in range(len(v)):
                if i == len(v)-1:
                  w.write(f'{v[i]}')
                else:
                  w.write(f'{v[i]},')
              w.write('\n')

      with open(f'{args.output_AA}', 'a') as w:
            for k,v in dict_of_subsets1.items():
              w.write(f'{k}\t')
              for i in range(len(v)):
                if i == len(v)-1:
                  w.write(f'{v[i]}')
                else:
                  w.write(f'{v[i]},')
              w.write('\n')

      #if resolution <= 1.5:
      # path=os.path.join(os.path.abspath(os.path.dirname(__file__)), '../pdb_one_halide.txt')
      # os.remove(path)
  #print(f'get {all_sites} {halide_type} sites, {all_site_deleted} sites were deleted, protein sites: {all_number_PROTEIN_sites}, DNA sites: {all_number_DNA_sites}, RNA sites: {all_number_RNA_sites}, other sites: {all_number_OTHER_sites}')
  with open ('stat.txt', 'a') as w:
           w.write(f'Halide: {halide_type}\n')
           with open(f'data/info/info_{halide_type}.txt', 'r') as f:
              text1=f.readlines()
              number_inputs=len(text1)
              base_list1=[]
              for line in text1:
                    base1=list(line.split('!'))
                    base_list1+=base1
              total_2=0
              for i in range(1,len(base_list1),4): 
                    if float(base_list1[i])<=2:
                       total_2+=1
              x_ray1=base_list1.count('X-RAY DIFFRACTION\n')
              NMR1=base_list1.count('NMR\n')
           Total_entries_after_homologous_filter=int(len(base_list)/4)
           w.write('Total_statistics\n')
           w.write(f'Total_entries: {number_inputs}\n')
           w.write(f'Entries_with_resolution_less_(and_equal)_then_2: {total_2}\n')
           w.write(f'X_RAY_entries: {x_ray1}\n')
           w.write(f'NMR_entries: {NMR1}\n')
           Other_entries=int(number_inputs)-int(x_ray1)-int(NMR1)
           w.write(f'Other_entries: {Other_entries}\n')
           w.write('Statistics_after_homologous_filter\n')
           w.write(f'Total_entries_after_homologous_filter: {Total_entries_after_homologous_filter}\n')
           w.write(f'Entries_with_resolution_less_(and_equal)_then_2: {high_resolution}\n')
           #x_ray=base_list.count('X-RAY DIFFRACTION\n')
           #NMR=base_list.count('NMR\n')
           #w.write(f'X_RAY_entries: {x_ray}\n')
           #w.write(f'NMR_entries: {NMR}\n')
           #Other_entries1=int(Total_entries_after_homologous_filter)-int(x_ray)-int(NMR)
           #w.write(f'Other_entries: {Other_entries1}\n')
           w.write(f'Total_sites: {all_sites}\n')
           w.write(f'Total_site_deleted: {all_site_deleted}\n')
           w.write(f'Protein_sites: {all_number_PROTEIN_sites}\n')
           w.write(f'DNA_sites: {all_number_DNA_sites}\n')
           w.write(f'RNA_sites: {all_number_RNA_sites}\n')
           w.write(f'Other_sites: {all_number_OTHER_sites}\n\n')
if __name__=='__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('-input_filter', type=str)
    parser.add_argument('-input_struct', type=str,
#                         help='Path to input file.')
                        help='PDB id or PDB structure in .ent format.')
    parser.add_argument('-input_type', type=str, help='Pass your input type.')
    
    # parser.add_argument('-halide', type=str, default='F', help='Type of halide')
    parser.add_argument('-angstrem_radius', type=int, default=5, help='Threshold radius in Å.')
    parser.add_argument('-output_full', type=str, 
                        help='Name of output file with full data.')
    parser.add_argument('-output_AA', type=str, 
                        help='Name of output file with AA composition ')
    # parser.add_argument('-output_dir', type=str, 
    #                 help='Name of output dir.')
    parser.add_argument('-ligands', type=str, help='y-do not ignore HETATM field in processing; n - ignore')
    parser.add_argument('-C', type=int, help='1-all atoms; 2-no C,H atoms; 3 - no C; 4 - no C=0 no C except CA')
    parser.add_argument('-prot', type=str, help='y-only protein sites; n- add DNA,RNA,ligands to result')
    parser.add_argument('-sec_struct', type=str, help='y - collect secondary scructure; n - do not collect')
    parser.add_argument('-water', type=str, help='y - collect water; n - do not collect')
    args = parser.parse_args()

    main(args)
