Notes
1.  What if user enters uncomplexedproteinpdbname, binding, ligandxyzfilename
   uncomplexedproteinpdbname=Q00534_fill.BL00020001.fixed.pdb

Minimum Input Example Binding Free Energy
complexedproteinpdbname=anilinecomp.pdb 
binding
keyfilename=aniline.key
ligandxyzfilename=aniline.xyz


Minimum Input Example Binding Free Energy
uncomplexedproteinpdbname=anilinecomp.pdb 
binding
keyfilename=aniline.key
ligandxyzfilename=aniline.xyz

from openbabel import openbabel

temp=open("complexed.Q00534_fill.BL00020001.fixed.obabel.connect.pdb",'r')
/storage1/fs1/bolton/Active/projects/BWILEYtest/Ribociclib_Bind/Q00534_fill.BL00010001_ligand_conect.fixed.out.pdb

/storage1/fs1/bolton/Active/data/hg38/test/db/KMT2A/4GQ6.pdb

from openbabel import openbabel

uncomplexedatomindices=[]
indextocoordinates={}
temp=open("/Users/brian/Bolton/Docker/openbabel/4gq6.pdb",'r')
results=temp.readlines()
temp.close()
for line in results:
   if 'ATOM' in line:
      linesplit=line.split()
      index=int(linesplit[1])
      uncomplexedatomindices.append(index) 

pdbmol=openbabel.OBMol()
obConversion = openbabel.OBConversion()
obConversion.SetInFormat('pdb')
obConversion.ReadFile(pdbmol,"/Users/brian/Bolton/Docker/openbabel/4gq6.pdb")
totalatoms=pdbmol.NumAtoms()
atomiter=openbabel.OBMolAtomIter(pdbmol)
for atom in atomiter:
   atomidx=atom.GetIdx()
   if atomidx in uncomplexedatomindices:
      coords=[atom.GetX(),atom.GetY(),atom.GetZ()]
      indextocoordinates[atomidx]=coords


indexestodelete=[]
for i in range(1,totalatoms+1):
   if i not in uncomplexedatomindices:
      indexestodelete.append(i)

indexestodelete.sort(reverse=True)

for idx in indexestodelete[0:30]:
   pdbmol.DeleteAtom(pdbmol.GetAtom(idx))
   print(pdbmol.GetAtom(682).GetResidue().GetName())

for idx in indexestodelete[30:31]:
   pdbmol.DeleteAtom(pdbmol.GetAtom(idx))
   print(pdbmol.GetAtom(848).GetResidue().GetName())


for atom in   pdbmol.GetResidue(109) pdbmol.GetAtom(848).GetType()


minimized: Q00534_fill.BL00010001_ligand_conect.fixed.out.xyz (now in side)




for idx in indexestodelete:
   atom=pdbmol.GetAtom(idx)
   pdbmol.DeleteAtom(atom)

obConversion.SetOutFormat('pdb')
obConversion.WriteFile(pdbmol,"/Users/brian/Bolton/Docker/openbabel/4gq6.obabel.pdb")


print(idx, pdbmol.GetAtom(682).GetResidue().GetName())

pdbmol.GetAtom(682).GetType()
pdbmol.GetAtom(5239).GetType()



# pdbmol.DeleteAtom(pdbmol.GetAtom(5239))
# pdbmol.GetAtom(1604).GetResidue().GetName()
# pdbmol.GetAtom(5239).GetType()
pdbmol.GetAtom(666).GetType()
pdbmol.GetAtom(666).GetResidue().GetName()
# for idx in indexestodelete[30:31]:
#    pdbmol.DeleteAtom(pdbmol.GetAtom(idx))


# for idx in indexestodelete:
#    atom=pdbmol.GetAtom(idx)
#    pdbmol.DeleteAtom(atom)

obConversion.SetOutFormat('pdb')
obConversion.WriteFile(pdbmol,"uncomplexed2.pdb")
return indextocoordi



bsub -n16 -oo gpu.log -G compute-timley -g /bwileytest -q general -M 196G -R "gpuhost rusage[mem=196GB] span[hosts=1]" -gpu "num=1:gmodel=TeslaV100_SXM2_32GB" -a "docker($POLTYPE_G)" python $POLRUN


/scratch1/fs1/bolton/brian/poltype2/PoltypeModules/poltype.py


while read -r line
do
  line=$(echo "$line" | cut -d= -f2 | sed 's/ --numproc//')
  echo 
done < "Ribociclib_Bind_proddynamicsjobs.txt"

while read -r line
do
  line=$(echo "$line" | cut -d= -f2 | sed 's/ --numproc//')
  chg_dir=$(echo "$line" | cut -d';' -f1) 
  $chg_dir
  log=$(basename $(pwd))
  com=$(echo "$line" | cut -d';' -f2)
  bsub -n4 -oo $log.log -G compute-timley -g /bwileytest -q general -M 196G -R "gpuhost rusage[mem=196GB] span[hosts=1]" -gpu "num=1:gmodel=TeslaV100_SXM2_32GB" -a "docker($POLTYPE_G)" /bin/bash -c "$com"
done < "Ribociclib_Bind_proddynamicsjobs2.txt"