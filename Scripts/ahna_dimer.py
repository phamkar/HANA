xplor.requireVersion("2.26") # requiring Xplor NIH version 2.26
xplor.parseArguments() # check for typos on the command-line # to parse is to break down command lines and structured form that can be easily used by a program
quick = False ## quick is iterated but does not appear for the rest of the script
#=================== Load Topology and Parameters ==============================
command = xplor.command
import protocol
protocol.initRandomSeed(46191)
protocol.initTopology('protein')
protocol.initParams('protein')
#========================== Import a Sequence ==================================
import psfGen
#============================== Oligomer ======================================
num_oligomer=2
#======================== Create the PSF File ==================================
import os, sys, subprocess
#==================== Generate Symmetric Dimer Structure =======================
f = open('finalx.pdb', 'r')
lines = f.readlines()
f.close()
#print (lines)

content = ''
for line in lines:
  if len(line) < 66:
    content += line
    continue
  if line[:4] != 'ATOM':
    content += line
    continue
  content += line[:21] + ' ' + line[22:]

f = open('temp.pdb', 'w')
f.write(content)
f.close()

protocol.loadPDB('temp.pdb')

from symSimulation import SymSimulation
symSim = SymSimulation ('symSim', subSel='not PSEUDO') # atoms to duplicate

from math import pi
from vec3 import Vec3
from mat3 import rotVector
symSim.addCopy(rotVector(Vec3(0,0,1),pi), Vec3(10,0,0), segidSuffix='B')
symSim.symOp(0).segidSuffix='A'

import regularize
regularize.addUnknownAtoms()
#====================== Output PSF and PDB Files ===============================
#xplor.command("write psf output=thread.psf end")
if num_oligomer >= 2:
  os.system("rm -rf duplicate.pdb*")
  protocol.writePDB("duplicate.pdb", selection=AtomSel('all',symSim))
  xplor.command('delete selection=(all) end')
  #protocol.writePDB("duplicate.pdb", selection=AtomSel('all',symSim))
  protocol.loadPDB('duplicate.pdb')
else:
   subprocess.run(['python3', '/Users/Karen/ahna/run/BestEvaluated/2n74_demo/ahna.py'])
#============================= Annealing from here ============================#
#==================== Initialize List of Potentials ===========================#
from potList import PotList
potList = PotList()
crossTerms=PotList('cross terms') # can add some pot terms which are not
                                  # refined against- but included in analysis
#=============== Initialize Simulated Annealing Parameters=====================#
from simulationTools import MultRamp, StaticRamp
rampedParams=[]
highTempParams=[]

#======================= Residual Dipolar Coupling ===========================#
# moved upto here to make this work with CDIH
# make a copy A and B
iNumOfRdcOrientations = 0
segids = 'ABCD'
if os.path.exists("rdc_ponderosa.cfg"):
  from varTensorTools import create_VarTensor
  f = open( 'rdc_ponderosa.cfg' , 'r')
  rdclines = f.readlines()
  f.close()

  iNumOfRdcOrientations = int(rdclines[0].strip().split()[1]) # NUM_OF_ORIENTATIONS 2

# check the RDC orientation sequence number used: 360
if iNumOfRdcOrientations > 0:
    rdc_tensor = {}
    rdc_weight = {}
    tensorCalc_list = []
    for i in range(1, iNumOfRdcOrientations+1):
      rdcline = rdclines[i]
      print(rdcline) #output '1 999 999 360'
      rdcline_splited = rdclines[i].strip().split()

      for j in range(1, num_segid+1): # j is 1, 2
        rdc_ori_id = rdcline_splited[0] # output: '1'
        ##oTensor = create_VarTensor(rdc_ori_id)
        oTensor = create_VarTensor(rdc_ori_id+'_'+segids[j-1]) # 1_A, 1_B, 2_A, 2_B, 3_A, 3_B - create twice
        try:
            rdclines3 = ''.join(rdcline_splited[3]) # rdcline_splited[3] from list to str
                                                      # output: '360'
            fDa = float(rdcline_splited[1]) # 999
            fRh = float(rdcline_splited[2]) # 999
            oTensor.setDa(fDa)
            oTensor.setRh(fRh)
            #oTensor.setFreedom ( 'varyDa , varyRh ') # allow Da, Rh to vary
            ##rdc_tensor[rdc_ori_id] = oTensor
            rdc_tensor[rdc_ori_id+'_'+segids[j-1]] = oTensor #rdc_tensor['1_A'] = oTensor
            #rdc_weight[rdc_ori_id+'_'+segids[j-1]] = float(rdcline_splited[3]) # output: rdc_weight['1_A'] = '360'
            rdc_weight[rdc_ori_id+'_'+segids[j-1]] = int(rdclines3) + int(j-1) #output 360 361...
            ##rdc_weight[rdc_ori_id] = float(rdclines3) + float(j-1)
            tensorCalc_list.append(0)
            print (rdc_weight)
        except:
            oTensor.setDa(10.0) # arbitrary
            oTensor.setRh(0.5) # arbitrary
            #oTensor.setFreedom ( 'varyDa , varyRh ') # allow Da, Rh to vary
            ##rdc_tensor[rdc_ori_id] = oTensor
            rdc_tensor[rdc_ori_id+'_'+segids[j-1]] = oTensor
            rdc_weight[rdc_ori_id+'_'+segids[j-1]] = int(rdclines3) + int(j-1)  #output: {'1_A': 360.0, '1_B': 361.0}
            ##rdc_weight[rdc_ori_id] = float(rdclines3) + float(j-1)
            tensorCalc_list.append(1)

    # files
    # iNumOfRdcTables -> 1
    iNumOfRdcTables = int(rdclines[iNumOfRdcOrientations+1].strip().split()[1])  # NUM_OF_RDC_FILES 1
    rdc_table = []
    # create two rdc files

    # double this part
    for i in range(len(rdclines)):
      if i < iNumOfRdcOrientations+2: continue
      rdcline = rdclines[i]
      rdcline_splited = rdclines[i].strip().split() # ['rdc_ponderosa_1.tbl', '1', '360', '1', '1.000000',],

      # make rdc_ponderosa_1_A.tbl, rdc_ponderosa_1_B.tbl
      f = open(rdcline_splited[0], 'r')
      content = f.readlines()
      f.close()


      # write rdc_ponderosa_1_1.tbl, rdc_ponderosa_1_2.tbl
      for j in range(1, num_segid+1): # j is 1 2 ...

        new_content = ''


        # rdcline_splited[0][:-4] <- rdc_ponderosa_1
        new_rdctable_file = rdcline_splited[0][:-4] + '_' + segids[j-1] + '.tbl'

        f = open(new_rdctable_file, 'w') #rdc_ponderosa_1_A.tbl, rdc_ponderosa_1_B.tbl
        for line in content:
          temp = line.replace ('resid  ' + rdcline_splited[2],
                               'resid  ' + str(int(rdcline_splited[2]) + j-1) )
          new_content += temp.replace('(', '(segid %s and ' % (segids[j-1]))
        f.write(new_content)
        f.close()
        new_rdcline_splited = list(rdcline_splited) # output: [rdc_ponderosa_1_B.tbl', '2', '360', '1', '1.000000']
        new_rdcline_splited[0] = new_rdctable_file
        new_rdcline_splited[1] = str(j)
        new_rdcline_splited[2] = int(rdcline_splited[2]) + (j-1)# 360...
        new_rdcline_splited = new_rdcline_splited + [['A', 'B', 'C', 'D'][j-1],]
        rdc_table.append(new_rdcline_splited) # output: ['rdc_ponderosa_1_1.tbl', '1', '360', '1', '1.000000','A'],


    # up to here
    # Double the RDC orientations
    iNumOfRdcOrientations = 2 * iNumOfRdcOrientations


    from rdcPotTools import create_RDCPot, scale_toNH
    rdcs = PotList('rdc')
    for (rdc_file, rdc_orient, rdc_ori, rdc_type, rdc_weightStr,segids) in rdc_table:
        # ['rdc_ponderosa_1_1.tbl', '1', '360', '1', '1.000000','A']
      rdc_weight = float(rdc_weightStr)
    #rdc = create_RDCPot("%s_%s"%(rdc_orient,rdc_type),
    #              rdc_file,
    #              rdc_tensor[rdc_ori_id])
      rdc = create_RDCPot("%s_%s"%(rdc_orient,rdc_type),
                          rdc_file,
                          rdc_tensor[rdc_ori_id+'_'+segids]) # "1_1"
      rdc.setScale(rdc_weight) # 1.0
      rdc.setShowAllRestraints(1) #all restraints are printed during analysis
      rdc.setThreshold(2)       # in Hz
      rdcs.append(rdc)
      pass
    potList.append(rdcs)
    #hipotList.append(rdcs)
    #rampedParams.append( MultRamp(  0.05, 0.25, "rdcs.setScale( VALUE )") )


    rampedParams.append( MultRamp( 0.01, 2.00, "rdcs.setScale( VALUE )") )
                                 #0.01,1.00,

    #highTempParams.append( StaticRamp("rdcs.setScale( 0.05 )"))
    # calc. initial tensor orientation
    from varTensorTools import calcTensorOrientation, calcTensor
#  for i in range(len(rdc_tensor.keys())):
#    medium = rdc_tensor.keys()[i]
#    if tensorCalc_list[i] == 1:
#      calcTensorOrientation(rdc_tensor[medium])
#      rampedParams.append( StaticRamp("calcTensor(rdc_tensor['%s'])" % medium) )
#      print('Xplor-NIH calculated tensor will be used.')
#    pass
  #except Exception as e:
   # print(e)
    #pass
    rdc_key_list = list(rdc_tensor.keys())
    for i in range(len(rdc_key_list)):
      medium = rdc_key_list[i]
      if tensorCalc_list[i] == 1:
        calcTensorOrientation(rdc_tensor[medium])

        rampedParams.append( StaticRamp("calcTensor(rdc_tensor['%s'])" % medium) )
        print('Xplor-NIH calculated tensor will be used.')
      pass
#    except Exception as e:
#       print(e)
#    pass

#======================= Pseudo Contact Shift ===========================#
# moved upto here to make this work with CDIH
# make a copy A and B
iNumOfPcsOrientations = 0
#if os.path.exists("pcs_ponderosa.cfg"):
if os.path.exists('pcs_ponderosa.cfg'):
  from varTensorTools import create_VarTensor
  f = open( 'pcs_ponderosa.cfg' , 'r')
  pcslines = f.readlines()
  f.close()


  # orientations -> need to create tensor
  try:
    iNumOfPcsOrientations = int(pcslines[0].strip().split()[1]) # NUM_OF_ORIENTATIONS 2
    if iNumOfPcsOrientations == 0:
        pass
    else:
        pcs_tensor = {}
        pcs_weight = {}
        tensorPcsCalc_list = []
        for i in range(1, iNumOfPcsOrientations+1):
          pcsline = pcslines[i]
          pcsline_splited = pcslines[i].strip().split()
          oTensor = create_VarTensor(pcsline_splited[0])
          try:
              fDa = float(pcsline_splited[1])
              fRh = float(pcsline_splited[2])
              oTensor.setDa(fDa)
              oTensor.setRh(fRh)
              oTensor.setFreedom ( 'varyDa , varyRh ') # allow Da, Rh to vary
              pcs_tensor[pcsline_splited[0]] = oTensor
              pcs_weight[pcsline_splited[0]] = float(pcsline_splited[3])
              tensorPcsCalc_list.append(0)
          except:
              oTensor.setDa(10.0) # arbitrary
              oTensor.setRh(0.5) # arbitrary
              #oTensor.setFreedom ( 'varyDa , varyRh ') # allow Da, Rh to vary
              pcs_tensor[pcsline_splited[0]] = oTensor
              pcs_weight[pcsline_splited[0]] = float(pcsline_splited[3])
              tensorPcsCalc_list.append(1)
          pass
        # files
        iNumOfPcsTables = int(pcslines[iNumOfPcsOrientations+1].strip().split()[1])  # NUM_OF_PCS_TABLES 10
        pcs_table = []
        for i in range(len(pcslines)):
          if i < iNumOfPcsOrientations+2: continue
          pcsline = pcslines[i]
          pcsline_splited = pcslines[i].strip().split()
          pcs_table.append(pcsline_splited)

        from rdcPotTools import create_RDCPot, scale_toNH
        pcss = PotList('pcs')

        for (pcs_file, pcs_orient, pcs_ori, pcs_type, pcs_weightStr) in pcs_table:
          pcs_weight = float(pcs_weightStr)
          pcs = create_RDCPot("%s_%s"%(pcs_orient,pcs_type),pcs_file,pcs_tensor[pcs_orient])
          pcs.setScale(pcs_weight)
          pcs.setShowAllRestraints(1) #all restraints are printed during analysis
          pcs.setThreshold(2)       # in Hz
          pcs.setUseDistance(True)
          pcss.append(rdc)
          pass
        potList.append(pcss)
        #hipotList.append(rdcs)
        #rampedParams.append( MultRamp( 0.05, 0.25, "rdcs.setScale( VALUE )") )
        rampedParams.append( MultRamp( 0.01, 1.00, "pcss.setScale( VALUE )") )
        #highTempParams.append( StaticRamp("rdcs.setScale( 0.05 )"))

        # calc. initial tensor orientation
        from varTensorTools import calcTensorOrientation, calcTensor
        for i in range(len(pcs_tensor.keys())):
          medium = pcs_tensor.keys()[i]
          if tensorPcsCalc_list[i] == 1:
            calcTensorOrientation(pcs_tensor[medium])
            rampedParams.append( StaticRamp("calcTensor(pcs_tensor['%s'])" % medium) )
          pass
  except:
    pass
#===================== Distance Restraint Potentials ==========================#
from noePotTools import create_NOEPot

if os.path.exists('noe.tbl'):

    # Copy the content of source to destination
    f = open("noe.tbl", 'r')
    lines = f.readlines()
    f.close()
    f = open("noe_dimer.tbl", 'w')
    contentA, contentB = '', ''
    for line in lines:
        # add segid A and B
        new_line_with_segidA = line.replace('( ', '( segid A and ')
        new_line_with_segidB = line.replace('( ', '( segid B and ')

        contentA += new_line_with_segidA
        contentB += new_line_with_segidB

    f.write(contentA)
    f.write(contentB)
    f.close()

    noe = create_NOEPot('noe', 'noe_dimer.tbl')
    noe.setPotType("soft")
    potList.append(noe)

    #f = open('scale.txt', 'r')
    #data_list = f.readlines()
    #f.close()
    #data = float((''.join(data_list)))[1:10]
    rampedParams.append (MultRamp(1.0, 15.0  , "noe.setScale( VALUE )"))
                                #1.0, 15.0

if os.path.exists('alt.tbl'):
# alt.tbl are the distant restraints
    alt = create_NOEPot('alt', 'alt.tbl')
    noe.setPotType("soft")
    potList.append(alt)
    rampedParams.append (MultRamp( 0.5, 30.0,  "alt.setScale( VALUE )"))
                                #0.5, 20.0
if os.path.exists('dist_ca.tbl'):

  f = open("dist_ca.tbl", 'r')
  lines = f.readlines()
  f.close()
  f = open("dist_ca_dimer.tbl", 'w')
  contentA, contentB = '',''
  for line in lines:

    i = line.split()
    #max number
    upperLim = (i[15])
    #print (upperLim)
    #reference number
    referenceNum = (i[16])
    #print (referenceNum)

    lowerLim = float(upperLim) - float(referenceNum)
    averageLim = (float(upperLim)- lowerLim)/2.00

    newUpperA = float(upperLim) - averageLim
    newUpper = round(newUpperA, 2)
    newReference = float(referenceNum)/2.0
    #print (newUpper)
    i[15] = str(newUpper)
    i[16] = str(newReference)

    i_modified = ' '.join(i)
    # add segid A and B
    new_line_with_segidA = i_modified.replace('( ', '( segid A and ')
    new_line_with_segidB = i_modified.replace('( ', '( segid B and ')

    contentA += new_line_with_segidA
    contentB += new_line_with_segidB


  f.write(contentA)
  f.write(contentB)
  f.close()

  distca = create_NOEPot('distca', 'dist_ca_dimer.tbl')
  distca.setPotType("soft")
  potList.append(distca)

  rampedParams.append(MultRamp( 0.25, 5.0, "distca.setScale( VALUE )"))
#========================= Dihedral Potential =================================#
# Set up dihedral angles
from xplorPot import XplorPot
    # make a copy for A and B

if os.path.exists("aco.tbl"):

    f = open("aco.tbl", 'r')
    lines = f.readlines()
    f.close()
    f = open("aco_dimer.tbl", 'w')

    segidA, segidB = '', ''

    for line in lines:
      segidA += line.replace('(', '(segid A and ')
      segidB += line.replace('(', '(segid B and ')

    f.write(segidA)
    f.write(segidB)
    f.close()
    protocol.initDihedrals("aco_dimer.tbl")
    potList.append( XplorPot('CDIH') )
    potList['CDIH'].setThreshold( 5 )
    rampedParams.append( StaticRamp("potList['CDIH'].setScale(200)") )

#========================= H-BOND Potential ===================================#
if os.path.exists("hbond.tbl"):
    f = open("hbond.tbl", 'r')
    lines = f.readlines()
    f.close()

    f = open("hbond_dimer.tbl", 'w')

    contentA, contentB = '', ''
    for line in lines:
        # add segid A and B
        new_line_with_segidA = line.replace('( ', '( segid A and ')
        new_line_with_segidB = line.replace('( ', '( segid B and ')

        #
        contentA += new_line_with_segidA
        contentB += new_line_with_segidB

    f.write(contentA)
    f.write(contentB)
    f.close()
    hbond = create_NOEPot('hbond','hbond_dimer.tbl')
    hbond.setPotType("soft")
    potList.append(hbond)
    rampedParams.append(StaticRamp("potList['hbond'].setScale(30)"))
if os.path.exists("hbda.tbl"):

    # Copy the content of source to destination
    f = open("hbda.tbl", 'r')
    lines = f.readlines()
    f.close()

    f = open("hbda_dimer.tbl", 'w')

    contentA, contentB = '', ''
    for line in lines:
        # add segid A and B
        new_line_with_segidA = line.replace('( ', '( segid A and ')
        new_line_with_segidB = line.replace('( ', '( segid B and ')


        contentA += new_line_with_segidA
        contentB += new_line_with_segidB

    f.write(contentA)
    f.write(contentB)
    f.close()
    protocol.initHBDA("hbda_dimer.tbl")
    potList.append( XplorPot('HBDA') )

#protocol.initHBDB()
#potList.append( XplorPot('HBDB') )
#========================= Torsion Potential ==================================#
from torsionDBPotTools import create_TorsionDBPot
torsionDBPot = create_TorsionDBPot('tDB')
potList.append( torsionDBPot )
rampedParams.append( MultRamp( 0.002, 4, "torsionDBPot.setScale(VALUE)") )
                            #0.002, 2
#=============== Bonds, Angles, and Impropers Potentials ======================#
#covalent terms
from xplorPot import XplorPot
for term in ('BOND', 'ANGL', 'IMPR'):
    potList.append( XplorPot(term) )
    pass
# Set threshold for terms in potList to allow violation analysis.

potList['ANGL'].setThreshold(5.0)   # default is 2.0
potList['IMPR'].setThreshold(5.0)   # default is 2.0
# Use default values for the rest (bond: 0.05, cdih: 5.0, noe: 0.5).

rampedParams.append( MultRamp( 0.4, 2.0, "potList['ANGL'].setScale(VALUE)"))
                            #0.4, 1.0
rampedParams.append( MultRamp( 0.1, 2.0, "potList['IMPR'].setScale(VALUE)"))
                            #0.4, 1.0
#=================================== EEFX =====================================#
#=========================== NOT USED FOR FOLDING =============================#
#from eefxPotTools import create_EEFxPot, param_LK
#eefxpot=create_EEFxPot("eefxpot","not name H*",paramSet=param_LK)
#eefxpot.setVerbose(False)
#potList.append(eefxpot)
#rampedParams.append(MultRamp( 0.1,1.0,"eefxpot.setScale(VALUE)"))
#rampedParams.append(StaticRamp("potList['VDW'].setScale(0)"))
#use repel at high temp
#potList.append( XplorPot("VDW") )
#highTempParams.append(StaticRamp("""protocol.initNBond(cutnb=100,rcon=0.04,
#                                                     tolerance=45,nbxmod=5,
#                                                     repel=1.2,onlyCA=True)"""))
#highTempParams.append(StaticRamp("potList['VDW'].setScale(1)"))
#highTempParams.append(StaticRamp("eefxpot.setScale(0)"))
#====================== Van Der Waal's Potential ==============================#
potList.append( XplorPot('VDW') )
#hipotList.append( XplorPot('VDW') )
if xplor.version_info[0] * 10000 + xplor.version_info[1] < 20043:
    rampedParams.append( StaticRamp("protocol.initNBond()") )
    rampedParams.append( MultRamp( 0.9, 0.8,
                              "command('param nbonds repel VALUE end end')") )
else:
    rampedParams.append( StaticRamp("protocol.initNBond(repel=0.9)") )
rampedParams.append( MultRamp( .004, 4,
                              "command('param nbonds rcon VALUE end end')") )
## Nonbonded interaction occurs only between CA atoms at high temperatures
highTempParams.append( StaticRamp("""protocol.initNBond(cutnb=100,
                                                        rcon=0.004,
                                                        tolerance=45,
                                                        repel=1.2,
                                                        nbxmod=5,
                                                        onlyCA=1)""") )
# ==========================Gyration Volume Term==========================#
if os.path.exists(".use_gyr"):
  from gyrPotTools import create_GyrPot
  gyr = create_GyrPot("Vgyr","resid 3:15")
  potList.append(gyr)
  #doubled
  rampedParams.append( MultRamp( .004, 4, "gyr.setScale(VALUE)") )

#===================== Small Angle X-ray Scattering ===========================#
if os.path.exists("saxs_ponderosa.tbl") and os.path.exists(".use_saxs"):
  import solnXRayPotTools
  from solnXRayPotTools import create_solnXRayPot
  xray = create_solnXRayPot('xray',experiment='saxs_ponderosa.tbl',
                      aSelection='not PSEUDO and not name H* and (resid 3:8 or resid 15)',
                      numPoints=50,
                      useInternalSpline=True,
                      normalizeIndex=-3,preweighted=False)

  xrayCorrect = create_solnXRayPot('xray-c',experiment='saxs_ponderosa.tbl',
                      aSelection='not PSEUDO and not name H* and (resid 3:8 or resid 15)',
                      numPoints=50,
                      useInternalSpline=True,
                      normalizeIndex=-3,preweighted=False)

  solnXRayPotTools.useGlobs(xray)

  xray.setNumAngles(50)
  xrayCorrect.setNumAngles(500)
  xray.setScale(40)
  xray.setCmpType("plain")

  potList.append(xray)
  #hipotList.append(xray)
  crossTerms.append(xrayCorrect)

  print(xray.calcEnergy())
  from solnScatPotTools import fitParams
  rampedParams.append( StaticRamp("fitParams(xrayCorrect);xray.calcGlobCorrect(xrayCorrect.calcd())", stride=10))
  #rampedParams.append( StaticRamp("fitParams(xrayCorrect);xray.calcGlobCorrect(xrayCorrect.calcd())", stride=10))
  #rampedParams.append( StaticRamp("xray.setScale(5000)") )
  #highTempParams.append( StaticRamp("fitParams(xrayCorrect);xray.calcGlobCorrect(xrayCorrect.calcd())", stride=10))
  #highTempParams.append( StaticRamp("xray.setScale(2000)") )

#============================== Density Map ===================================#
try:
  if os.path.exists("_em.map") and os.path.exists(".use_em"):
    from densityGrid import DensityGrid
    dmap = DensityGrid()
    dmap.readCCP4("_em.map",verbose=True)

    from probDistPotTools import create_probDistPot
    mapProb = create_probDistPot("mapProb", dmap, potType="cross_correlation",
                                            sel="not PSEUDO and not name H* and (resid 3:8 or resid 15)")
#                                                    sel="not pseudo and not name H*")
    #mapProb.setScale(10)  # use default 1
    #rampedParams.append( MultRamp( 0.1, 10.0, "mapProb.setScale(VALUE)") )
    potList.append(mapProb)
except:
  print("WARNING: DensityGrid could not be applied properly.")
  pass
#========================== Initial Minimization ==============================#
# Give atoms uniform weights, except for the anisotropy axis.
protocol.massSetup()

#### Set up IVM object(s).
#### IVM object for torsion-angle dynamics.
from ivm import IVM
dyn = IVM()
### Initially minimize the structure in cartesian space with just the covalent
### constraints, because torsion dynamics can't change bond lengths, bond
### angles, and many improper angles.
protocol.cartesianTopology(dyn)
protocol.initMinimize(dyn,
                      potList=[XplorPot(name)
                               for name in ('BOND','ANGL','IMPR')],
                      numSteps=1000)
dyn.run()

### Resets topology for torsion angle dynamics
dyn.reset()
protocol.torsionTopology(dyn)

#### IVM object for final Cartesian minimization.
minc = IVM()
protocol.cartesianTopology(minc)

#### The object that performs simulated annealing
from simulationTools import AnnealIVM
init_t  = 3500
final_t = 25
cool = AnnealIVM(initTemp = init_t,
                 finalTemp= final_t,
                 tempStep = 12.5,
                 ivm=dyn,
                 rampedParams = rampedParams)



#======================= Structure Loop ==============================#
### Generate a new structure with randomized torsion angles
def calcOneStructure(loopInfo):
    ### Generate a new structure with randomized torsion angles
    from monteCarlo import randomizeTorsions
    randomizeTorsions(dyn)
    protocol.fixupCovalentGeom(maxIters=100,useVDW=1)

    ### High Temperature Dynamics Stage.
    ### Set torsion angles from restraints.
    if os.path.exists("aco_dimer.tbl"):
        from torsionTools import setTorsionsFromTable
        setTorsionsFromTable("aco_dimer.tbl")

    ### Initialize parameters for high temperature dynamics.
    from simulationTools import InitialParams
    InitialParams( rampedParams )
    InitialParams( highTempParams )

    ### Torsion angle minimization.
    #protocol.initMinimize(dyn,potList=potList,
    #                      numSteps=50,printInterval=10)
    #dyn.run()

    protocol.writePDB("step1.pdb")
    ### High temperature dynamics
    protocol.initDynamics(dyn,potList=potList,bathTemp=init_t,
                          initVelocities=1,finalTime=1000,numSteps=10000,printInterval=2000)
    dyn.setETolerance(init_t/100) # used to det. stepsize. default: temp/1000
    dyn.run()
    protocol.writePDB("step2.pdb")
    ### Initialize parameters for simulated annealing
    InitialParams( rampedParams )

    ### Initializes integrator for simulated annealing
    protocol.initDynamics(dyn,potList=potList,
                  finalTime=10,numSteps=100,printInterval=100)

    ### Start simulated annealing
    cool.run()
    protocol.writePDB("step3.pdb")

    ### The final torsion angle minimization
    protocol.initMinimize(dyn,printInterval=400)
    dyn.run()

    ### The final cartesian coordinate minimization
    protocol.initMinimize(minc,potList=potList,dEPred=10)
    minc.run()
    protocol.writePDB("step4.pdb")

    os.system("rm -rf step*.pdb_*")
    pass

#====================== Loop Control and Output ===============================#
from simulationTools import StructureLoop, FinalParams
try:
    StructureLoop(numStructures=2,
              pdbTemplate="ahna_dimer_STRUCTURE.pdb",
              ### Executes the structure calculation
              structLoopAction=calcOneStructure,
              ### Generate structure files
              doWriteStructures=True,
              ### Generate violation summary file
              genViolationStats=True,
              ### Potentials to use for calculating average structure
              averagePotList=potList,
              ### The fraction of lowest energy structures used for average
              averageTopFraction=1,
              ### Minimizes average structure
              averageContext=FinalParams(rampedParams),
              ### Generate a regularized average structure
              averageFilename="ahna_ave.pdb",
              ### Structures are aligned by the selected atoms
              #averageFitSel="name CA or name C or name N or name O",
              averageFitSel="resid 4:14",
              ### RMSDs are calculated with all non-heavy atoms
              #averageCompSel="not name H* and not PSEUDO"
              averageCompSel="resid 4:14"
              ).run()


except Exception as e:
  print(e)
#os.system("rm -rf ahna_*.pdb_*")