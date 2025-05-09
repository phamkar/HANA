



############################################################################################################################################################
#                                                     NOESY Simulation: Setting Simulated NOESY Peak List                                                  #                                                               
############################################################################################################################################################

def distance3D(atom1, atom2):
    """ 
    _________________________________________________________
    takes two coordinates. ex: ((26.266, 25.413, 2.842),
                             (26.913, 26.639, -3.51))
    returns the distance

    Euclidean calculation
    _________________________________________________________
    """
    return sqrt((atom1[0] - atom2[0]) ** 2 +
              (atom1[1] - atom2[1]) ** 2 +
              (atom1[2] - atom2[2]) ** 2)

def createDistanceMatrix(pdb_list):
    '''
    _________________________________________________________________________
    input: parsed PDB list
    output: distance matrix and key list

    c_list = 'ABCD'

    calculating the distance between each and every atom

    setting the 'keylist': chainID + '_'+ amino acid + sequence number + atom
    _________________________________________________________________________
    '''
    distMat = np.zeros( (len(pdb_list), len(pdb_list), 4, 4 ) )
    keyList = []
    for i in range(len(pdb_list)):
        nSeq, x, y, z, aaa, a, atm, c = pdb_list[i]
        keyList.append(c+'_'+a+str(nSeq)+atm)
        for j in range(i+1, len(pdb_list)):
            nSeq2, x2, y2, z2, aaa2, a2, atm2, c2 = pdb_list[j]
            dist = distance3D( (x, y, z), (x2, y2, z2))
            distMat[i, j, c_list.index(c), c_list.index(c2)] = \
            distMat[j, i, c_list.index(c2), c_list.index(c)] = dist

    return distMat, keyList

def parse_key(key):
    '''
    ____________________________________________________________________
    input: keys in KeyList
    output: breakdown of keys: chainID, amino acid, residue number, atom
    ____________________________________________________________________
    '''
    c = key.split('_')[0]
    a = key.split('_')[1][0]
    nSeq = int(re.search(regex, key.split('_')[1]).group())
    sep = c+'_'+a+str(nSeq)
    atm = key.split(sep)[1]
    return c, a, nSeq, atm

def create_noesy_peak_list(keyList,noesy_type,cs_list,distMat):
  '''
    keyList: c+'_'+a+str(nSeq)+atm
    cs_shifts: [nseq, atom, chemical shift]
    c_list = 'ABCD'

  '''
  # create header
  lines = ''
  header = '  %20s    %5s   %5s   %5s   %20s' % ('Assignments','w1','w2','w3', 'Peak Height')
  lines += header + '\n'
  lines += '\n'

  # parse keyList 
  for i in range(len(keyList)):
    key = keyList[i]
    # break down keyList; contains atoms for chain A and B
    c, a, nSeq, atm = parse_key(key) 
    for j in range (i+1,len(keyList)):
      key2 = keyList[j]     # break down keyList; contains atoms for chain A and B
      c2, a2, nSeq2, atm2 = parse_key(key2)
      dist = distMat[nSeq, nSeq2, c_list.index(c), c_list.index(c2)] # possible index combo: distMat[nSeq, nSeq2, 0, 0] --> dist == value 
                                                                     #                       distMat[nSeq, nSeq2, 0, 1] --> dist == 0's 
                                                                     #                       distMat[nSeq, nSeq2, 1, 0] --> dist == 0's
                                                                     #                       distMat[nSeq, nSeq2, 0, 0] --> dist == 0's
  return 