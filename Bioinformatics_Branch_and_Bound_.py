import time

def ParentMass(spectrum):
    #spec = spectrum.split(" ")
    return int(max(spectrum))


def Amino_Acid_Mass():
    dict = {}
    dict['A'] = 71; dict['I'] = 113; dict['N'] = 114; dict['D'] = 115; dict['C'] = 103; dict['Q'] = 128
    dict['E'] = 129; dict['G'] = 57; dict['H'] = 137; dict['L'] = 113; dict['K'] = 128; dict['M'] = 131
    dict['F'] = 147; dict['P'] = 97; dict['S'] = 87; dict['T'] = 101; dict['W'] = 186; dict['Y'] = 163
    dict['V'] = 99; dict['R'] = 156
    return dict

def CycloSpectrum(peptide):
    l = []
    sum = 0
    total_sum = 0
    string = ''
    dict = Amino_Acid_Mass()
    #l.append(str(0))
    l.append(0)
    for amino_acid in peptide:
        elem = dict.get(amino_acid)
        total_sum += elem
        #l.append(str(elem))
        l.append(elem)
    #l.append(str(total_sum))
    l.append(total_sum)
    ln = len(peptide)
    for i, val in enumerate(peptide):
        for j in range(1, ln - 1):
            sum = 0
            if (i+j > (ln - 1)):
                tmp = peptide[i:ln] + peptide[:j - (ln - i) + 1]
            else:
                tmp = peptide[i:i+j+1]
            for char in tmp:
                sum += dict.get(char)
            #l.append(str(sum))
            l.append(sum)
    
    l.sort()
    return l
    #string += str(l[0])
    #for i in range(1,len(l)):
    #    string += ' ' + str(l[i])
    #return string

def peptideMass(pept):
    dict = Amino_Acid_Mass()
    sum = 0
    for item in pept:
        sum += dict.get(item)
    return sum

def Expand(peptides,spectrum):
    dict = Amino_Acid_Mass()
    newlist = list()
    for item in peptides:
        for key in dict:
            if not consistent(list(key), spectrum):
                continue
            newSubList = list(item)
            newSubList.append(key)
            newlist.append(newSubList)
            #newlist.append(item+key)
    return newlist


def consistent(peptide, spectrum):
    dict = Amino_Acid_Mass()
    #pept = []
    #for i in range(len(peptide)):
    #    pept.append(peptide[i:i+1])
    #for j in range(len(pept)):
    #    pept[j] = str(dict.get(pept[j]))
    #spec = spectrum.split(" ")
    #coincidence = 0
    pept = list(peptide)
    spec = list(spectrum)
    for i,val in enumerate(pept):
        #pept[i] = str(dict.get(val))
        pept[i] = dict.get(val)

    for item in pept:
        if item in spec:
            spec.remove(item)
        else:
            return False
    return True


def peptide_to_mass(peptide):
    dict = Amino_Acid_Mass()
    res = ''
    tmp = ''
    st = set()
    for item in peptide:
        for char in item:
            tmp += str(dict.get(char)) + '-'
        tmp = tmp[:len(tmp)-1]
        st.add(tmp)
        tmp = ''
    for item in st:
        res += item + ' '
    return res


def CyclopeptideSequencing(spectrum):
   #peptide = ''
   peptides = ['']
   parentmass = ParentMass(spectrum)
   result = []
   while len(peptides) != 0:
           peptides = Expand(peptides,spectrum)
           #ln = (len(peptides))
           #tmp = 0
           #to_delete = []
           #for peptide in peptides:
           ln = (len(peptides))
           tmp = 0
           for i in range(ln):
           #for item in peptides:
               i -= tmp
               #for j in peptides[i]:
               #   peptide += j
               if peptideMass(peptides[i]) == parentmass:
                   if CycloSpectrum(peptides[i]) == spectrum:
                        #print(peptide)
                        result.append(peptides[i])
                        peptides.remove(peptides[i])
                        tmp += 1
                        #to_delete.append(peptide[i])
               elif not consistent(peptides[i], spectrum):
                   #to_delete.append(peptides[i])
                   peptides.remove(peptides[i])
                   tmp += 1
               #peptide = ''
           #for item in to_delete:
            #   peptides.remove(item)
   return result

def Cyclopeptide_Scoring(peptide, spectrum):
    theory_spec = CycloSpectrum(peptide)
    theory_spec = theory_spec.split(" ")
    spectrum = spectrum.split(" ")
    score = 0
    for item in spectrum:
        if item in theory_spec:
            theory_spec.remove(item)
            score += 1
    print(score)


start_time = time.time()

#Cyclopeptide_Scoring('NQEL','0 99 113 114 128 227 257 299 355 356 370 371 484')

#spectrum = '0 113 128 186 241 299 314 427'.split(" ")
#spectrum = list(map(int,spectrum))
spectrum = CycloSpectrum('LIWTVEG')
spectrum.sort()
res = CyclopeptideSequencing(spectrum)
print(peptide_to_mass(res))

print("--- %s seconds ---" % (time.time() - start_time))

