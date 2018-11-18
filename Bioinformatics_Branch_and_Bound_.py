def ParentMass(spectrum):
    spec = spectrum.split(" ")
    return int(spec.pop())


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

    l.append(0)
    for i in range(len(peptide)):
        total_sum += dict.get(peptide[i:i+1])
        l.append(dict.get(peptide[i:i+1]))
    l.append(total_sum)
    for i in range(len(peptide)):
        for j in range(1, len(peptide) - 1):
            sum = 0
            if (i+j > (len(peptide) - 1)):
                tmp = peptide[i:len(peptide)] + peptide[:j - (len(peptide) - i) + 1]
            else:
                tmp = peptide[i:i+j+1]
            for k in range(len(tmp)):
                sum += dict.get(tmp[k:k+1])
            l.append(sum)
    
    l.sort()
    string += str(l[0])
    for i in range(1,len(l)):
        string += ' ' + str(l[i])
    return string

def peptideMass(str):
    dict = Amino_Acid_Mass()
    sum = 0
    for k in range(len(str)):
        sum += dict.get(str[k:k+1])
    return sum

def Expand(peptides):
    dict = Amino_Acid_Mass()
    newlist = list()
    for item in peptides:
        for key in dict:
            newSubList = list(item)
            newSubList.append(key)
            newlist.append(newSubList)
    return newlist


def consistent(peptide, spectrum):
    dict = Amino_Acid_Mass()
    #pept = peptide.split(" ")
    pept = []
    for i in range(len(peptide)):
        pept.append(peptide[i:i+1])
    for j in range(len(pept)):
        pept[j] = str(dict.get(pept[j]))
    spec = spectrum.split(" ")
    coincidence = 0
    for i in pept:
        if i in spec:
            coincidence += 1
            spec.remove(i)
    if coincidence == len(pept):
        return True
    else:
        return False


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
   peptide = ''
   peptides = ['']
   parentmass = ParentMass(spectrum)
   result = []
   while len(peptides) != 0:
           peptides = Expand(peptides)
           ln = (len(peptides))
           tmp = 0
           for i in range(ln):
           #for item in peptides:
               i -= tmp
               for j in peptides[i]:
                   peptide += j
               if peptideMass(peptide) == parentmass:
                   if CycloSpectrum(peptide) == spectrum:
                        #print(peptide)
                        result.append(peptide)
                        peptides.remove(peptides[i])
                        tmp += 1
               if not consistent(peptide, spectrum):
                   peptides.remove(peptides[i])
                   tmp += 1
               peptide = ''
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





Cyclopeptide_Scoring('NQEL','0 99 113 114 128 227 257 299 355 356 370 371 484')


#strin = CycloSpectrum('LIQKWTV')
#res = CyclopeptideSequencing(strin)
#print(peptide_to_mass(res))



