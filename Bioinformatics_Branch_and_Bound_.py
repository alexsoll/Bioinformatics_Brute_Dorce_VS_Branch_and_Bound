import time

def ParentMass(spectrum):
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
    l.append(0)
    for amino_acid in peptide:
        elem = dict.get(amino_acid)
        total_sum += elem
        l.append(elem)
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
            l.append(sum)
    
    l.sort()
    return l


def peptideMass(pept):
    dict = Amino_Acid_Mass()
    sum = 0
    for item in pept:
        sum += dict.get(item)
    return sum

def Linear_Spectrum(peptide):
    spectrum = [0]
    for i, item in enumerate(peptide):
        for window in range(len(peptide) - i):
            spectrum.append(peptideMass(peptide[i:i + window + 1]))
    return spectrum

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
    return newlist


def consistent(peptide, spectrum):
    dict = Amino_Acid_Mass()
    pept = Linear_Spectrum(peptide)
    spec = list(spectrum)
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
           ln = (len(peptides))
           tmp = 0
           for i in range(ln):
               i -= tmp
               if peptideMass(peptides[i]) == parentmass:
                   if CycloSpectrum(peptides[i]) == spectrum:
                        result.append(peptides[i])
                        peptides.remove(peptides[i])
                        tmp += 1
               elif not consistent(peptides[i], spectrum):
                   peptides.remove(peptides[i])
                   tmp += 1
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

spectrum = '0 113 128 186 241 299 314 427'.split(" ")
spectrum = list(map(int,spectrum))
#spectrum = CycloSpectrum('LIWTVEG')
spectrum.sort()
res = CyclopeptideSequencing(spectrum)
print(peptide_to_mass(res))

#consistent('IIEGTWW', spectrum)
#print(CycloSpectrum('IIEGTWW') == spectrum)

print("--- %s seconds ---" % (time.time() - start_time))

