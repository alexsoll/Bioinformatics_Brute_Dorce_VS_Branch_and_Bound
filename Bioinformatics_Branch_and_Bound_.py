import time
from collections import Counter

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


#####################################################################


def Score(peptide, spectrum):
    theory_spec = Linear_Spectrum(peptide)
    experimental_spectrum = list(spectrum)
    score = 0
    for item in experimental_spectrum:
        if item in theory_spec:
            theory_spec.remove(item)
            score += 1
    return score

def Trim(LeaderBoard,spectrum, N):
    if len(LeaderBoard) < N:
        return LeaderBoard
    score_list = []
    for peptide in LeaderBoard:
        score = Score(peptide,spectrum)
        score_list.append(score)
    score_list, LeaderBoard = (list(t) for t in zip(*sorted(zip(score_list, LeaderBoard), reverse = True)))
    new_leaderBoard = LeaderBoard[:N]
    for i in range(N,len(LeaderBoard)):
        if (score_list[i] == score_list[N-1]):
            new_leaderBoard.append(LeaderBoard[i])
        else:
            break
    return new_leaderBoard 

def LeaderBoardCyclopeptideSequencing(spectrum, N):
    leaderBoard = ['']
    LeaderPeptide = []
    parentmass = ParentMass(spectrum)
    while len(leaderBoard) != 0:
        leaderBoard = Expand(leaderBoard,spectrum)
        ln = len(leaderBoard)
        tmp = 0
        for i in range(ln):
            i -= tmp
            if peptideMass(leaderBoard[i]) == parentmass:
                if Score(leaderBoard[i], spectrum) > Score(LeaderPeptide, spectrum):
                    LeaderPeptide = list(leaderBoard[i])
            elif peptideMass(leaderBoard[i]) > parentmass:
                leaderBoard.remove(leaderBoard[i])
                tmp += 1
        leaderBoard = Trim(leaderBoard,spectrum, N)
    return LeaderPeptide

def print_mass(peptide):
    dict = Amino_Acid_Mass()
    result = ''
    for amino_acid in peptide:
        result += str(dict.get(amino_acid)) + '-'
    result = result[:len(result)-1]
    return result

start_time = time.time()

#Cyclopeptide_Scoring('NQEL','0 99 113 114 128 227 257 299 355 356 370 371 484')

#spectrum = '0 113 128 186 241 299 314 427'.split(" ")
#spectrum = list(map(int,spectrum))
#spectrum = CycloSpectrum('LIWTVEG')
#spectrum.sort()
#res = CyclopeptideSequencing(spectrum)
#print(peptide_to_mass(res))

spectrum = '0 71 101 103 113 114 128 131 156 156 172 199 232 242 259 269 270 287 300 303 313 372 372 373 388 398 400 414 431 459 469 486 501 501 503 528 545 570 572 572 587 604 614 642 659 673 675 685 700 701 701 760 770 773 786 803 804 814 831 841 857 874 901 917 917 942 945 959 960 970 972 1002 1073'.split(" ")


N = 10
#spectrum = '0 71 113 129 147 200 218 260 313 331 347 389 460'.split(" ")
spectrum = list(map(int, spectrum))
result = LeaderBoardCyclopeptideSequencing(spectrum, N)
print(print_mass(result))

#consistent('IIEGTWW', spectrum)
#print(CycloSpectrum('IIEGTWW') == spectrum)

print("--- %s seconds ---" % (time.time() - start_time))