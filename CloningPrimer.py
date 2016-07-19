#input sequences
gene_name ='ZF974'
forward_sequence = 'gagcgccccttccagtgtcgcatttgcatgcggaactttt'
reverse_sequence = 'aaaaactggtttgagagttcacctaaaaacccacctgaga'
restriction1 = 'CTCGAGAAAAAAATG' #XhoI-Kozak-ATG
restriction2 = 'ggagatggtgctggattgattgatATCGAT' #linker1-ClaI
flanking1 = 'ATATTA'
flanking2 = 'ATATTA'

#initialization
import numpy as np
from selenium import webdriver
import re
import time
import itertools

# reverse complement taking function
def ReverseComplement(seq):
    for base in seq:
        if base not in 'ATCGatcg':
            print ("Error: NOT a DNA sequence")
            return None
    seq1 = 'ATCGTAGCatcgtagc'
    seq_dict = { seq1[i]:seq1[i+4] for i in range(16) if i < 4 or 8<=i<12 }
    return "".join([seq_dict[base] for base in reversed(seq)])
	
#GC content calculating function
def GCcontent(seq):
    GC = 0
    n = 0
    for base in seq:
        n = n + 1
        if base not in 'ATCGatcg':
            print ("Error: NOT a DNA sequence")
            return None
        if base in 'CGcg':
            GC = GC + 1
    return GC / n * 100

def TmcalculatorNEB(seq1, seq2):
    common_Tm=0
    pr1_Tm =0
    pr2_Tm=0
    elem_NEBpr1 = browser1.find_element_by_css_selector('#p1')
    elem_NEBpr1.send_keys(seq1)
    result_NEBpr1= browser1.find_element_by_css_selector('#tm1 > div:nth-child(2) > strong:nth-child(5)').text

    elem_NEBpr2 = browser1.find_element_by_css_selector('#p2')
    elem_NEBpr2.send_keys(seq2)
    result_NEBpr2= browser1.find_element_by_css_selector('#tm2 > div:nth-child(2) > strong:nth-child(5)').text

    result_NEB= browser1.find_element_by_css_selector('#ta').text

    match = re.search(r'Anneal at\n(\d+) °C', result_NEB)
    if match:
        common_Tm = float(match.group(1))
    match = 0
    match = re.search(r'Tm: (\d+) °C', result_NEBpr1)
    if match:
        pr1_Tm = float(match.group(1))
    match = 0
    match = re.search(r'Tm: (\d+) °C', result_NEBpr2)
    if match:
        pr2_Tm = float(match.group(1))
    elem_NEBpr1.clear()
    elem_NEBpr2.clear()
    return (common_Tm, pr1_Tm, pr2_Tm)

def deltaGcalculatorIDT(seq1, seq2):
#get deltaG for primer1 and primer2 from IDT's website
    delay = 3 # seconds
    elem = browser2.find_element_by_css_selector('textarea.form-control')
    elem.send_keys(seq1)
    button = browser2.find_element_by_css_selector('button.btn:nth-child(3)')
    button.click()
    time.sleep(delay)
    result= browser2.find_element_by_css_selector('div.well:nth-child(5) > span:nth-child(2)').text

    match = re.search(r'(\+?-?\d+\.?\d*) kcal\/mole', result)
    if match:
        Delta_G_primer1 = float(match.group(1))

    elem.clear()
    elem = browser2.find_element_by_css_selector('textarea.form-control')
    elem.send_keys(seq2)
    button = browser2.find_element_by_css_selector('button.btn:nth-child(3)')
    button.click()
    time.sleep(delay)
    result= browser2.find_element_by_css_selector('div.well:nth-child(5) > span:nth-child(2)').text
    match = re.search(r'(\+?-?\d+\.?\d*) kcal\/mole', result)
    if match:
        Delta_G_primer2 = float(match.group(1))
    else:
        time.sleep(3*delay)
        result= browser2.find_element_by_css_selector('div.well:nth-child(5) > span:nth-child(2)').text
        match = re.search(r'(\+?-?\d+\.?\d*) kcal\/mole', result)
        Delta_G_primer2 = float(match.group(1))
    return (Delta_G_primer1, Delta_G_primer2)
	
def homodimer_dG(seq):
#get homodimer deltaG of seq from IDT's website
    delay = 3 # seconds
    match = 0
    elem = browser2.find_element_by_css_selector('textarea.form-control')
    elem.send_keys(seq)
    button = browser2.find_element_by_css_selector('button.btn:nth-child(3)')
    button.click()
    time.sleep(delay)
    for i in range(1,300,1):
        time.sleep(1)
        result= browser2.find_element_by_css_selector('div.well:nth-child(5) > span:nth-child(2)').text
        match = re.search(r'(\+?-?\d+\.?\d*) kcal\/mole', result)
        if match:
            Delta_G = float(match.group(1))
            break
    if (i == 300):
        print('Timeout on IDTs website')
    
    elem.clear()
    return (Delta_G)

def heterodimer_dG(seq1, seq2):
#get heterodimer deltaG of seq from IDT's website, takes reverse complement of seq2
    delay = 3 # seconds
    elem = browser2.find_element_by_css_selector('textarea.form-control')
    #seq1 as input primer
    elem.send_keys(seq1)
    button = browser2.find_element_by_css_selector('button.btn:nth-child(4)')
    button.click()
    time.sleep(delay)
    elem_pair =  browser2.find_element_by_css_selector('#OAResults > div:nth-child(2) > div:nth-child(4) > div:nth-child(3) > div:nth-child(4) > div:nth-child(1) > div:nth-child(1) > textarea:nth-child(2)')
    #reverse complement of seq2 as second primer
    elem_pair.send_keys(ReverseComplement(seq2))
    button2 = browser2.find_element_by_css_selector('#OAResults > div:nth-child(2) > div:nth-child(4) > div:nth-child(3) > div:nth-child(5) > div:nth-child(1) > div:nth-child(1) > button:nth-child(2)')
    button2.click()
    time.sleep(delay)
    
    result= browser2.find_element_by_css_selector('#OAResults > div > div:nth-child(4) > div:nth-child(6) > span:nth-child(2)').text

    match = re.search(r'(\+?-?\d+\.?\d*) kcal\/mole', result)
    if match:
        Delta_G = float(match.group(1))
    elem.clear()
    return (Delta_G)
	
def unique(iterable):
    seen = set()
    for x in iterable:
        if x in seen:
            continue
        seen.add(x)
        yield x
	
	
#initialize browser and get NEB's Tm calculator
browser1 = webdriver.Chrome()
browser1.get('http://tmcalculator.neb.com/#!/')

#initialize browser and get IDT's OligoAnalyzer to calculate self dimers and heterodimerization
browser2 = webdriver.Chrome()
browser2.get('http://www.idtdna.com/calc/analyzer')


#optimize forward primer length to get Tm = 60
pr1_Tm = 0
#i = len(forward_sequence)
i = 10
while (pr1_Tm < 60 and i < len(forward_sequence)):
    elem_NEBpr1 = browser1.find_element_by_css_selector('#p1')
    elem_NEBpr1.send_keys(forward_sequence[0:i])
    result_NEBpr1= browser1.find_element_by_css_selector('#tm1 > div:nth-child(2) > strong:nth-child(5)').text
    match = 0
    match = re.search(r'Tm: (\d+) °C', result_NEBpr1)
    if match:
        pr1_Tm = float(match.group(1))
    i=i+1
    elem_NEBpr1.clear()
i=i-1
#print(i) #length of forward oligo with Tm = 60
#print(forward_sequence[0:i])

#optimize reverse primer length to get Tm = 60
pr2_Tm = 0
j = len(reverse_sequence)-10
while (pr2_Tm < 60 and j >= 0):
    elem_NEBpr2 = browser1.find_element_by_css_selector('#p2')
    elem_NEBpr2.send_keys(reverse_sequence[j:len(reverse_sequence)])
    result_NEBpr2= browser1.find_element_by_css_selector('#tm2 > div:nth-child(2) > strong:nth-child(5)').text
    match = 0
    match = re.search(r'Tm: (\d+) °C', result_NEBpr2)
    if match:
        pr2_Tm = float(match.group(1))
    j=j-1
    elem_NEBpr2.clear()
j=j+1
#print(len(reverse_sequence)-j) # length of reverse oligo with Tm = 60
#print(reverse_sequence[j:len(reverse_sequence)])

#read out primers
primer1_hom = forward_sequence[0:i]
primer2_hom = ReverseComplement(reverse_sequence[j:len(reverse_sequence)])
#read out homology Tms and homodimer, heterodimer deltaGs
common_Tm_hom, pr1_Tm_hom, pr2_Tm_hom = TmcalculatorNEB(primer1_hom, primer2_hom)
Delta_G_primer1_hom, Delta_G_primer2_hom = deltaGcalculatorIDT(primer1_hom, primer2_hom)

#add restriction enzyme sites and 6 bp AT-rich sequence at 5' end
primer1 =  restriction1 + primer1_hom
primer2 = ReverseComplement(restriction2) + primer2_hom

#change A to C T to G for lower Tm primer until Tms match
#compare Tms of homology region + restriction site + flanking sequence
#if difference is larger than 1C and lower Tm primer GC content is less than 100%
#then go through flanking sequence and replace the first A to C or the first T to G
#recalculate Tms to check if they equal -> if not, repeat

i=1
common_Tm, pr1_Tm, pr2_Tm = TmcalculatorNEB(flanking1 + primer1, flanking2 + primer2)

while abs(pr1_Tm - pr2_Tm)>1:
    if (pr1_Tm - pr2_Tm) < -1 : # flanking1 needs changes
        if GCcontent(flanking1) < 100:
            if flanking1.count('A')>flanking1.count('T'):
                flanking1 = flanking1.replace('A','C',1)
            else:
                flanking1 = flanking1.replace('T','G',1)
        else:
            print('GC content reached 100%')
    if (pr1_Tm - pr2_Tm) > 1 : #flanking2 needs changes
        if GCcontent(flanking2) < 100:
            if flanking2.count('A')>flanking2.count('T'):
                flanking2 = flanking2.replace('A','C',1)
            else:
                flanking2 = flanking2.replace('T','G',1)
        else:
            print('GC content reached 100%')
    
    common_Tm, pr1_Tm, pr2_Tm = TmcalculatorNEB(flanking1 + primer1, flanking2 + primer2)
    i=i+1
    if i > 10: break

#calculate homodimer deltaG for all possible permutations of flanking1, and store the best 5 sequences (lowest homodimer deltaG)
flanking1_list = [] #empty list
for base in unique(itertools.permutations(flanking1)):
    dG_homo = homodimer_dG(''.join(base)+ primer1)
    flanking1_list.append((''.join(base),-dG_homo))
#sort the list by dG_homo
flanking1_list.sort(key=lambda t: t[1]) 
flanking1_list = flanking1_list[0:5]

#calculate homodimer deltaG for all possible permutations of flanking2, and store the best 5 sequences (lowest homodimer deltaG)
flanking2_list = []
for base in unique(itertools.permutations(flanking2)):
    dG_homo = homodimer_dG(''.join(base)+ primer2)
    flanking2_list.append((''.join(base),-dG_homo))
#sort the list by dG_homo
flanking2_list.sort(key=lambda t: t[1]) 
flanking2_list = flanking2_list[0:5]
	
#calculate heterodimer deltaG for all 25 combinations of flanking1 and flanking2
#until get a list of 10 pairs with lower than 1C difference in Tm and lower than -10 kcal/mole heterodimer deltaG
flanking_list=[]
for i in range(0,5,1):
	for j in range(0,5,1):
		Tm, pr1_Tm, pr2_Tm=TmcalculatorNEB(flanking1_list[i][0] + primer1 ,flanking2_list[j][0] + primer2)
		#if abs(pr1_Tm - pr2_Tm) < 2 and len(flanking_list)<=10 :
		dG_hetero = heterodimer_dG(flanking1_list[i][0] + primer1, flanking2_list[j][0] + primer2)
		if dG_hetero > -10:
			flanking_list.append((flanking1_list[i][0], flanking2_list[j][0],flanking1_list[i][1], flanking2_list[j][1], dG_hetero))
	if len(flanking_list) > 10:
                break
flanking_list.sort(key=lambda t: t[4]) 

if len(flanking_list) < 1:
	print('No primer pairs found')
#read out full primer Tms for the top hit
Tm, pr1_Tm, pr2_Tm=TmcalculatorNEB(flanking_list[0][0] + primer1 ,flanking_list[0][1] + primer2)

	
	

print('forward primer: ',flanking_list[0][0]+primer1,' homology Tm', pr1_Tm_hom,' C, full primer Tm: ', pr1_Tm, ' C, homodimer dG:',flanking_list[0][2],'kcal/mole')
print('forward primer: ',flanking_list[0][1]+primer2,' homology Tm', pr2_Tm_hom,' C, full primer Tm: ', pr2_Tm, ' C, homodimer dG:',flanking_list[0][3],'kcal/mole')
print('common homology Tm: ', common_Tm_hom, 'C, full primer Tm: ', Tm, 'C, heterodimer deltaG: ', flanking_list[0][4],'kcal/mole')

#write output into file with gene name as filename
with open('%s.txt' % gene_name, 'w') as file:
	file.write('forward primer: ' + flanking_list[0][0] + primer1 + '\n')
	file.write('homology Tm: ' + str(pr1_Tm_hom) + '\n')
	file.write('homology region dG: ' + str(Delta_G_primer1_hom) + '\n')
	file.write('full primer Tm: ' + str(pr1_Tm) + '\n')
	file.write('homodimer dG: ' + str(flanking_list[0][2]) + ' kcal/mole\n')
	file.write('reverse primer: ' + flanking_list[0][1] + primer2 + '\n')
	file.write('homology Tm: ' + str(pr2_Tm_hom) + '\n')
	file.write('homology region dG: ' + str(Delta_G_primer2_hom) + '\n')
	file.write('full primer Tm: ' + str(pr2_Tm) + '\n')
	file.write('homodimer dG: ' + str(flanking_list[0][3]) + ' kcal/mole\n')
	file.write('common homology Tm: ' + str(common_Tm_hom) + '\n')
	file.write('full primer Tm: ' + str(Tm) + '\n')
	file.write('heterodimer deltaG: ' + str(flanking_list[0][4]))
