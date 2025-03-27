import sys
import re
import pandas as pd
import numpy as np
from scipy.stats import beta
from scipy import stats
import os
import subprocess

def translate_bases(ref, depth, match):
    ref_lc = ref.lower()
    bases = ['A','a','T','t','C','c','G','g', "I", "i", "D", "d"]
    pos_n = {}
    count = {}
    for base in bases:
        count[base] = 0
    for base in bases:
        pos_n[base] = []
    I = 0; D = 0
    in_base= 0; del_base = 0; #mis_base = 0
    i = 0
    pos_q = 0
    end_reads = {}
    while(i < len(match)):
        if match[i] == '^':
            i += 2
        elif match[i] == "$":
            i += 1
        elif match[i] == ".":
            count[ref] += 1
            pos_n.setdefault(ref,[]).append(pos_q)
            pos_q += 1
            i += 1
        elif match[i] == ",":
            count[ref_lc] += 1
            pos_n.setdefault(ref_lc,[]).append(pos_q)
            pos_q += 1
            i += 1
        elif match[i] in bases:
            count[match[i]] += 1
            pos_n.setdefault(match[i],[]).append(pos_q)
            pos_q += 1
            i += 1
        elif match[i] == '+':
            I += 1
            if re.match(r"\d",match[i+2]):
                if match[i+3].isupper():
                    count["I"] += 1
                else: 
                    count["i"] +=1 
                c = 10*int(match[i+1]) + int(match[i+2])
                i += c+3
                pos_q += c+3
            else:
                if match[i+2].isupper():
                    count["I"] += 1
                else: 
                    count["i"] +=1 
                c = int(match[i+1])
                i += c+2
                pos_q += c+2
            in_base += c
        elif match[i] == '-':
            D += 1
            if re.match(r'\d',match[i+2]):
                if match[i+3].isupper():
                    count["D"] += 1
                else: 
                    count["d"] +=1 
                c = 10*int(match[i+1])+ int(match[i+2])
                i += c+3
            else:
                if match[i+2].isupper():
                    count["D"] += 1
                else: 
                    count["d"] +=1 
                c = int(match[i+1])
                i += c+2
            del_base += c
        elif match[i] == "*":
            D += 1
            i += 1
        else:
            sys.stderr.write("ERROR: " + str(i) + " " + str(match[i]) + "\n")
            sys.exit(0)

    return count, pos_n, I, D



def translate_qualities(ref, alt, pos_n, qualities): 
    #print(pos_n)
    quality_ref = ""
    quality_alt = ""
    for position in pos_n[ref]+pos_n[ref.lower()]:
        quality_ref += qualities[position]
    for position in pos_n[alt]+pos_n[alt.lower()]:
        quality_alt += qualities[position]
    quality_ref = re.sub(r'!','',quality_ref)
    quality_ref = re.sub(r'\\','',quality_ref)
    quality_alt = re.sub(r'!','',quality_alt)
    quality_alt = re.sub(r'\\','',quality_alt)
    return quality_ref, quality_alt


def compute_table(ref, alt, ref_qualities, alt_qualities):
    qualities = ref_qualities + alt_qualities
    ref_length = len(ref_qualities)
    alt_length = len(alt_qualities)
    length = ref_length + alt_length
    #print(length, ref_length)
    matrix = np.zeros((length + 1, length + 1))
    #compute matrix
    matrix[0,0] = 1
    for i in range(length):
        prob = 1- 10**-(ord(qualities[i])-33)
        #print(prob)
        if i >= ref_length:
            prob = 1-prob
        for j in range(i+1):
            matrix[j, i+1] += matrix[j, i] * prob
            matrix[j+1, i+1] += matrix[j, i] * (1-prob)
    scores = matrix[:,-1]
    return scores

def compute_MAF_and_CI(scores):
    length = len(scores)
    MAF = np.argmax(scores)/length
    scores_index = np.nonzero(scores)
    nonzero_scores = list(scores[scores_index])
    start_index = scores_index[0][0]
    if scores_index[0][0] != 0:
        start_index = start_index -1
        nonzero_scores = [0] + nonzero_scores
    if scores_index[0][-1] != len(scores) -1:
        nonzero_scores = nonzero_scores + [0]
    #expanding scores
    #print(sum(nonzero_scores))
    new_scores = []
    for i in range(len(nonzero_scores)-1):
        first = nonzero_scores[i]
        second = nonzero_scores[i+1]
        new_scores.append(first)
        new_scores.append(first + (second-first)/5)
        new_scores.append(first + (second-first)/5*2)
        new_scores.append(first + (second-first)/5*3)
        new_scores.append(first + (second-first)/5*4)
    new_scores.append(second)
    sum_scores = sum(new_scores)
    #print(sum_scores)
    #compute ci by adding up the scores from left and from right
    end_index = start_index + len(new_scores) - 1
    agg = 0
    for i in range(len(new_scores)):
        agg += new_scores[i]
        if agg >= 0.025*sum_scores:
            lower_CI = (start_index + i/5)/length
            break
    agg = 0
    for i in range(len(new_scores)):
        agg += new_scores[-(i+1)]
        if agg >= 0.025*sum_scores:
            upper_CI = (start_index+(end_index-i)/5)/length
            break
    return MAF, lower_CI, upper_CI


def clopper_binom_interval(success, total, confint=0.95):
    quantile = (1 - confint) / 2.
    lower = beta.ppf(quantile, success, total-success+1)
    upper = beta.ppf(1 - quantile, success+1, total-success)
    if np.isnan(lower):
        lower = 0
    return lower, upper

def wilson_binom_interval(success, total, alpha = 0.05):
    q_ = success / total
    crit = stats.norm.isf( alpha / 2.)
    crit2 = crit**2
    denom = 1 + crit2 / total
    center = (q_ + crit2 / (2 * total)) / denom
    dist = crit * np.sqrt(q_ * (1. - q_) / total + crit2 / (4. * total**2))
    dist /= denom
    ci_low = center - dist
    ci_upp = center + dist
    return ci_low, ci_upp
    


def parse_samtools(output):
    f = open(file, "r")
    gt_pos = int(gt_position.split(":")[-1])       
    for line in f:
        #print(line.split("\t"))
        if len(line.split("\t")) > 0:
            items = line.strip().split("\t")
            chrom = items[0]
            pos = int(items[1])
            #print(pos) 
            if pos >  gt_pos and pos < gt_pos + len(gt_ref):
                continue
            position = str(chrom) + ":" + str(pos)
            #print(position)
            ref = items[2].upper()
            depth = int(items[3])
            if depth > 0:
                match = items[4]
                quality = items[5]
                count, pos_n, in_base, del_base = translate_bases(ref, depth, match)
                #+- 5bp or this bp
                if position == gt_position:
                    if len(gt_ref) > len(gt_alt):
                        num_alt = del_base
                        num_ref = count[ref] + count[ref.lower()]
                    elif len(gt_ref) < len(gt_alt):
                        num_alt = in_base
                        num_ref = count[ref] + count[ref.lower()]
                    else:
                        alt = gt_alt
                        num_ref = count[ref] + count[ref.lower()]
                        num_alt = count[alt] + count[alt.lower()]
                    gt_ref_count = num_ref
                    gt_alt_count = num_alt
                    if gt_ref_count + gt_alt_count == 0:
                        return
                    else:
                        gt_maf = gt_alt_count / (gt_ref_count + gt_alt_count)
                        gt_lower, gt_upper = wilson_binom_interval(num_alt, num_alt + num_ref, alpha = 0.05)
                else:
                    allbases = {}
                    for base in ['A','T','C','G']:
                        count.setdefault(base, 0)
                        count.setdefault(base.lower(),0)
                        allbases[base] = count[base]+count[base.lower()] 
                    allbases["in"] = in_base
                    allbases["del"] = del_base
                    keys = sorted(allbases,key=lambda k:allbases[k],reverse=True)
                    if keys[0] == ref:
                        alt1 = keys[1]
                    else:
                        alt1 = keys[0]
                    num_ref = allbases[ref]
                    num_alt = allbases[alt1]                
                    ci_lower, ci_upper = wilson_binom_interval(num_alt, num_alt + num_ref, alpha = 0.05)
                    if pos < gt_pos:
                        pre_upper_cis.append(ci_upper)
                    else:
                        pos_upper_cis.append(ci_upper)
    if gt_maf is None:
        return
    f.close()
    return gt_ref_count, gt_alt_count, gt_maf, gt_lower, gt_upper, pre_upper_cis, pos_upper_cis



def main(argv):
    if len(argv) != 5:
        sys.stderr.write("usage: " + argv[0] + "<samtools_output> <chrom:pos> <REF> <ALT>\n")
        sys.exit(2)
    file = argv[1]
    gt_position = argv[2]
    gt_ref = argv[3]
    gt_alt = argv[4]
    output =  parse_samtools(gt_position, gt_ref, gt_alt, file)
    if output is None:
        gt_ref_count = np.nan
        gt_alt_count = np.nan
        gt_maf = np.nan
        gt_lower = np.nan
        gt_upper = np.nan
        greater = np.nan
    else: 
        gt_ref_count, gt_alt_count, gt_maf, gt_lower, gt_upper, pre_upper_cis, pos_upper_cis = output
        greater = ""
        for value in pre_upper_cis:
            if gt_lower <= value:
                greater += "F"
            else:
                greater += "P"
        greater += "_"
        for value in pos_upper_cis:
            if gt_lower <= value:
                greater += "F"
            else:
                greater += "P"
    #writing  
    print(gt_ref_count, gt_alt_count, gt_maf, gt_lower, gt_upper, greater)


#data = pd.read_csv('libraries.csv')
#data = pd.read_csv('random96.csv', dtype=str)
data = pd.read_csv('spotcheck.csv', dtype=str)
#data = pd.read_csv('jareds.csv')

df_tas = data
#df_tas = pd.DataFrame(data)

#TAS_list = df_tas.values.tolist()


#df_positions = pd.read_csv('snp_locations_truesnps.csv')
#df_positions = pd.read_csv('snp_locations_justdel.csv')
#df_positions = pd.read_csv('snp_locations.csv')
#df_positions = pd.read_csv('Clinvar1-96.csv')
#df_positions = pd.read_csv('clinvar2_96_random.csv', dtype=str)
df_positions = pd.read_csv('clinvar2_spotcheck.csv', dtype=str)

chromosome = df_positions['CHROM'].astype(str)

gt_position = df_positions['POS'].astype(str)

ref_gt = df_positions['REF']

alt_gt = df_positions['ALT']

fam = df_positions['FAMILY']

#dm_maf = df_positions['maf_from_DM']

#chromosome=['chr11', 'chr11', 'chr5', 'chr22', 'chr11', 'chr15', 'chr2', 'chr20', 'chr5', 'chr10', 'chr9', 'chr12', 'chr2']
#gt_position = ['1172511', '119475487', '141441861', '20602713', '180318', '42528253', '100107183', '5066284', '139051170', '49739461', '35682258', '11120964', '31217463']
#ref_gt = [ "A", "T", "C", "T", "C", "C", "C", "G", "G", "C", "G", "G", "G", "G"]
#alt_gt = [ "C", "C", "T", "C", "T", "T", "G", "A", "A", "A", "A", "C", "C", "A"]


#chromosome=['chr16', 'chr19', 'chr3', 'chr4', 'chr10', 'chr4', 'chr2', 'chr17']
#gt_position = ['79595584', '54055546', '179021806', '177961089', '124446608', '47406821', '170815568', '1727101']
#ref_gt = ["A", "C", "C", "T", "G", "C", "C", "T"]
#alt_gt = ["G", "G", "G", "C", "A", "T", "G", "C"]


#chromosome=['chr6', 'chr16', 'chr19', 'chr3', 'chr4', 'chr10', 'chr4', 'chr2', 'chr17']
#gt_position = ['29268945', '79595584', '54055546', '179021806', '177961089', '124446608', '47406821', '170815568', '1727101']
#ref_gt = ["G", "A", "C", "C", "T", "G", "C", "C", "T"]
#alt_gt = ["A", "G", "G", "G", "C", "A", "T", "G", "C"]


sample = []

chromo = []

poso = []

reference_allele = []

alt_allele = []

low = []

ma_fract = []

up = []

family = []

variant_id = []

deletions = []

for i in range(len(df_tas['tas'])):
    for j in range(len(chromosome)):
        try:
        #subprocess.run(["samtools", "mpileup", "-r", chromosome[j]+":"+gt_position[j]+"-"+gt_position[j], "-f" ,"/pl/active/Breuss-Lab/reference_genomes/hg38.fa", "-Q", "15", "-q0", "-AB",
        #    "-d100000", "tas_recal_bams/"+df_tas['tas'][i]+".recal.bam", "-o" ,"pileup_tas_jan/"+chromosome[j]+"_"+gt_position[j]+".txt"])
            print('reading_pileup')

            pileup=open("pileup/"+df_tas['tas'][i]+"_"+chromosome[j]+"_"+gt_position[j]+".txt")
            content = pileup.read()

            if int(os.path.getsize("pileup/"+df_tas['tas'][i]+"_"+chromosome[j]+"_"+gt_position[j]+".txt")) == 0:
                print("PILEUP IS EMPTY")
                sample.append(df_tas['tas'][i])
                chromo.append(chromosome[j])
                poso.append(gt_position[j])

                num_ref = 'NA'
                reference_allele.append(num_ref)

                num_alt='NA'
                alt_allele.append(num_alt)

                ci_lower='NA'
                low.append(ci_lower)

                maf='NA'
                ma_fract.append(maf)

                ci_upper='NA'
                up.append(ci_upper)

                family.append(fam[j])

                deletion_AF = 'NA'
                deletions.append(deletion_AF)
                continue
                #print(chromosome[j]+"_"+gt_position[j])
                #raise Exception("Pileup is empty")

    #for j in range(len(chromosome)):
        #subprocess.run(["samtools", "mpileup", "-r", chromosome[j]+":"+gt_position[j]+"-"+gt_position[j], "-f" ,"/projects/jonpitsch@xsede.org/project2/pasm_pipe/data/hg38.fa", "-Q", "15", "-q0", "-AB", 
	 #   "-d100000", "j_bqsr/"+df_tas['jareds'][i]+".recal.bam", "-o" ,"pileup_jareds/"+chromosome[j]+"_"+gt_position[j]+".txt"])

        #pileup=open("j_pileup/"+df_tas['jareds'][i]+"_"+chromosome[j]+"_"+gt_position[j]+".txt")
        #content = pileup.read()
            items = content.rstrip().split("\t")

            print(df_tas['tas'][i])
            print(chromosome[j]+"_"+gt_position[j])

            chrom = chromosome[j]
            pos = int(gt_position[j])
            gt_ref = ref_gt[j]
            gt_alt = alt_gt[j]

            if len(gt_ref) > 1:
                gt_ref = ref_gt[j][0]
            else:
                gt_ref = ref_gt[j]

            if len(gt_alt) > 1:
                gt_alt = alt_gt[j][0]
            else:
                gt_alt = alt_gt[j]

            #vari = str(chrom+'_'+pos+'_'+gt_ref+'_'+gt_alt)
            #variant_id.append(vari)
            #print(variant_id)

            len_depth = len(items[3::3])
            ref = items[2].upper()

            for k in range(0,len_depth):
                try:
                    depth = int(items[3::3][k])
                    if depth > 0:
                        match = items[4::3][k]
                        quality = items[5::3][k]
                        count, pos_n, in_base, del_base = translate_bases(ref, depth, match)

            #num_ref = count[gt_ref] + count[gt_ref.lower()]
                    num_ref = count[ref] + count[ref.lower()]
                    num_alt = count[gt_alt] + count[gt_alt.lower()]

                    print('COUNTED_REF_AND_ALT')


            #print(num_ref, num_alt)
                    maf = num_alt / (num_alt+num_ref)
                    
                    deletion_AF = del_base/(num_ref+del_base)

                    deletions.append(deletion_AF)

                    print("DELETION_AF")
            #print(maf)
                    ci_lower, ci_upper = wilson_binom_interval(num_alt, num_alt + num_ref, alpha = 0.05)
            #print(ci_lower, ci_upper)

                    print("RAN_WILSON_BINOM_INTERVAL")

            #import scipy.stats as stats
                    num_ref_for = count[gt_ref]
            #print('This is count ref upper: {}'.format(count[gt_ref]))
                    num_ref_rev = count[gt_ref.lower()]
            #print('This is count ref lower: {}'.format(count[gt_ref.lower()]))
                    num_alt_for = count[gt_alt]
            #print('This is count alt upper: {}'.format(count[gt_alt]))
                    num_alt_rev = count[gt_alt.lower()]
            #print('This is count alt lower: {}'.format(count[gt_alt.lower()]))

                    print("CALCULATED_CI_VALUES")
            #oddsratio, pvalue = stats.fisher_exact([[num_ref_for, num_ref_rev], [num_alt_for, num_alt_rev]])

            #quality_ref, quality_alt = translate_qualities(gt_ref, gt_alt, pos_n, quality)

            #scores = compute_table(gt_ref, gt_alt, quality_ref, quality_alt)

            #x = compute_MAF_and_CI(scores)
            #print("This is the sample:{}".format(df_tas['tas'][i]))
            #print(chromosome[j]+":"+(gt_position[j]))
            #print("This is MAF and CI:{}".format(x))
            #print(num_ref, num_alt)
            #maf = num_alt / (num_alt+num_ref)
            #print(maf)
            #ci_lower, ci_upper = wilson_binom_interval(num_alt, num_alt + num_ref, alpha = 0.05)
            #print(ci_lower, ci_upper)

            #df_tas_results = pd.DataFrame()


            #df_tas_results.columns = ['Sample', 'chrom', 'pos', 'ref_count', 'alt_count', 'lower CI', 'maf', 'upper CI']

            #df_tas_results['Sample'] = [df_tas['tas'][i]]

                    sample.append(df_tas['tas'][i])

            #df_tas_results['chrom'] = [chromosome[j]]

                    chromo.append(chromosome[j])

            #df_tas_results['pos'] = [gt_position[j]]

                    poso.append(gt_position[j])

            #df_tas_results['ref_count'] = [num_ref]

                    reference_allele.append(num_ref)

            #df_tas_results['alt_count'] = [num_alt]

                    alt_allele.append(num_alt)

                    print("APPENDED_SAMPLE_CHROM_POS_REF_ALT_ARRAYS")

            #df_tas_results['lower CI'] = [ci_lower]
                    print(ci_lower)
                    low.append(ci_lower)

            #df_tas_results['maf'] = [maf]
                    print(maf)
                    ma_fract.append(maf)

            #df_tas_results['upper CI'] = [ci_upper]
                    print(ci_upper)
                    up.append(ci_upper)

                    family.append(fam[j])

                    print("FINISHED_APPENDING_ALL_ARRAYS")
            
                except:
                    sample.append(df_tas['tas'][i])
                    chromo.append(chromosome[j])
                    poso.append(gt_position[j])

                    num_ref = 'NA'
                    reference_allele.append(num_ref)

                    num_alt='NA'
                    alt_allele.append(num_alt)

                    ci_lower='NA'
                    low.append(ci_lower)

                    maf='NA'
                    ma_fract.append(maf)

                    ci_upper='NA'
                    up.append(ci_upper)

                    family.append(fam[j])

                    deletion_AF = 'NA'
                    deletions.append(deletion_AF)
                    pass

        except:
            pass

            #df_tas_results.to_csv('tas_csvs/'+df_tas['tas'][i]+"_"+chromosome[j]+":"+gt_position[j]+'_tas_results.csv')

df_tas_results = pd.DataFrame()

df_tas_results['Family'] = family

df_tas_results['Sample'] = sample

df_tas_results['chrom'] = chromo

df_tas_results['pos'] = poso

df_tas_results['ref_count'] = reference_allele

df_tas_results['alt_count'] = alt_allele

df_tas_results['lower CI'] = low

df_tas_results['maf'] = ma_fract

#df_tas_results['maf_from_DM'] = dm_maf

df_tas_results['upper CI'] = up

#df_tas_results['deletions_AF'] = deletions

#df_tas_results['variant_id'] = variant_id

df_tas_results.to_csv('Clinvar2_spotcheck_TAS_validation_redo.csv', index=False)
#df_tas_results.to_csv('Clinvar2_random_TAS_validation_redo.csv', index=False)























