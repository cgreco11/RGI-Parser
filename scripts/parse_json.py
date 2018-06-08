import sys,json, os, glob, argparse
from collections import defaultdict
from pprint import pprint
import pandas as pd
from shutil import copyfile

parser = argparse.ArgumentParser(description = "Parse Resistance Genome Identifier (RGI) Output ")
parser.add_argument('-t', '--target_genomes', help = 'File List (ls *.json -- except card.json file)', required = 'True')
parser.add_argument('-a', '--all_hits', help = 'Generate an all hits file in addition to best hits file: Default is False', action = 'store_true')
parser.add_argument('-m', '--mapping_file', help = 'Drug Class Mapping Term', required = 'True')
args = parser.parse_args()

target_genomes = args.target_genomes
all_hits = args.all_hits
mapping_file = args.mapping_file

directory = 'AMR_STATISTICS/'
if not os.path.exists(os.path.dirname(directory)):
    os.makedirs(os.path.dirname(directory))

pwd = os.getcwd()

genomes_included = [] # "ex.json"
with open(target_genomes, 'r') as genome_list:
    for genomes in genome_list:
        genomes_included.append(genomes.strip())


if 'card.json' in genomes_included:
	print
	print "ERROR: Remove card.json from the genome list file"
	print
	sys.exit()

def generate_all_hits(genome_list):
    hits_dict = defaultdict(lambda: defaultdict(dict))
    for i in genome_list:
        dir_name = i.split('.json')[0]
        if not os.path.exists(pwd + '/AMR_STATISTICS/' + dir_name):
            os.makedirs(pwd+ '/AMR_STATISTICS/' + dir_name)
            copyfile(pwd + '/' + i, pwd + '/AMR_STATISTICS/' + dir_name + '/' + i)
        filename = pwd + '/AMR_STATISTICS/' + dir_name + '/' + i
        with open(filename, 'r') as json_file:
            data = json.load(json_file)
        for k,v in data.items():
            count = 0
            if k != '_metadata':
                best_hit = {} #Best Hit for that Contig, not AMR Gene in that Genome
                for a,b in v.items():
                    if count == 0:
                        try:
                            best_hit[k] = [b['ARO_name'].encode('utf-8'), b['perc_identity'], b['orf_dna_sequence'].encode('utf-8'), b['bit-score']]
                        except:
                            print k
                        count = count + 1
                    else:
                        if b['bit-score'] > best_hit[k][-1]:
                            best_hit[k] = [b['ARO_name'].encode('utf-8'), b['perc_identity'], b['orf_dna_sequence'].encode('utf-8'), b['bit-score']]
                        else:
                            pass
                hits_dict[dir_name].update(best_hit)
    return hits_dict

def generate_best_hits(all_hits_dictionary):
    hits_dict = defaultdict(dict)
    for k,v in all_hits_dict.items():
        aro_dict = {}
        for a,b in v.items():
            aro_name =  b[0]
            if aro_name not in aro_dict:
                aro_dict[aro_name] = {'Percent Identity' : b[1],"Sequence" : b[2], "Bit Score": b[3]}
            else:
                if b[1] > float(str(aro_dict[aro_name]['Percent Identity']).replace("*", "")):
                    aro_dict[aro_name] = {'Percent Identity' : str(b[1]) + '*', 'Sequence' : b[2], "Bit Score" : b[3]}
        hits_dict[k] = aro_dict
    #print hits_dict
    return hits_dict

all_hits_dict = generate_all_hits(genomes_included)
best_hits_dict = generate_best_hits(all_hits_dict)

if all_hits:
    for k,v in all_hits_dict.items():
        all_hits_filename = pwd + '/AMR_STATISTICS/' + k + '/' + k +"_all_hits.fa"
        with open(all_hits_filename, 'w') as f:
            for a,b in v.items():
                f.write('>' + b[0] + '_' + str(b[1]) + '\n')
                f.write(b[2].strip() + '\n')

for k,v in best_hits_dict.items():
    best_hits_filename = pwd + '/AMR_STATISTICS/' + k + '/' + k + "_best_hits.fa"
    with open(best_hits_filename, 'w') as f:
        for a,b in v.items():
            f.write(">" + a + '_' + str(b['Percent Identity']) + '\n')
            f.write(b['Sequence'].strip() + '\n')


def top_hits_generator(best_hits):
    top_hits = defaultdict(lambda: defaultdict(dict))
    for k,v in best_hits.items():
        for a,b in v.items():
            if b['Percent Identity'] != 100.0 and b["Percent Identity"] != "100.0*":
                b['Percent Identity'] = '[' + str(b['Percent Identity']) + ']'
            top_hits[k][a] = b['Percent Identity']
    return top_hits

top_hits_dict = top_hits_generator(best_hits_dict)

def create_top_hits_file(top_hits, map_file):
    main_class = {}
    sub_class = {}
    Type = {}
    classification = {}
    curated = {}
    mapping_file = pd.read_csv(map_file, sep = '\t', index_col = 0)
    map_dict = pd.DataFrame.to_dict(mapping_file, orient = 'index')
    for k,v in map_dict.items():
        main_class[k] = v['Main Class']
        sub_class[k] = v['Sub Class']
        Type[k] = v['Type']
        classification[k] = v['Classification']
        curated[k] = v['Curated']
    df = pd.DataFrame.from_dict(top_hits)
    df = pd.DataFrame.fillna(df, value = 0)
    df.index.name = 'Gene'
    df = df.reset_index()
    df['Type'] = df['Gene'].map(Type)
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    df['Classification'] = df['Gene'].map(classification)
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    df['Sub Class'] = df["Gene"].map(sub_class)
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    df['Main Class'] = df['Gene'].map(main_class)
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    df['Curated Name'] = df['Gene'].map(curated)
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return df

number_of_genomes = len(best_hits_dict.keys())

def status_of_arg(row):
    if row['Count'] == number_of_genomes:
        return 'Core'
    if row['Count'] == 1:
        return 'Singleton'
    else:
        return 'Shared'

top_hits_file_df = create_top_hits_file(top_hits_dict, mapping_file)
top_hits_file_df = top_hits_file_df.fillna("NA")
top_hits_file_df = top_hits_file_df.set_index('Gene') # RESET THIS AND REMOVE GENE/ REPLACE WITH CURATED
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#Genome Hits
genome_hits_series = top_hits_file_df.astype(bool).sum(axis = 0)[5:]
genome_hits_df = genome_hits_series.to_frame()
genome_hits_df.index.name = 'Genome'
genome_hits_df.columns = ['Total Count']

main_class_values = top_hits_file_df['Main Class'].unique()
main_class_df = pd.DataFrame(genome_hits_series.to_frame())
for i in main_class_values:
    main_class_df[i] = top_hits_file_df.loc[top_hits_file_df['Main Class'] == i].astype(bool).sum(axis = 0)[5:].to_frame()
main_class_df.columns.values[0] = 'Main Class Total'

sub_class_values = top_hits_file_df['Sub Class'].unique()
sub_class_df = pd.DataFrame(genome_hits_series.to_frame())
for i in sub_class_values:
    sub_class_df[i] = top_hits_file_df.loc[top_hits_file_df['Sub Class'] == i].astype(bool).sum(axis = 0)[5:].to_frame()
sub_class_df.columns.values[0] = 'Sub Class Total'

classification_values = top_hits_file_df['Classification'].unique()
classification_df = pd.DataFrame(genome_hits_series.to_frame())
for i in classification_values:
    classification_df[i] = top_hits_file_df.loc[top_hits_file_df['Classification'] == i].astype(bool).sum(axis = 0)[5:].to_frame()
classification_df.columns.values[0] = 'Classification Total'

type_values = top_hits_file_df['Type'].unique()
type_df = pd.DataFrame(genome_hits_series.to_frame())
for i in type_values:
    type_df[i] = top_hits_file_df.loc[top_hits_file_df['Type'] == i].astype(bool).sum(axis = 0)[5:].to_frame()
type_df.columns.values[0] = 'Type Total'

'''with open('Genome_Stats.txt', 'w') as handle:
    genome_hits_df.to_csv(handle, sep = '\t')
    main_class_df.to_csv(handle, sep = '\t')
    sub_class_df.to_csv(handle, sep = '\t')
    classification_df.to_csv(handle, sep = '\t')
    type_df.to_csv(handle, sep = '\t')
'''
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#ARG Hits
amr_hits = top_hits_file_df.iloc[:,5:] #Remove MetaData so we can sum across rows
amr_hits_series = amr_hits.astype(bool).sum(axis = 1)
amr_hits_df = amr_hits_series.to_frame()
amr_hits_df.columns = ['Count']
amr_hits_df['Curated Name'] = top_hits_file_df['Curated Name']
amr_hits_df['Main Class'] = top_hits_file_df['Main Class']
amr_hits_df['Sub Class'] = top_hits_file_df['Sub Class']
amr_hits_df['Classification'] = top_hits_file_df['Classification']
amr_hits_df['Type'] = top_hits_file_df['Type']
amr_hits_df['ARG Status Across Genomes'] = amr_hits_df.apply(lambda row: status_of_arg (row),axis=1)

cols = list(amr_hits_df)
cols[0], cols[1] = cols[1], cols [0]
amr_hits_df = amr_hits_df.ix[:,cols]
singleton_genome = {}
for idx,row in amr_hits_df.iterrows():
    if row['ARG Status Across Genomes'] == 'Singleton':
        arg_series = top_hits_file_df.loc[idx][5:]
        singleton_genome[idx] = arg_series[arg_series != 0].index.tolist()[0]

top_hits_file_df = top_hits_file_df.reset_index()
top_hits_file_df['Singleton Genome'] = top_hits_file_df['Gene'].map(singleton_genome)
top_hits_file_df = top_hits_file_df.set_index("Gene")
top_hits_file_df = top_hits_file_df.sort_values(['Main Class', 'Sub Class', 'Classification', 'Type', 'Curated Name'])

amr_hits_df['Singleton Genome'] = top_hits_file_df['Singleton Genome']
amr_hits_df = amr_hits_df.sort_values(['Main Class', 'Sub Class', 'Classification', 'Type', 'Curated Name'])
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
#Top AMR hits
del top_hits_file_df['Singleton Genome']
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# combined_top_hits.txt
combined_df = top_hits_file_df.copy()

del combined_df['Classification']
del combined_df['Type']

for idx, row in combined_df.iterrows():
    if row['Main Class'] == 'Beta-lactams':
        combined_df.loc[idx, 'Parent Allele'] = row['Curated Name'].rsplit('-', 1)[0]
    else:
        combined_df.loc[idx, 'Parent Allele'] = row['Curated Name']

cols = combined_df.columns.tolist()
cols.insert(3, cols[-1])
del cols[-1]
cols.insert(4, cols[0])
del cols[0]
combined_df = combined_df[cols]

info_df = combined_df.copy()

del info_df['Main Class']
del info_df['Sub Class']
del info_df['Parent Allele']
del info_df['Curated Name']



info_df = info_df.transpose()
info_dict = pd.DataFrame.to_dict(info_df, orient = 'index')
for k,v in info_dict.items():
    for a,b in v.items():
        if str(b) == '0':
            info_dict[k][a] = ''
        elif str(b) == '100.0':
            info_dict[k][a] = a
        elif b == '100.0*':
            info_dict[k][a] = a + "*"
        else:
            info_dict[k][a] = '[' + a + ']'

info_df_2 = pd.DataFrame.from_dict(info_dict, orient = 'index')
info_df_2 = info_df_2.transpose()

info_df_2['Parent Allele'] = combined_df['Parent Allele']
cols = info_df_2.columns.tolist()
cols = cols[-1:] + cols[:-1]
info_df_2 = info_df_2[cols]
#FIX WHEN SUB CLASS INFORMATION IS FINALIZED
#info_df_2['Sub Class'] = combined_df['Sub Class']
#cols = info_df_2.columns.tolist()
#cols = cols[-1:] + cols[:-1]
#info_df_2 = info_df_2[cols]
info_df_2['Main Class'] = combined_df['Main Class']
cols = info_df_2.columns.tolist()
cols = cols[-1:] + cols[:-1]
info_df_2 = info_df_2[cols]

indices = info_df_2.index.tolist()
for i in indices:
    if isinstance(info_df_2.loc[i, 'Parent Allele'], basestring):
        pass
    else:
        info_df_2.loc[i, 'Parent Allele'] = i
        info_df_2.loc[i, "Main Class"] = i


info_df_2 = info_df_2.groupby(["Parent Allele", 'Main Class']).agg(' '.join)
info_df_2 = info_df_2.sort_index(axis = 0, level = 1, sort_remaining = True)
info_df_2 = info_df_2.transpose()

cols = info_df_2.columns.tolist()
for i in cols:
    info_df_2[i] = info_df_2[i].str.strip()

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
missing_hits = []
for idx, row in top_hits_file_df.iterrows():
    if row['Curated Name'] == 'NA':
        missing_hits.append(idx)


with open('AMR_STATISTICS/Genome_Stats.txt', 'w') as handle:
    genome_hits_df.to_csv(handle, sep = '\t')
    main_class_df.to_csv(handle, sep = '\t')
    sub_class_df.to_csv(handle, sep = '\t')
    classification_df.to_csv(handle, sep = '\t')
    type_df.to_csv(handle, sep = '\t')

amr_hits_df.to_csv("AMR_STATISTICS/ARG_Stats.txt", sep = '\t')


top_hits_file_df.to_csv('AMR_STATISTICS/Top_AMR_Hits.txt', sep = '\t')
info_df_2.to_csv('AMR_STATISTICS/Combined_Top_Hits.txt', sep = '\t')

with open('AMR_STATISTICS/Missing_Genes.txt', 'w') as f:
    for line in missing_hits:
        f.write(line.strip() + '\n')
