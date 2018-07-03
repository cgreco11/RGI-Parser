import sys,json, os, glob, argparse, math
from collections import defaultdict
from pprint import pprint
import pandas as pd
from shutil import copyfile

def make_dir(directory):
    if not os.path.exists(os.path.dirname(directory)):
        os.makedirs(os.path.dirname(directory))

def get_json_files(target_genomes):
    genomes_included = [] # "ex.json"
    with open(target_genomes, 'r') as genome_list:
        for genomes in genome_list:
            genomes_included.append(genomes.strip())

    if 'card.json' in genomes_included:
    	print
    	print "ERROR: Remove card.json from the genome list file"
    	print
    	sys.exit()
    return genomes_included

def generate_all_hits(genome_list):
    pwd = os.getcwd()
    hits_dict = defaultdict(lambda: defaultdict(dict))
    contig_dict = defaultdict(lambda: defaultdict(dict))
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
                contig_hit = {}
                for a,b in v.items():
                    if count == 0:
                        try:
                            best_hit[k] = [b['ARO_name'].encode('utf-8'), b['perc_identity'], b['orf_dna_sequence'].encode('utf-8'), b['bit-score']]
                            contig_hit[b['orf_From']] = {b["ARO_name"] : b['perc_identity']}
                        except:
							#print k
                            pass
                        count = count + 1
                    else:
                        if b['bit-score'] > best_hit[k][-1]:
                            best_hit[k] = [b['ARO_name'].encode('utf-8'), b['perc_identity'], b['orf_dna_sequence'].encode('utf-8'), b['bit-score']]
                            contig_hit[b['orf_From']] = {b["ARO_name"] : b['perc_identity']}
                        else:
                            pass
                hits_dict[dir_name].update(best_hit)
                contig_dict[dir_name].update(contig_hit)
    return hits_dict, contig_dict

def generate_best_hits(all_hits_dictionary):
    hits_dict = defaultdict(dict)
    for k,v in all_hits_dictionary.items():
        aro_dict = {}
        for a,b in v.items():
            aro_name =  b[0]
            if aro_name not in aro_dict:
                aro_dict[aro_name] = {'Percent Identity' : b[1],"Sequence" : b[2], "Bit Score": b[3]}
            else:
                if b[1] > float(str(aro_dict[aro_name]['Percent Identity']).replace("*", "")):
                    aro_dict[aro_name] = {'Percent Identity' : str(b[1]) + '*', 'Sequence' : b[2], "Bit Score" : b[3]}
        hits_dict[k] = aro_dict
    return hits_dict

def write_all_hits(hits_dict):
    for k,v in hits_dict.items():
        filename = os.getcwd() + "/AMR_STATISTICS/" + k + "/" + k + "_all_hits.fa"
        with open(filename, 'w') as f:
            for a,b in v.items():
                f.write(">" + b[0] + "_" + str(b[1]) + "\n")
                f.write(">" + b[2].strip() + "\n")

def write_best_hits(best_hits_dict):
    for k,v in best_hits_dict.items():
        best_hits_filename = os.getcwd() + '/AMR_STATISTICS/' + k + '/' + k + "_best_hits.fa"
        with open(best_hits_filename, 'w') as f:
            for a,b in v.items():
                f.write(">" + a + '_' + str(b['Percent Identity']) + '\n')
                f.write(b['Sequence'].strip() + '\n')

def write_contig_hits(contig_hits_dict):
    for k,v in contig_hits_dict.items():
        with open("AMR_STATISTICS/" + k + "/query_hits.txt", 'w') as f:
            for contig, hit_info in v.items():
                for hit, perc_id in hit_info.items():
                    f.write("{}\t{}\t{}\n".format(contig, hit, perc_id))


def top_hits_generator(best_hits):
    top_hits = defaultdict(lambda: defaultdict(dict))
    for k,v in best_hits.items():
        for a,b in v.items():
            if float(str(b['Percent Identity']).replace("*", "")) < 100.0:
                b['Percent Identity'] = '[' + str(b['Percent Identity']) + ']'
            top_hits[k][a] = b['Percent Identity']
    return top_hits

def create_mapping_dict(mapping_file):
    mapping_df = pd.read_csv(mapping_file, sep = '\t', index_col = 0).fillna("")
    mapping_dict = pd.DataFrame.to_dict(mapping_df, orient = 'index')
    return mapping_dict

def map_gene_info(row, mapping_dict):
    main_class = mapping_dict[row.name]['Main Class']
    row['Main Class'] = main_class
    sub_class = mapping_dict[row.name]['Sub Class']
    row['Sub Class'] = sub_class
    type = mapping_dict[row.name]['Type']
    row["Type"] = type
    curated = mapping_dict[row.name]['Curated']
    if curated == "":
        curated = row.name
    row["Curated"] = curated
    classification = mapping_dict[row.name]['Classification']
    row["Classification"] = classification
    return row

def create_top_hits_file(top_hits_dict, mapping_dict):
    top_hits_df = pd.DataFrame.from_dict(top_hits_dict).fillna("")
    top_hits_df.index.name = "Gene"
    pre_column_names = top_hits_df.columns.values
    after_column_names = ["Curated", "Main Class", "Sub Class", "Type", "Classification"]
    for i in pre_column_names:
        after_column_names.append(i)
    out_df = top_hits_df.apply(lambda x: map_gene_info(x, mapping_dict), axis =1)
    out_df = out_df[after_column_names]
    out_df.to_csv("AMR_STATISTICS/Top_AMR_Hits.txt", sep = '\t')
    return out_df


def arg_status(row, genomes):
    count = 0
    hit_genomes = []
    for i in genomes:
        if row[i] != "":
            count += 1
            hit_genomes.append(i)
    arg_status = ""
    singleton_genome = ""
    if count == len(genomes):
        arg_status = "Core"
    elif count == 1:
        arg_status = "Singleton"
        singleton_genome = hit_genomes[0]
    else:
        arg_status = "Shared"
    return pd.Series({"ARG Status": arg_status, "Count": count, "Singleton Genome": singleton_genome})

def create_args_file(top_hits_df, genomes):
    arg_stats_df = top_hits_df.copy()
    arg_stats_df[["ARG Status", "Count", "Singleton Genome"]] = arg_stats_df.apply(lambda x: arg_status(x, genomes), axis = 1)
    arg_stats_df = arg_stats_df.drop(genomes, axis = 1)
    arg_stats_df.columns = ["Curated", "Main Class", "Sub Class", "Classification", "Type", "ARG Status", "ARG Count", "Singleton Genome"]
    arg_stats_df.to_csv("AMR_STATISTICS/ARG_Stats.txt", sep = '\t')

def create_genome_stats_file(top_hits_df, genomes):
    genome_stats_df = top_hits_df.copy()
    counts_series = genome_stats_df.astype(bool).sum(axis = 0)[5:]

    total_counts_df = counts_series.to_frame()
    total_counts_df.index.name = "Genome"
    total_counts_df.columns = ['Total Count']

    main_class_values = genome_stats_df['Main Class'].unique()
    main_class_df = pd.DataFrame(counts_series.to_frame())
    for i in main_class_values:
        main_class_df[i] = genome_stats_df.loc[genome_stats_df['Main Class'] == i].astype(bool).sum(axis = 0)[5:].to_frame()
    main_class_df.columns.values[0] = 'Main Class Total'

    sub_class_values = genome_stats_df['Sub Class'].unique()
    sub_class_df = pd.DataFrame(counts_series.to_frame())
    for i in sub_class_values:
        sub_class_df[i] = genome_stats_df.loc[genome_stats_df['Sub Class'] == i].astype(bool).sum(axis = 0)[5:].to_frame()
    sub_class_df.columns.values[0] = 'Sub Class Total'
    column_headers = []
    for i in sub_class_df.columns:
        if i == "":
            i = "N/A"
        else:
            i = i
        column_headers.append(i)
    sub_class_df.columns = column_headers

    writer = pd.ExcelWriter("AMR_STATISTICS/Genome_Stats.xlsx")
    total_counts_df.to_excel(writer, sheet_name = "Total")
    main_class_df.to_excel(writer, sheet_name = "Main Class")
    sub_class_df.to_excel(writer, sheet_name = "Sub Class")
    writer.save()

def missing_genes(top_hits_df):
    missing = []
    for idx, row in top_hits_df.iterrows():
        try:
            if math.isnan(row['Main Class']):
                missing.append(idx)
        except:
            pass
    with open("AMR_STATISTICS/Missing_Genes.txt", 'w') as f:
        for i in missing:
            f.write(i + "\n")

def map_names_to_cells(row, genomes):
    for i in genomes:
        if row[i] == 100 or row[i] == '100.0*':
            row[i] = row['Curated']
        elif row[i] != "":
            row[i] = "[" + row["Curated"] + "]"
    return row


def combined_top_hits(top_hits_df, genomes):
    combined_hits_df = top_hits_df.copy()
    for idx, row in combined_hits_df.iterrows():
        if (row['Main Class'] == 'Beta-lactam' or row["Main Class"] == 'Aminoglycoside') and 'ef-tu' not in row['Curated'].lower():
            combined_hits_df.loc[idx, "Parent Allele"] = row['Curated'].rsplit('-', 1)[0]
        else:
            combined_hits_df.loc[idx, 'Parent Allele'] = row['Curated']
    combined_hits_df = combined_hits_df.drop(['Sub Class', "Classification", 'Type'], axis = 1)
    combined_hits_df = combined_hits_df.reset_index(drop = True)
    combined_hits_df = combined_hits_df.apply(lambda x: map_names_to_cells(x, genomes), axis = 1)
    combined_hits_df = combined_hits_df.drop(['Curated'], axis = 1)
    combined_hits_df = combined_hits_df.groupby(["Parent Allele", "Main Class"]).agg("".join)
    combined_hits_df[combined_hits_df.columns] = combined_hits_df.apply(lambda x: x.str.strip())
    combined_hits_df = combined_hits_df.transpose()
    combined_hits_df.to_csv("AMR_STATISTICS/Combined_Top_Hits.txt", sep = '\t', header = 0)

def main():
    parser = argparse.ArgumentParser(description = "Parse Resistance Genome Identifier (RGI) Output ")
    parser.add_argument('-t', '--target_genomes', help = 'File List (ls *.json -- except card.json file)', required = 'True')
    parser.add_argument('-a', '--all_hits', help = 'Generate an all hits file in addition to best hits file: Default is False', action = 'store_true')
    parser.add_argument('-m', '--mapping_file', help = 'Drug Class Mapping Term', default = '/usr/local/devel/ANNOTATION/cgreco/CARD/term_drug_class_mapping.txt')
    args = parser.parse_args()

    target_genomes = args.target_genomes
    all_hits = args.all_hits
    mapping_file = args.mapping_file

    directory = 'AMR_STATISTICS/'
    make_dir(directory)
    genomes_included = get_json_files(target_genomes)
    all_hits, contig_hits = generate_all_hits(genomes_included)
    best_hits = generate_best_hits(all_hits)
    write_all_hits(all_hits)
    write_best_hits(best_hits)
    write_contig_hits(contig_hits)

    top_hits_dict = top_hits_generator(best_hits)
    mapping_dict = create_mapping_dict(mapping_file)
    top_hits_df = create_top_hits_file(top_hits_dict, mapping_dict)
    genomes = list(top_hits_df.columns.values)[5:]

    create_args_file(top_hits_df, genomes)
    create_genome_stats_file(top_hits_df, genomes)
    missing_genes(top_hits_df)
    combined_top_hits(top_hits_df, genomes)

if __name__ == '__main__':
    main()
