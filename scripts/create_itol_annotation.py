import sys, os, argparse, re
import pandas as pd
from collections import defaultdict, Counter, OrderedDict
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

pd.options.display.max_rows = 999
parser = argparse.ArgumentParser(description = "Create ITOL Annotation based upon CARD RGI output")
parser.add_argument('-t', '--top_hits_file', help = 'Top AMR Hits File (in AMR_STATISTICS directory)', required = 'True')
parser.add_argument("-a", '--annotation_type', help = "Which type of ITOL annotation file would you like generated? Options are 'Bar' and 'Binary'", default = 'Bar')
parser.add_argument("-m", '--main_class', help = "Choose one Class to plot (e.g. Beta-lactam, Aminoglycoside, etc. Case Sensitive")
parser.add_argument("-o", '--out_file', help = 'Name of the file you want Annotations written to. Default is ITOL_OUT.txt', default = 'ITOL_OUT.txt')
parser.add_argument("-r", '--remove_multidrug_genes', help = 'Remove genes that confer non specific resistance (Efflux Pumps, etc.)', action = 'store_true')
parser.add_argument("-c", '--carbapenemase_highlight', help = 'If main class is Beta-lactam, then this flag highlights those AMR genes that confer carbapenem resistance', action = 'store_true')
parser.add_argument("-g", '--group_alleles', help = 'Group all alleles of a gene together (With -a of Bar-- KPC-2 and KPC-3 become KPC with a count of 2; With -a of Binary-- KPC-2 and KPC-3 becomes KPC with a shading of the highest %%id', action = 'store_true')
args = parser.parse_args()
top_amr_hits = args.top_hits_file
annotation_type = args.annotation_type.lower()
class_subset = args.main_class
out_file = args.out_file
rm_multidrug = args.remove_multidrug_genes
carb_highlight_yes = args.carbapenemase_highlight
group_alleles = args.group_alleles
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

df = pd.read_csv(top_amr_hits, sep = '\t', header = 0, index_col = "Gene", low_memory = False)

amr_dict = df.to_dict()
main_class_dict = amr_dict['Main Class']

color_palette = ["#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921",
				"#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E",
				"#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B",
				"#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0",
				"#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861"]



def create_class_counts_not_grouped(top_hits_dict, class_mapping_dict, multidrug_option):
	genome_to_class  = defaultdict(list)
	for k,v in top_hits_dict.items():
		if k != 'Main Class' and k != 'Sub Class' and k != 'Type' and k != 'Classification' and k != 'Curated Name':
			if multidrug_option:
				for a,b in v.items():
					if b == '0':
						pass
					else:
						category = str(class_mapping_dict[a])
						if category.lower().startswith('multi'):
							pass
						else:
							genome_to_class[k].append(category)
			else:
				for a,b in v.items():
					if b == '0':
						pass
					else:
						genome_to_class[k].append(class_mapping_dict[a])

	out_dict = {}
	for k, v in genome_to_class.items():
		cnt = Counter()
		for i in v:
			cnt[i] += 1
		out_dict[k] = cnt

	return out_dict

def create_color_dictionary(top_hits_dict, class_mapping_dict, color_palette_list):
	color_dictionary = {}

	main_class_set = set()
	for k,v in top_hits_dict.items():
		for a,b in v.items():
			main_class_set.add(a)

	main_class_list = list(main_class_set)

	for i in range(0, len(main_class_list)):
		color_dictionary[main_class_list[i]] = color_palette[i]

	return color_dictionary

def create_sorted_counts(class_counts, color_mapping):
	amr_classes = color_dictionary.keys()

	for k,v in class_counts.items():
		present_classes = []
		for a,b in v.items():
			present_classes.append(a)

		for i in amr_classes:
			if i not in present_classes:
				class_counts[k][i] = 0

	out_dict = {}
	for k,v in class_counts.items():
		out_dict[k] = OrderedDict(sorted(v.items(), key=lambda t: t[0]))

	return out_dict

def create_annotation_parts(color_mapping, sorted_counts_dict):
	amr_classes = sorted(color_mapping.keys())
	amr_colors = sorted(color_mapping.values())

	try:
		field_labels = 'FIELD_LABELS' + "," + ','.join(amr_classes)
		field_colors = "FIELD_COLORS" + "," + ','.join(amr_colors)
	except:
		print "Check for any missing genes. These missing AMR Classes cause the script to fail."

	annotation_list = []

	for k,v in sorted_counts_dict.items():
		sorted_counts = []
		for a,b in v.items():
			sorted_counts.append(str(b))
		annotation_list.append(k + "," + ",".join(sorted_counts))

	return field_labels, field_colors, annotation_list



def write_multibar_file(field_labels, field_colors, amr_counts, out_file):
	legend_labels = ",".join(field_labels.split(",")[1:])
	legend_colors = ",".join(field_colors.split(",")[1:])
	legend_shapes = "LEGEND_SHAPES"
	for i in range(0, len(legend_labels.split(","))):
		legend_shapes = legend_shapes + ",1"
	with open(out_file, 'w') as f:
		f.write("DATASET_MULTIBAR" + "\n")
		f.write("SEPARATOR COMMA"+ "\n")
		f.write("DATASET_LABEL,AMR Bar Chart"+ "\n")
		f.write("LEGEND_TITLE,AMR Classes" + "\n")
		f.write(legend_shapes + "\n")
		f.write("LEGEND_COLORS," + legend_colors + "\n")
		f.write("LEGEND_LABELS," + legend_labels + '\n')
		f.write("COLOR,#ff0000"+ "\n")
		f.write("HEIGHT_FACTOR,2" + "\n")
		f.write(field_labels+ "\n")
		f.write(field_colors+ "\n")
		f.write("ALIGN_FIELDS,1"+ "\n")
		f.write("DATA"+ "\n")
		for i in amr_counts:
			f.write(i+ "\n")


def subset_multibar_to_single_class(field_labels, field_colors, amr_counts, class_of_interest):
	label_list = field_labels.split(",")
	for i in range(0, len(label_list)):
		if label_list[i] == class_of_interest:
			index = i

	single_field_labels = field_labels.split(",")[0] + "," + field_labels.split(",")[index]
	single_field_colors = field_colors.split(",")[0] + "," + field_colors.split(",")[index]
	new_amr_counts = []
	for i in amr_counts:
		new_amr_counts.append(i.split(",")[0] + "," + i.split(",")[index])

	return single_field_labels, single_field_colors, new_amr_counts

def write_singlebar_file(single_label, single_color, single_counts, out_file, subset_class):
	legend_label = single_label.split(",")[1]
	legend_color = single_color.split(",")[1]
	legend_shape = '1'
	with open(out_file, 'w') as f:
		f.write("DATASET_SIMPLEBAR" + "\n")
		f.write("SEPARATOR COMMA"+ "\n")
		f.write("DATASET_LABEL," + subset_class + "\n")
		f.write("LEGEND_TITLE," + subset_class + "\n")
		f.write("LEGEND_SHAPES," + legend_shape + "\n")
		f.write("LEGEND_COLORS," + legend_color + "\n")
		f.write("LEGEND_LABELS,Count of " + legend_label + "\n")
		f.write("COLOR,#00b2ff"+ "\n")
		f.write("HEIGHT_FACTOR,2" + "\n")
		f.write(single_label+ "\n")
		f.write(single_color+ "\n")
		f.write("DATA"+ "\n")
		for i in single_counts:
			f.write(i+ "\n")

#----------------------------------------------------------------------------------------------------------------#
#Binary Dataset

def select_main_class_perc_ids(amr_dict, class_subset, main_class_dict, group_alleles):
	subset_dict = defaultdict(lambda: defaultdict(list))
	for k,v in amr_dict.items():
		if k != 'Main Class' and k != 'Sub Class' and k != 'Type' and k != 'Classification' and k != 'Curated Name':
			for a,b in v.items():
				if main_class_dict[a] == class_subset:
					if group_alleles:
						gene_name = a.rsplit("-", 1)[0]
						subset_dict[k][gene_name].append(b)
					else:
						subset_dict[k][a].append(b)

	node_dict = defaultdict(dict)

	for k,v in subset_dict.items():
		for a,b in v.items():
			if len(b) == 1:
				gene_value = b[0]
			else:
				gene_value = max(b)
			if gene_value == '0':
				node_dict[k][a] = -1
			elif gene_value == '100.0' or gene_value == '100' or gene_value == '100*':
				node_dict[k][a] = 1
			else:
				node_dict[k][a] = 0

	out_dict = {}
	for k,v in node_dict.items():
		out_dict[k] = OrderedDict(sorted(v.items(), key=lambda t: t[0]))
	return out_dict

def binary_annotations(color_palette, binary_dict):
	binary_colors = []
	for k,v in binary_dict.items():
		for i in range(0, len(v)):
			if i < len(color_palette):
				binary_colors.append(color_palette[i])
			else:
				try:
					binary_colors.append(color_palette[i - len(color_palette)])
				except:
					color = i % len(color_palette)
					binary_colors.append(color_palette[color])
		break

	binary_color = "FIELD_COLORS," + ",".join(binary_colors)

	binary_shapes = []
	for k,v in binary_dict.items():
		for i in range(1, len(v)+1):
			if i < 7:
				binary_shapes.append(str(i))
			else:
				binary_shapes.append(str(i % 6 + 1))
		break

	binary_shape = "FIELD_SHAPES," + ",".join(binary_shapes)

	labels = []
	for k,v in binary_dict.items():
		for a,b in v.items():
			labels.append(a)
		break

	binary_counts = defaultdict(list)
	for k,v in binary_dict.items():
		value = ''
		for a,b in v.items():
			value = value + "," + str(b)
		binary_counts[k] = value
	binary_labels = "FIELD_LABELS," + ",".join(labels)

	return binary_color, binary_shape, binary_labels, binary_counts

def carbapenemase_highlighting(shapes, labels, amr_dict):
	split_labels =  labels.split(",")[1:]
	split_shapes = shapes.split(",")[1:]
	count = 1
	for i in range(0, len(split_labels)):
		sub_class = amr_dict['Sub Class'][split_labels[i]]
		if isinstance(sub_class, basestring):
			if 'carbapenemase' in sub_class.lower():
				split_shapes[i] = '1'
			else:
				split_shapes[i] = '2'
		else:
			split_shapes[i] = '2'

		count = count + 1


	carb_shapes = "FIELD_SHAPES," + ','.join(split_shapes)
	carb_labels = "FIELD_LABELS," + ','.join(split_labels)

	return carb_shapes, carb_labels


def write_binary_file(colors, shapes, labels, counts, out_file, subset_class):
	legend_shapes = ",".join(shapes.split(",")[1:])
	legend_colors = ",".join(colors.split(",")[1:])
	legend_labels = ",".join(labels.split(",")[1:])
	with open(out_file, 'w') as f:
		f.write("DATASET_BINARY" + "\n")
		f.write("SEPARATOR COMMA"+ "\n")
		f.write("DATASET_LABEL," + subset_class + "\n")
		f.write("LEGEND_TITLE," + subset_class + "\n")
		f.write("LEGEND_SHAPES," + legend_shapes + "\n")
		f.write("LEGEND_COLORS," + legend_colors + "\n")
		f.write("LEGEND_LABELS," + legend_labels + "\n")
		f.write("COLOR,#00b2ff"+ "\n")
		f.write(shapes+ "\n")
		f.write(labels+ "\n")
		f.write(colors+ "\n")
		f.write("DATA"+ "\n")
		for k,v in counts.items():
			f.write(k + v + "\n")


class_counts_dict = create_class_counts_not_grouped(amr_dict, main_class_dict, rm_multidrug)
color_dictionary = create_color_dictionary(class_counts_dict, main_class_dict, color_palette)
sorted_dict = create_sorted_counts(class_counts_dict, color_dictionary)
labels, colors, annotation_counts = create_annotation_parts(color_dictionary, sorted_dict)
binary_dict = select_main_class_perc_ids(amr_dict, class_subset, main_class_dict, group_alleles)
binary_color, binary_shape, binary_label, binary_count = binary_annotations(color_palette, binary_dict)

if annotation_type.lower() == 'binary':
	if class_subset:
		if class_subset not in color_dictionary.keys():
			print "Your chosen Main Class (-m option) is not valid. Here are the available choices: "
			for i in sorted(color_dictionary.keys()):
				print i
			sys.exit()
		elif class_subset == 'Beta-lactam' and carb_highlight_yes:
			carb_binary_shapes, carb_binary_labels = carbapenemase_highlighting(binary_shape, binary_label, amr_dict)
			write_binary_file(binary_color, carb_binary_shapes, carb_binary_labels, binary_count, out_file, class_subset)
		elif class_subset != 'Beta-lactam' and carb_highlight_yes:
			print
			print "ERROR: You can't highlight Carbapenems without choosing Beta-lactam as the main class"
			print
			sys.exit()
		else:
			write_binary_file(binary_color, binary_shape, binary_label, binary_count, out_file, class_subset)
	else:
		print
		print "ERROR: You must specify a main class to create a binary type plot"
		print
		sys.exit()
else:
	if class_subset:
		single_labels, single_colors, single_amr_counts = subset_multibar_to_single_class(labels, colors, annotation_counts, class_subset)
		write_singlebar_file(single_labels, single_colors, single_amr_counts, out_file, class_subset)
	elif carb_highlight_yes:
		print
		print "ERROR: You can't highlight Carbapenems in a bar graph"
		print
		sys.exit()
	else:
		write_multibar_file(labels, colors, annotation_counts, out_file)
