import re
import pandas as pd
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sample', help = 'The name of the sample, e.g. Br61')
parser.add_argument('-i', '--input', help = 'The name of the tree in graphviz format')
parser.add_argument('--strict',action='store_true', help = 'Apply strict filtering of variants')
args = parser.parse_args()


print('Loading the graphviz-tree:')
with open(Path(f'/Users/jgawron/Documents/projects/CTC_backup/input_folder/{args.sample}/{args.input}'), 'r') as file:
    read_data = file.readlines()



pattern = re.compile(r'([0-9]{1,2})\[(color=gray93|color=ghostwhite)')

print('Removing single tumor cell and single WBCs from  the tree')
matches = []
remove = []
for it,line in enumerate(read_data):
    match = pattern.search(line)
    if match:
        matches.append(match.group(1))
        remove.append(it)

for it,match in enumerate(matches):
    stringToSearch = fr'([0-9]*) -> {match};'
    pattern2 = re.compile(stringToSearch)
    for it,line in enumerate(read_data):
        match2 = pattern2.search(line)
        if(match2):
            remove.append(it)

filtered_data = [read_data[it] for it in range(len(read_data)) if it not in remove]

print('Removing low impact variants from the tree')
### Removing variants from the tree that are either annotated as MODERATE or NONE
annotations = pd.read_csv(Path(f'/Users/jgawron/Documents/projects/CTC_backup/input_folder/{args.sample}/{args.sample}_variants_annotations.csv'))

if(args.strict == True):
    filter_variants = annotations.loc[(annotations['relevant'] == 'NONE') | (annotations['relevant'] == 'MODERATE'), 'variantName']
else:
    filter_variants = annotations.loc[(annotations['relevant'] == 'NONE') , 'variantName']

filtered_data2 = []
for var in filter_variants:
    filtered_data2 = []
    pattern = fr'\b{var}\w*\b'
    for line in filtered_data:
        filtered_data2.append(re.sub(pattern, '', line))
    filtered_data = filtered_data2


filtered_data2 = [line for line in filtered_data if line.strip()]
filtered_data = []

for line in filtered_data2:
    if line.strip().endswith('label="'):
        filtered_data.append(line.rstrip('\n'))
    else:
        filtered_data.append(line)

print('Removing intermediate empty nodes')
###Removing intermediate empty nodes
### First, as a cleanup, make sure that each line of the file corresponds to exactly one entry of the vector.

full_string = ''.join(filtered_data)
filtered_data = full_string.splitlines()

filtered_data = [line + '\n' for line in filtered_data]

###Next, identify the node number of the root:
#### It is the one specified in the line right before the following: 
pattern = re.compile(r'node \[fontname=helvetica')
for it, line in enumerate(filtered_data):
    matching = pattern.search(line)
    if matching:
        pattern = re.compile(r'([0-9]*) \[')
        while True:
            it = it - 1
            matching = pattern.search(filtered_data[it])
            if matching:
                rootNode = matching.group(1)
                break
        break


matches = []
remove = []
pattern = re.compile(r'([0-9]*) \[label=""\];')

### Now identify and remove all empty nodes
leave = 0
while True:
    for it,line in enumerate(filtered_data):
        matching = pattern.search(line)
        if matching:
            match = matching.group(1)
            if match == rootNode:
                continue
            remove=[it]
            stringToSearch = fr'([0-9]*) -> {match};'
            pattern2 = re.compile(stringToSearch)
            stringToReplace = fr'{match} -> ([0-9]*);'
            pattern3 = re.compile(stringToReplace)


            for it2,line2 in enumerate(filtered_data):
                match2 = pattern2.search(line2)
                if match2:
                    ancestor = match2.group(1)
                    break
            remove.append(it2)

            for it3, line3 in enumerate(filtered_data):
                match3 = pattern3.search(line3)
                if match3:
                    remove.append(it3)
                    filtered_data.append(str(ancestor) + ' -> ' + match3.group(1) + ';\n')



            for idx in sorted(remove, reverse=True):
                del filtered_data[idx]
            break
        if it == len(filtered_data)-1:
            leave = 1
    if leave == 1:
        break    


        
 #       matches.append(match.group(1))
 #       remove.append(it)

#print(matches)
#print(remove)


#for it,match in enumerate(matches) :
    

#filtered_data2 = [line for line in filtered_data if line.strip()]
filtered_data2 = filtered_data
for it,line in enumerate(filtered_data2):
    if line == '}\n':
        del filtered_data2[it]

filtered_data2.append('}')

print('Removing nodes without descending cells')
### Next we need to remove nodes that have no descendants
remove = []
pattern = re.compile(r'([0-9]*) \[label;')

for it,line in enumerate(filtered_data):
    matching = pattern.search(line)
    if matching:
        match = matching.group(1)
        stringToSearch = fr'{match} -> ([0-9]*);'
        pattern2 = re.compile(stringToSearch)
        any_matching = 0
        for it2,line2 in enumerate(filtered_data):
                match2 = pattern2.search(line2)
                if match2:
                    any_matching = 1
        if any_matching == 0:
            remove.append(it)

for idx in sorted(remove, reverse=True):
    del filtered_data2[idx]


print('Removing WBCs that are part of cell clusters')
### Next, we remove all WBCs that are part of cell clusters. We first load a cell description file and select the relevant cells based on this.

cellDescription = pd.read_table(Path(f'/Users/jgawron/Documents/projects/CTC_backup/input_folder/{args.sample}/{args.sample}_cellDescription.tsv') , sep = '\t')
cells = [value for value in cellDescription['color'].unique() if value not in ['gray93', 'ghostwhite']]


for color in cells:
    print(color)
    pattern = re.compile(rf'([0-9]*)\[color={color}')
    occurrences = sum(cellDescription.loc[cellDescription['color'] == color,'WBC'] == 1)
    
    for i in range(occurrences):
        print(i)
        remove = []
        for it,line in enumerate(filtered_data):
            matching = pattern.search(line)
            if matching:
                match = matching.group(1)
                remove.append(it)
                pattern2 = re.compile(rf'([0-9]*) -> {match};')
                for it2, line2 in enumerate(filtered_data):
                    matching2 = pattern2.search(line2)
                    if matching2:
                        remove.append(it2)
                break
        print(remove)
        for idx in sorted(remove, reverse=True):
            del filtered_data2[idx]

print(filtered_data2)
outpath = Path(f'/Users/jgawron/Documents/projects/CTC_backup/input_folder/{args.sample}/{args.input}')
if args.strict == True:
    outfile = f'{outpath.stem}_modified_strict.gv'
else:
    outfile = f'{outpath.stem}_modified.gv'

with open(outpath.parents[0] / outfile, 'w') as file:
    file.writelines(line for line in filtered_data2)
