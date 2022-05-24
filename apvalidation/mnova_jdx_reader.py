import os

def add_indents(line_list):    
    tabs=0
    for index, line in enumerate(line_list):
        if line.startswith("##END="):
            tabs-=1
        if line.startswith("##TITLE="):
            tabs+=1
        line_list[index]=tabs*"\t"+line
    return line_list

def count_tabs(input_string):
    return len(input_string) - len(input_string.lstrip("\t"))

def find_deep_groups(line_level_list, max_depth):
    groups = []
    agroup = []
    grouping = False
    for line in line_level_list:
        if line["level"] == max_depth or line['level'] == max_depth-1:
            line["value"] = line["value"].replace("\t", "")
            line["value"] = line["value"].replace("\n", "")
            agroup.append(line)
            grouping = True
        elif grouping == True and (line["level"] != max_depth or line["level"] != max_depth-1):
            grouping = False
            groups.append(agroup)
            agroup = []
        
    return groups

def make_filename(item, group):
    for item in group:
        if item["value"].startswith("##TITLE="):
            group_title = item["value"].split("=")[1] + ".jdx"
            group_title = group_title.split("/")[-1]
            return group_title
        else:
            return "no_title_found"

def save_separate_files(merged_groups, save_path):

    if not os.path.exists(save_path):
        os.makedirs(save_path)
    
    for group in merged_groups:
        for item in group:
            group_title = make_filename(item, group)
        
        single_file = open(f"{save_path}/{group_title}", "w+")
        for index, item in enumerate(group):
            print(f"writing to {group_title} item #{index} being {item}....")
            single_file.write(item["value"]+"\n")
        single_file.close()
    
def separate_mnova_jdx(filepath, save_location):
    with open(f"{filepath}", "r") as f:
        line_list = f.readlines()
        line_list_tabs = add_indents(line_list)
    
    line_level_list = []
    max_depth = 0
    for line in line_list_tabs:
        depth = count_tabs(line)
        if max_depth < depth:
            max_depth = depth
        line_level_list.append({"value": line, "level": depth})

    deep_groups = find_deep_groups(line_level_list, max_depth)
    
    save_separate_files(deep_groups, f"{save_location}")
    print(len(deep_groups))
    print(max_depth)
    
    return deep_groups

