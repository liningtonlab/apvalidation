import os, shutil

def add_indents(line_list):
    """
        add indents to a list of string lines read from a JDX file based on the ##TITLE and 
        ##END tags. The purpose is to define the nested layers of the JDX file to ensure
        easy parsing

        :param line_list: a list of string lines read from the JDX file in question.
        :return: a list of string lines with the specified number of tabs
    """

    tabs=0
    for index, line in enumerate(line_list):
        if line.startswith("##END="):
            tabs-=1
        if line.startswith("##TITLE="):
            tabs+=1
        line_list[index]=tabs*"\t"+line
    return line_list

def count_tabs(input_string):
    """
        count the number of tabs in a given line and return it.
        :param input_string: a string line with tabs
        :return: an int representing the amount of tabs
    """
    return len(input_string) - len(input_string.lstrip("\t"))

def find_deep_groups(line_level_list, max_depth):
    """
        extract the two deepest levels of the nested jdx file as groups. The result will
        be a list of groupings where each group is the largest set contiguous lines at the two deepest
        levels of the file. The purpose of this is to remove the meta data tagging that comes with the combined
        jdx file.
        :param line_level_list: a list of dict objects that contain the value for each line as well as the depth the line is
                                nested at.
        :param max_depth: an int representing the max nesting in the file
        :return: a list containing groupings of deep lines in the file. This list is separated by lines that were
                contiguous in the original file.
    """
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

def make_filename(item, group, group_num):
    """
        determine the name of the single experiment jdx file being saved.
        :param item: a single line 
    """

    for item in group:
        if item["value"].startswith("##TITLE="):
            group_title = item["value"].split("=")[1]
            group_title = group_title.split("/")[-1]
            group_title = f"{group_title}_exp_{group_num+1}.jdx"
            return group_title
        else:
            return "titlenotfound"

def save_separate_files(merged_groups, save_path):

    if not os.path.exists(save_path):
        os.makedirs(save_path)
    else:
        for filename in os.listdir(save_path):
            file_path = os.path.join(save_path, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))


    for index, group in enumerate(merged_groups):
        for item in group:
            group_title = make_filename(item, group, index)
        
        single_file = open(f"{save_path}/{group_title}", "w+")
        for index, item in enumerate(group):
            print(f"writing to {group_title} item #{index} being {item}....")
            single_file.write(item["value"]+"\n")
        single_file.close()
    return save_path

    
def separate_mnova_jdx(input_filepath, save_location):
    with open(f"{input_filepath}", "r") as f:
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
    
    saved_folder = save_separate_files(deep_groups, f"{save_location}")
    
    return saved_folder

