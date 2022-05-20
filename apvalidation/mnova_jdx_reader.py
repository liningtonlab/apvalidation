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

def find_maxdepth_groups(line_level_list, max_depth):
    groups = []
    agroup = []
    grouping = False
    for line in line_level_list:
        if line["level"] == max_depth or line['level'] == max_depth-1:
            line["value"] = line["value"].replace("\t", "")
            agroup.append(line)
            grouping = True
        elif grouping == True and (line["level"] != max_depth or line["level"] != max_depth-1):
            grouping = False
            groups.append(agroup)
            agroup = []
        
    return groups

def save_separate_files(max_depth_groups, save_path):

    if not os.path.exists(save_path):
        os.makedirs(save_path)
    

    for group in max_depth_groups:
        filename = group[0]["value"]
        


        

def separate_mnova_jdx(filepath):
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
    max_depth_groups = find_maxdepth_groups(line_level_list, max_depth)
    print(max_depth_groups)
    return max_depth_groups

    

    
























# def read_mnova_jdx(filepath):
#     def recursive_function(line_list, index):
#         curr_dict = {"values": [], "nested_dicts": []}
#         # TRY USING A WHILE LOOP THAT ALLOWS THE INDEX VALUE TO CHANGE
#         print(f"length of list: {len(line_list)}")
#         try:
#             while index < len(line_list[index:]):
#                 curr_dict['values'].append(line_list[index])
#                 if line_list[index].startswith("##END="):
#                     # print(f"RETURN {line_list[index]}")
#                     return curr_dict, index
#                 if line_list[index].startswith("##TITLE="):
#                     # print(f"CALL {line_list[index]}")
#                     nested_dict, output_index = recursive_function(line_list, index+1)
#                     curr_dict["nested_dicts"].append(nested_dict)
#                     index = output_index
#                 print(index)
#                 index+=1
#         except TypeError:
#             print(curr_dict)

#         # for index, line in enumerate(line_list):
#         #     curr_dict['values'].append(line)
#         #     if line.startswith("##END="):
#         #         print(f"RETURN {line}")
#         #         return curr_dict
#         #     if line.startswith("##TITLE="):
#         #         print(f"CALL {line}")
#         #         nested_dict = recursive_function(line_list[index+1:])
#         #         curr_dict["nested_dicts"].append(nested_dict)

#     def split_on_items_function(line_list):
#         get_next_title = False
#         get_lines = False
#         for line in line_list:
#             if line.startswith("##$MNOVA/LINK_BLOCK_TYPE=	ITEM	$$ NMR Spectrum"):
#                 get_next_title = True

#             if line.startswith("##TITLE=") and get_next_title==True:
#                 get_lines = True
#                 get_next_title = False

#             if line.startswith("##END="):
#                 get_lines = False

#             if get_lines:
#                 pass


                
        

#     with open(f"{filepath}", "r") as f:
#         line_list = f.readlines()
#         # print(len(line_list))
#         output, index = recursive_function(line_list, 0)
#         return output
        

#     # with open(f"{filepath}", "r") as f:
#     #     line_list = f.readlines()
#     #     tabs=0
#     #     with open("test_files/MNOVA_jdx/WriteFiles/tab_jdx.txt", "w") as fi:
#     #         nested_dict = {}
#     #         depth = 0
#     #         for line in line_list:
#     #             if line.startswith("##TITLE="):
#     #                 curr_key = 'temp_key'
#     #                 nested_dict[curr_key] = {'values': [], 'nested_lines': {}}
#     #             if line.startswith("##$MNOVA/LINK_BLOCK_TYPE="):
#     #                 new_key = line.split("=")[-1]
#     #                 nested_dict[new_key] = nested_dict.pop(curr_key)
#     #                 curr_key = new_key
#     #             nested_dict[curr_key]['values'].append(line)
                    




















# # def read_varian_combined(filepath):
# #     with open(filepath, 'r') as file:
# #         line_list = file.readlines()

# #     fid_list = []
# #     procpar_list = []
# #     for line in line_list:
# #         if line.startswith("##$PARAMETER FILE=") and "procpar" in line:
# #             procpar_list.append(line)
# #         if line.startswith("##$PARAMETER FILE=") and "text" in line:
# #             fid_list.append(line)

# #     experiment_procpars = {}
# #     experiment_fids = {}
# #     current_fid_name = None
# #     current_fid_lines = []
# #     current_procpar_name = None
# #     current_procpar_lines = []

# #     for line in line_list:
# #         if line in procpar_list:
# #             current_procpar_name = line
# #         elif current_procpar_name is not None:
# #             current_procpar_lines.append(line.replace("\n", ""))

# #         if line in fid_list:
# #             current_fid_name = line

# #             experiment_procpars[current_procpar_name] = current_procpar_lines
# #             current_procpar_name = None
# #             current_procpar_lines = []
# #         elif current_fid_name is not None:
# #             current_fid_lines.append(line.replace("\n", ""))

# #         if "##END=" in line and current_fid_name is not None:
# #             experiment_fids[current_fid_name] = current_fid_lines
# #             current_fid_name = None
# #             current_fid_lines = []

    
# #     complete_experiment_list = {}
# #     for i in range(len(experiment_procpars)):
# #         key_procpar = list(experiment_procpars.keys())[i]
# #         key_fid = list(experiment_fids.keys())[i]
# #         title = key_procpar.split("/")[-2]
# #         complete_experiment_list[title] = [experiment_procpars[key_procpar], experiment_fids[key_fid]]
    
# #     # with open("test_files/read_varian_jdx_output", "w") as f:
# #     #     f.write(str(complete_experiment_list))
# #     return complete_experiment_list