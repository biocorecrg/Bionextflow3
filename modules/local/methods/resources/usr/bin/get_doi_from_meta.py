#!/usr/bin/env python

import argparse
import glob
import os
import yaml

parser = argparse.ArgumentParser(
    description="Retrieve information from meta.yml of each module for methods description.",
    epilog="Example usage: python methods.py -p /path/to/repo",
    add_help=True,
)
parser.add_argument("-p", "--path", help="Path to the repository", required=True, type=os.path.abspath)
args = parser.parse_args()

path = args.path


## Find meta.yaml files in a given path and stores in a ditionary dirpath


def find_meta_yaml_files(path):
    """find meta.yaml recursively in a given directory"""
    meta_yaml_files = []
    path = os.path.abspath(path) 
    for dirpath, dirnames, filenames in os.walk(path, followlinks=True):
        if "meta.yml" in filenames:
            meta_yaml_files.append(
                os.path.join(dirpath, filenames[filenames.index("meta.yml")])
            )
    return meta_yaml_files
    


def get_tool_info_from_meta(meta_yaml):
    """Reading a given yaml file, and retrieving the doi and homepage of each tool, returning a list of dictionaries for each tool found"""
    tool_list = list()
    try:
        with open(meta_yaml, "r") as file:
            meta = yaml.safe_load(file)
            if "tools" in meta:
                tools = meta["tools"][0]
                for k, v in tools.items():
                    doi = v.get("doi", "")
                    homepage = v.get("homepage", "")
                    tool_list.append({k: {"doi": doi, "homepage": homepage}})

        return tool_list
    except yaml.scanner.ScannerError:
        return ["Error reading {meta_yaml} file, please check your yaml syntax"]


def read_mqc_versions(path):
    """Reading the nf_core_pipeline_software_mqc_versions.yml file, returning a dictionary with the tool name and version"""
    mqc_versions = glob.glob("**/*mqc_versions.yml", recursive=True)
    used_tools_ver = {}
    
    if mqc_versions:
        with open(mqc_versions[0], "r") as file:
            mqc_data = yaml.safe_load(file)
            workflow = mqc_data["Workflow"]
            for values in mqc_data.values():
                used_tools_ver.update(values)

        return [used_tools_ver, workflow]
    else:
        return ["No mqc_versions.yml file found, please check your pipeline path contains one"]


def match_used_tools(path):
    """Matching the used tools found in mqc_versions and the tools found in meta.yaml, to create a final dictionary with doi, homepage and version"""


    ## Tools and versions used in the pipeline
    used_tools_ver = read_mqc_versions(path)[0]
    workflow = read_mqc_versions(path)[1]

    meta_files = find_meta_yaml_files(path)
    if len(meta_files) == 0:
        return "No meta.yaml files found, please check the path"
    else:

        final_tools = {}
        for meta_file in meta_files:
            # print(f"\nProcessing {meta_file}")
            tools_info = get_tool_info_from_meta(meta_file)
            print(f"Tools info from {meta_file}: {tools_info}")
            for tool in tools_info:
                # print(f"Tool info: {tool}")
                try:
                    common_keys = tool.keys() & used_tools_ver.keys()
                    joined = {
                        k: {**tool[k], "version": used_tools_ver[k]} for k in common_keys
                    }
                    final_tools.update(joined)
                except AttributeError:
                    pass

        # print(f"\nThe tools used in mqc_version.yml with meta.yml information are:\n{final_tools}")
        tools_and_workflow = {"workflow": workflow, "tools": final_tools}

        return tools_and_workflow

    # mqc_versions = read_mqc_versions(path)
    # print(mqc_versions)


# find_meta_yaml_files(path)
# dir1_meta_files = find_meta_yaml_files(path)
# # print(dir1_meta_files[1])

# #print(f"{list(dir1_meta_files.values())[0]}")

# # print(get_tool_info_from_meta(dir1_meta_files[1]))

# tools_info = [get_tool_info_from_meta(f) for f in dir1_meta_files]
# print(tools_info)


# ## Matching tools used with the total found in the repository

with open("dois.yml", "w") as outfile:
    yaml.dump(match_used_tools(path), outfile, default_flow_style=False)