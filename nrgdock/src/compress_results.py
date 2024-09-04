import os
import shutil
import sys
from pathlib import Path
import numpy as np
from main_processed_target import get_params_dict
import concurrent.futures
from itertools import repeat, groupby
import argparse
# import psutil
import time
import polars as pl
# from query_db import main as query_db


def load_df(csv_path, ligand_type=None):
    df = pl.scan_csv(csv_path)
    result_df = df.sort('CF').group_by('Name').first().sort('CF')
    result_df = result_df.collect()
    if ligand_type is None:
        ligand_type = result_df.row(0, named=True)['Type']
        if ligand_type == 'ligand':
            need_for_ef = False
        else:
            need_for_ef = True
    else:
        need_for_ef = False
    return result_df, need_for_ef


def get_params(cfg_path):
    suffix = "no_conf"
    params_dict = get_params_dict(cfg_path)
    conf_num = params_dict["CONFORMER_NUMBER"]
    if conf_num != 0:
        suffix = f"{conf_num}_conf"
    return params_dict, suffix


def get_output_name(suffix, n_rotations, folder_end):
    base_output_path = f"./results_processed"
    if not os.path.isdir(base_output_path):
        os.mkdir(base_output_path)
    output_path = f"./results_processed/{n_rotations}_rotations_{suffix}"
    if folder_end is not None:
        output_path += f'_{folder_end}'
    if os.path.exists(output_path):
        number = 2
        while os.path.exists(f"{output_path}_{str(number)}/"):
            number += 1
        output_path = f"{output_path}_{str(number)}"
    os.mkdir(output_path)
    return output_path


def reset_result_folder(results, params_dict):
    target_names = next(os.walk(results), (None, None, []))[1]
    if params_dict["CLEAN"]:
        shutil.rmtree(results)
        os.mkdir("./temp/results")
        for target in target_names:
            output_path = os.path.join(results, target)
            if not os.path.exists(output_path):
                os.makedirs(output_path)


def get_good_ligands(number_pdb_to_keep, ligand_list):
    good_ligands = []
    if len(ligand_list) < number_pdb_to_keep + 1:
        number_pdb_to_keep = len(ligand_list) - 1
    for ligand_counter, ligand in enumerate(ligand_list):
        if ligand_counter > int(number_pdb_to_keep) - 1:
            return good_ligands
        else:
            try:
                good_ligands.append(ligand['Name']+"_" + ligand["Conformer_number"])
            except KeyError:
                good_ligands.append(ligand['Name'])


def write_remove_duplicate(target, remove_duplicate_path, df):
    categories_of_interest = ['Name', 'CF']
    write_path = os.path.join(remove_duplicate_path, target + ".csv")
    df.select(categories_of_interest).write_csv(write_path, separator=",")


def calculate_ef(df, active_total_counter, decoy_total_counter, ef_percentage):
    list_size_percent_decoy = round(df.height * ef_percentage)
    top_1_percent_df = df.head(list_size_percent_decoy)
    active_ligands = top_1_percent_df.filter(pl.col('Type') == 'active').height
    proportion = active_ligands/(df.height*ef_percentage)
    try:
        EF = str(round(proportion/(active_total_counter/(active_total_counter + decoy_total_counter)), 2))
    except ZeroDivisionError:
        EF = "Error"
    return EF


def delete_ligands(number, ligand_list, output_path, dir):
    number = int(number)
    good_ligands = get_good_ligands(number, ligand_list)
    base_ligand_output_path = os.path.join(output_path, "poses")
    if not os.path.isdir(base_ligand_output_path):
        try:
            os.mkdir(base_ligand_output_path)
        except:
            pass
    new_dir = os.path.join(base_ligand_output_path, dir + "_ligand_poses")
    init_dir = None
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
    for filename in good_ligands:
        init_dir = os.path.join(".", "temp", "ligand_poses", dir)
        try:
            shutil.copyfile(os.path.join(init_dir, filename), os.path.join(new_dir, filename))
        except:
            print("error copying ligand")
            continue
    try:
        shutil.rmtree(init_dir)
        os.mkdir(init_dir)
    except FileNotFoundError:
        print("Could not delete ligand folder")
        return


def save_ef_csv(list_of_dict, result_save_path):
    organised_results_sorted = sorted(list_of_dict, key=lambda x: x.get('Name'))
    with open(result_save_path, "w") as f:
        f.write('Target,EF\n')
        for result in organised_results_sorted:
            try:
                f.write(f'{result["Name"]},{float(result["EF"]):.2f}\n')
            except ValueError:
                f.write(f'{result["Name"]},{result["EF"]}\n')


def merge_csv(target_path, processed_result_path):
    current_info = os.path.join(target_path, 'info.txt')
    lines_to_skip = 1
    filenames = sorted(next(os.walk(target_path), (None, None, []))[2])
    filenames = [filename for filename in filenames if filename.endswith('.txt') and 'info' not in filename]
    if os.path.exists(current_info):
        new_info_path = os.path.join(processed_result_path, 'target_info', f'{os.path.basename(target_path)}.txt')
        try:
            shutil.copy(current_info, new_info_path)
        except FileNotFoundError:
            print("Could not copy info file")
    else:
        lines_to_skip = 0
        single_result_path = os.path.join(target_path, filenames[0])
        with open(single_result_path, "r") as f:
            for line in f:
                if line.startswith('REMARK') or line.startswith('HEADER'):
                    lines_to_skip += 1
                if line.startswith('RESULT'):
                    break
        new_info_path = None
    concatenated_file_path = os.path.join(target_path, 'concatenated.csv')
    with open(concatenated_file_path, "wb") as fout:
        with open(os.path.join(target_path, filenames[0]), "rb") as f:
            if not os.path.exists(current_info):
                for i in range(lines_to_skip-1):
                    next(f)
            fout.write(f.read())
        for filename in filenames[1:-1]:
            with open(os.path.join(target_path, filename), "rb") as f:
                for i in range(lines_to_skip):
                    next(f)  # skip the header
                fout.write(f.read())
    return concatenated_file_path, new_info_path


def compress_one_target(target, output_path, output_pose, keep_pdb_number, result_path, default_cf, ef_percentage,
                        keep_only_best_conf_pose, htpvs, root_software_path, remove_duplicate_path, original_data_path,
                        number_ligands_to_keep=None, db_path=None, ligand_type=None):
    print(f"Processing results for target: {target}")
    target_result_path = os.path.join(result_path, target)
    result_csv_path, new_info_path = merge_csv(target_result_path, output_path)
    try:
        df, need_for_ef = load_df(result_csv_path, ligand_type)
    except:
        print("ERROR target: ", target)
    write_remove_duplicate(target, remove_duplicate_path, df)
    if df.height != 0:
        if need_for_ef:
            csv_path = os.path.join(root_software_path, "deps", f"dud_e_ligand_count.csv")
            ligand_count_df = pl.scan_csv(csv_path)
            result = ligand_count_df.filter(pl.col("Target name") == target).select(["actives", "decoys"]).collect()
            active_ligand_count = result.row(0)[0]
            decoy_ligand_count = result.row(0)[1]
            first_score = df.row(1, named=True)['CF']
            if first_score == 0 or first_score == default_cf:
                ef = "Error"
            else:
                ef = calculate_ef(df, active_ligand_count, decoy_ligand_count, ef_percentage)
            if new_info_path:
                with open(new_info_path, 'a') as f:
                    f.write(f"REMARK EF = {ef}\n")
        if not htpvs:
            shutil.copyfile(result_csv_path, os.path.join(original_data_path, f"{target}.csv"))
        if output_pose:
            if keep_only_best_conf_pose:
                ligand_list = ligand_list_sorted
            else:
                ligand_list = result_list
            if keep_pdb_number.isdigit():
                delete_ligands(keep_pdb_number, ligand_list, output_path, target)
            else:
                new_dir = os.path.join(output_path, target + "_ligand_poses")
                shutil.make_archive(new_dir, 'tar', os.path.join(".", "temp", "ligand_poses", target))
                shutil.rmtree(os.path.join(".", "temp", "ligand_poses", target))
                os.mkdir(os.path.join(".", "temp", "ligand_poses", target))
        if need_for_ef:
            return {"Name": target, "EF": ef}, need_for_ef
        if number_ligands_to_keep is not None:
            table_name = os.path.splitext(os.path.basename(db_path))[0]
            name_top_n = df.select(pl.col('Name')).head(number_ligands_to_keep)['Name'].to_list()
            query_db(db_path, name_top_n, table_name, os.path.join(os.path.dirname(remove_duplicate_path), 'top_smiles.csv'))
        return None

    else:
        print("No results found to analyse")


def main(specific_target, htpvs, folder_end, top_n=None, db_path=None, custom_result_path=None):
    # TODO:do an error check on job slurm files for out of time
    root_software_path = Path(__file__).resolve().parents[1]
    os.chdir(root_software_path)
    result_path = "./temp/results/"
    if custom_result_path:
        result_path = os.path.dirname(custom_result_path)
        ligand_type = 'ligand'
    else:
        ligand_type = None
    print(result_path)
    config_file = "./deps/config.txt"
    ef_percentage = 0.01
    original_data_path = None

    params_dict, suffix = get_params(config_file)
    n_rotations = params_dict['ROTATIONS_PER_AXIS']
    output_pose = params_dict["OUTPUT_POSE"]
    number_to_keep = params_dict["KEPT_PDB_NUMBER"]
    default_cf = params_dict["DEFAULT_CF"]
    keep_only_best_conf_pose = params_dict['BEST_CONFORMER']

    output_path = get_output_name(suffix, n_rotations, folder_end)
    print("output path: ", output_path)
    info_output_path = os.path.join(output_path, "target_info")
    if not os.path.exists(info_output_path):
        os.mkdir(info_output_path)
    remove_duplicate = os.path.join(output_path, "duplicate_ligands_removed")
    if not os.path.isdir(remove_duplicate):
        os.mkdir(remove_duplicate)
    if not htpvs:
        original_data_path = os.path.join(output_path, 'original_data')
        if not os.path.exists(original_data_path):
            os.mkdir(original_data_path)
    if specific_target:
        targets = [specific_target]
    elif custom_result_path:
        targets = [os.path.basename(custom_result_path)]
    else:
        targets = sorted(next(os.walk(result_path))[1])
    result_list = []

    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     ef_list = list(executor.map(compress_one_target, targets, repeat(output_path), repeat(output_pose),
    #                                 repeat(number_to_keep), repeat(result_path), repeat(default_cf),
    #                                 repeat(ef_percentage), repeat(keep_only_best_conf_pose), repeat(htpvs),
    #                                 repeat(root_software_path), repeat(remove_duplicate), repeat(original_data_path),
    #                                 repeat(top_n), repeat(db_path), repeat(ligand_type)))

    ef_list = compress_one_target(targets[0], output_path, output_pose, number_to_keep, result_path, default_cf,
                                  ef_percentage, keep_only_best_conf_pose, htpvs, root_software_path, remove_duplicate,
                                  original_data_path, top_n, db_path, ligand_type)
    ef_list = [ef_list]

    all_none = all(item is None for item in ef_list)
    if not all_none:
        for ef in ef_list:
            if ef is not None:
                result_list.append(ef[0])
        save_ef_csv(result_list, os.path.join(output_path, "EF.csv"))
    reset_result_folder(result_path, params_dict)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--specific_target', default=None,
                        help='specify one target to be compressed and analysed, otherwise None')
    parser.add_argument('-htpvs', '--high_throughput', action='store_true',
                        help='will only output duplicate_ligands_removed')
    parser.add_argument('-z', '--zip_result', default=None, action='store_true',
                        help='zip the output')
    parser.add_argument('-f', '--folder_name', default=None,
                        help='add name of target folder in output')
    parser.add_argument('-csv', '--csv_top_n', default=None, type=int,
                        help='will get the smiles for the top n ligands')
    parser.add_argument('-db', '--db_path', default=None, help='path to db for -csv option')
    parser.add_argument('-cr', '--custom_result_path', default=None, help='path to results')

    args = parser.parse_args()
    spe_target = args.specific_target
    htpvs = args.high_throughput
    folder_end = args.folder_name
    top_n = args.csv_top_n
    db_path = args.db_path
    custom_result_path = args.custom_result_path
    if top_n is not None and db_path is None:
        parser.error('-csv flag also needs -db to specify path to database')
    main(spe_target, htpvs, folder_end, top_n, db_path, custom_result_path)


if __name__ == "__main__":
    get_args()


#TODO: see how to count ligands, seoarate case of large scale docking
# TODO: could save memory if no need for raw output