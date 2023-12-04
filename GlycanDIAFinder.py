from collections import defaultdict
from statistics import mean
import csv
from configparser import ConfigParser
import timeit
from matchms.importing import load_from_mzxml, load_from_mzml
from matchms import set_matchms_logger_level
import numpy as np

import matplotlib.pyplot as plt

from collections import defaultdict
# import argparse
import gc
from function import *



#%%
if __name__ == "__main__":
    set_matchms_logger_level("ERROR")


    # read configuration from config.ini
    cfg = ConfigParser()
    cfg.read("./config.ini")
    cfg_dict = dict(cfg.items("config"))
    print("cfg_dict", cfg_dict)

    input_path = cfg_dict["input_path"]
    output_path = cfg_dict["output_path"]
    ms_list_name = cfg_dict["ms_list_name"]
    polarity = cfg_dict["polarity"]
    charge = int(cfg_dict["max_charge"])
    ppm_ms1 = float(cfg_dict["ms1_mass_error_ppm"])
    ppm_ms2 = float(cfg_dict["ms2_mass_error_ppm"])
    adduct = cfg_dict["adduct"]
    match_count_ms2 = int(cfg_dict["min_matched_counts"])

    if "align_tolerance_snum" in cfg_dict.keys():
        align_delta = int(cfg_dict["align_tolerance_snum"])
    else:
        align_delta = 100


    # extract info from ms_list.csv
    dataset_list, cpdAddon_note_list, note_cpdAddon_list_dict, ms2_mass_dict = extrac_dataset_info(input_path, fn=ms_list_name)

    print("dataset_list", dataset_list)
    print("cpdAddon_note_list", cpdAddon_note_list)
    print("note_cpdAddon_list_dict", note_cpdAddon_list_dict)
    print("ms2_mass_dict", ms2_mass_dict)

    ### generate decoy library
    decoy = cfg_dict["decoy_mode"]
    if decoy:
        print('generating decoy library, this may be up to 100x slower...')
        target_seq_score = []
        decoy_seq_score = []
        target_int_score = []
        decoy_int_score = []
        num_decoy_sets = 100  # Number of decoy sets to generate
        decoy_spectra = defaultdict(list) # cpd_addon : list of spectrum array
        target_spectra = defaultdict(list) # cpd_addon : list of spectrum array
        for cpd_addon in ms2_mass_dict.keys():
            # convert string to list of float
            ms2_mass_list = ms2_mass_dict[cpd_addon].split(" ")
            ms2_mass_list = [float(ms2_mass) for ms2_mass in ms2_mass_list]
            target_spectra[cpd_addon] = ms2_mass_list
            for decoy_set in range(num_decoy_sets):
                decoy_spectra[cpd_addon].append(add_random_mass_shifts(ms2_mass_list, 0.5))





    comp_root_row_dict = defaultdict(list) # cpd-AddonMass : list of row 
    comp_root_head = ["Cpd-AddonMass", "Note"]
    subtype_root_row_dict = defaultdict(list) # subtype : list of row
    subtype_root_head = ["Subtype"]

    start = timeit.default_timer()

    # extract retention time
    if "min_time_min" in cfg_dict.keys():
        min_rt = float(cfg_dict["min_time_min"])
    else:
        min_rt = 0
    if "max_time_min" in cfg_dict.keys():
        max_rt = float(cfg_dict["max_time_min"])
    else:
        max_rt = float("inf")
    if "debug_mode" in cfg_dict.keys():
        debug_mode = True
    else:
        debug_mode = False

    for dataset in dataset_list:
        # input full path
        input_fn = input_path + "/" + dataset + ".mzXML"
        

        # processing the input file
        if debug_mode:
            print("starting processing the input file:", input_fn)
        # extract the format of the input file
        fn_list = input_fn.split(".")
        fn_format = fn_list[-1]
        
        # processing the MS1 file
        print("processing MS1...")
        if fn_format == "mzXML":
            spectrums_ms1_all_list = list(load_from_mzxml(input_fn, ms_level=1))
        else:
            raise("Unsupported file format.")

        # filter the spectrum_ms1 by rention time range
        spectrums_ms1_list = []
        for spectrum_ms1 in spectrums_ms1_all_list:
            rt_spec_all_ms1 = spectrum_ms1.metadata["retention_time"]
            if rt_spec_all_ms1 >= min_rt and rt_spec_all_ms1 <= max_rt:
                spectrums_ms1_list.append(spectrum_ms1)

        if debug_mode:
            print("Number of spectrums in MS1:", len(spectrums_ms1_list))
        
        
        # processing the MS2 file
        print("processing MS2...")
        if fn_format == "mzXML":
            spectrums_ms2_all_list = list(load_from_mzxml(input_fn, ms_level=2))
        else:
            raise("Unsupported file format.")

        spectrums_ms2_list = []
        for spectrum_ms2 in spectrums_ms2_all_list:
            rt_spec_all_ms2 = spectrum_ms2.metadata["retention_time"]
            if rt_spec_all_ms2 >= min_rt and rt_spec_all_ms2 <= max_rt:
                spectrums_ms2_list.append(spectrum_ms2)

        if debug_mode:
            print("Number of spectrums in MS2:", len(spectrums_ms2_list))
        
        # iteratively extract spectrum index from spectrums in MS2
        # list of spectrum_idx of a specific precursor_mz
        precMZ_spectID_dict = defaultdict(list) # precursor_mz : list of spectrum idx
        for spec_ms2_idx, spectrum_ms2 in enumerate(spectrums_ms2_list):
            spect_precursor_mz = spectrum_ms2.metadata["precursor_mz"]

            # precursor_mz_set.add(spectrum_ms2.metadata["precursor_mz"])
            precMZ_spectID_dict[spect_precursor_mz].append(spec_ms2_idx)
        
        if debug_mode:
            print("# of unique precursor_mz", len(precMZ_spectID_dict))
        

        comp_root_head += ["Height (MS1)_Data File " + dataset, "Relative Height (MS1 %)_Data File" + dataset, "Height (MS2)_Data File " + dataset, "Relative Height (MS2 %)_Data File" + dataset, "Max Cov. Sequence (%) " + dataset, "Avg Cov. Sequence (%) " + dataset, "Max Cov. Int (%) " + dataset, "Avg Cov. Int (%) " + dataset]
        subtype_root_head += ["Height (MS1)_Data File " + dataset, "Relative Height (MS1 %)_Data File" + dataset, "Height (MS2)_Data File " + dataset, "Relative Height (MS2 %)_Data File" + dataset, "Max Cov. Sequence (%) " + dataset, "Avg Cov. Sequence (%) " + dataset, "Max Cov. Int (%) " + dataset, "Avg Cov. Int (%) " + dataset]
        # aggregated csv information
        isomers_row_list = []
        # cpd-AddonMass : total height_ms1
        cpd_tot_h_ms1_dict = defaultdict(float)
        # cpd-AddonMass : total height_ms2
        cpd_tot_h_ms2_dict = defaultdict(float)
        # cpd-AddonMass : list of max match%
        cpd_match_perc_list_dict = defaultdict(list)
        cpd_int_perc_list_dict = defaultdict(list)
        # note : total height_ms1
        note_tot_h_ms1_dict = defaultdict(float)
        # note : total height_ms2
        note_tot_h_ms2_dict = defaultdict(float)

        for cpd_addon, note in cpdAddon_note_list:

            cpd = cpd_addon.split("-")[0]
            addon_mass = cpd_addon.split("-")[1]

            print("cpd_addon", cpd_addon)

            # create output paths if not existed
            output_fd = output_path + "/" + dataset + "/" + cpd_addon
            creat_path(output_fd)


            # create configuration list for an Arguments class
            if "min_rel_height" in cfg_dict.keys():
                threshold = cfg_dict["min_rel_height"]
            else:
                threshold = 0.001
            if "min_height" in cfg_dict.keys():
                min_height = cfg_dict["min_height"]
            else:
                min_height = 5000
            if "min_mass" in cfg_dict.keys():
                min_mass = cfg_dict["min_mass"]
            else:
                min_mass = 0
            if "max_mass" in cfg_dict.keys():
                max_mass = cfg_dict["max_mass"]
            else:
                max_mass = float("inf")


            if "flex_mode" in cfg_dict.keys():
                flex_mode = cfg_dict["flex_mode"]
            else:
                flex_mode = "false"
            if "charge_range" in cfg_dict.keys():
                charge_range = cfg_dict["charge_range"]
            else:
                charge_range = None


            if "max_aligned_record_ms2" in cfg_dict.keys():
                max_aligned_record_ms2 = int(cfg_dict["max_aligned_record_ms2"])
            else:
                max_aligned_record_ms2 = float("inf")
            
            args = Arguments(input_fn, output_fd, ms2_mass_dict[cpd_addon], min_height, threshold, charge, charge_range, polarity, ppm_ms1, ppm_ms2, cpd, adduct, addon_mass, match_count_ms2, note, min_mass, max_mass, min_rt, max_rt, flex_mode, debug_mode, max_aligned_record_ms2, align_delta)


            # os.system("python3 align_peaks_matchms_batch.py " + config + " | tee " + output_fd + "/debug.log")
            
            align_peaks_matchms_batch(args, spectrums_ms1_list, spectrums_ms2_list, precMZ_spectID_dict)
        
            plt.close('all')
            if decoy:
                tmp = get_target_decoy_scores(args, spectrums_ms1_list, spectrums_ms2_list, precMZ_spectID_dict,decoy_spectra,cpd_addon)
                target_seq_score += tmp[0]
                decoy_seq_score += tmp[1]
                target_int_score += tmp[2]
                decoy_int_score += tmp[3]

            # read Glycan_isomers.csv for each cpd_addon
            with open(output_fd + "/Glycan_isomers.csv", "r") as f:
                reader = csv.reader(f)
                # match_percent_max = 0
                # match_percent_tot = 0
                # peak_num = 0
                for line_idx, row in enumerate(reader):
                    # print(row, type(row))
                    # remove the header of the csv file
                    if line_idx == 0:
                        continue
                    
                    #isomers
                    isomers_row_list.append(row)

                    # composition
                    cpd_tot_h_ms1_dict[cpd_addon] += float(row[5])
                    cpd_tot_h_ms2_dict[cpd_addon] += float(row[6])
                    cpd_match_perc_list_dict[cpd_addon].append(float(row[8].split("%")[0]))
                    cpd_int_perc_list_dict[cpd_addon].append(float(row[9].split("%")[0]))
                    # subtype
                    note_tot_h_ms1_dict[note] += float(row[5])
                    note_tot_h_ms2_dict[note] += float(row[6])

        # delete ms data
        # print("freeing MS data...")
        savespace = True
        if savespace:
            
            del spectrums_ms1_all_list
            del spectrums_ms1_list
            del spectrums_ms2_all_list
            del spectrums_ms2_list
            del precMZ_spectID_dict
        gc.collect()

        # generate isomers form
        with open(output_path + "/" + dataset + "/Glycan_isomers.csv", "w") as f:
            writer = csv.writer(f)
            # write headline
            csv_head = ["Cpd-AddOnMass", "Note", "RT", "m/z", "charge", "Height(MS1)", "Height(MS2)", "Matched Counts","Cov. Sequence", "Cov. Int" ]
            writer.writerow(csv_head)
            # print(isomers_row_list)

            for row in isomers_row_list:
                writer.writerow(row)

        # generate composition form

        # calculate total height for all cpd of all AddOnMass in ms1 and ms2
        tot_h_ms1 = 0
        tot_h_ms2 = 0
        for cpd_addon, note in cpdAddon_note_list:
            tot_h_ms1 += cpd_tot_h_ms1_dict[cpd_addon]
            tot_h_ms2 += cpd_tot_h_ms2_dict[cpd_addon]
        
        with open(output_path + "/" + dataset + "/Glycan_composition.csv", "w") as f:
            writer = csv.writer(f)
            # write headline
            csv_head = ["Cpd-AddonMass", "Note", "Height (MS1)", "Relative Height (MS1 %)", "Height (MS2)", "Relative Height (MS2 %)", "Max Cov. Sequence (%)", "Avg Cov. Sequence (%)", "Max Cov. Int (%)", "Avg Cov. Int (%)"]
            writer.writerow(csv_head)

            for cpd_addon, note in cpdAddon_note_list:

                # print(dataset, cpd_addon)
                match_percent_max = 0
                match_percent_avg = 0
                int_percent_max = 0
                int_percent_avg = 0

                if cpd_match_perc_list_dict[cpd_addon]:
                    match_percent_max = max(cpd_match_perc_list_dict[cpd_addon])
                    match_percent_avg = mean(cpd_match_perc_list_dict[cpd_addon])
                    int_percent_max = max(cpd_int_perc_list_dict[cpd_addon])
                    int_percent_avg = mean(cpd_int_perc_list_dict[cpd_addon])


                writer.writerow([cpd_addon, note, cpd_tot_h_ms1_dict[cpd_addon], cpd_tot_h_ms1_dict[cpd_addon] / tot_h_ms1 * 100, cpd_tot_h_ms2_dict[cpd_addon], cpd_tot_h_ms2_dict[cpd_addon] / tot_h_ms2 * 100, match_percent_max, match_percent_avg, int_percent_max, int_percent_avg])

                comp_root_row_dict[cpd_addon] += [cpd_tot_h_ms1_dict[cpd_addon], cpd_tot_h_ms1_dict[cpd_addon] / tot_h_ms1 * 100, cpd_tot_h_ms2_dict[cpd_addon], cpd_tot_h_ms2_dict[cpd_addon] / tot_h_ms2 * 100, match_percent_max, match_percent_avg, int_percent_max, int_percent_avg]
        
        # subtype
        with open(output_path + "/" + dataset + "/Glycan_subtype.csv", "w") as f:
            writer = csv.writer(f)
            # write headline
            csv_head = ["subtype", "Height (MS1)", "Relative Height (MS1 %)", "Height (MS2)", "Relative Height (MS2 %)"]
            writer.writerow(csv_head)

            for note in note_tot_h_ms1_dict.keys():
                writer.writerow([note, note_tot_h_ms1_dict[note], note_tot_h_ms1_dict[note] / tot_h_ms1 * 100, note_tot_h_ms2_dict[note], note_tot_h_ms2_dict[note]/ tot_h_ms2 * 100])

                subtype_root_row_dict[note] += [note_tot_h_ms1_dict[note], note_tot_h_ms1_dict[note] / tot_h_ms1 * 100, note_tot_h_ms2_dict[note], note_tot_h_ms2_dict[note]/ tot_h_ms2 * 100]

    
    # Root Glycan_composition form
    with open(output_path + "/Glycan_composition_combined.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(comp_root_head)

        for cpd_addon, note in cpdAddon_note_list:
            row = [cpd_addon, note] + comp_root_row_dict[cpd_addon]
            writer.writerow(row)
    
    # Root Glycan_subtype form
    with open(output_path + "/Glycan_subtype_combined.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(subtype_root_head)

        for note in subtype_root_row_dict.keys():
            row = [note] + subtype_root_row_dict[note]
            writer.writerow(row)
    
    if decoy:
        # get the score index whose score is less than 1
        index = [i for i in range(len(decoy_int_score)) if decoy_int_score[i] < 1]
        decoy_int_score = [decoy_int_score[i] for i in index]
        decoy_seq_score = [decoy_seq_score[i] for i in index]
        # estimate threshold
        seq_threshold = get_threshold(target_seq_score, decoy_seq_score, fdr_threshold=0.02)
        int_threshold = get_threshold(target_int_score, decoy_int_score, fdr_threshold=0.02)
        print("seq_threshold@0.2fdr:", seq_threshold)
        print("int_threshold@0.2fdr:", int_threshold)

        import pickle
        print('saving scores in pickle...')
        with open(output_path + "target_decoy_score_batch2.pickle", "wb") as f:
            pickle.dump({'target_int':target_int_score,
                        "target_seq":target_seq_score,
                        "decoy_int":decoy_int_score,
                        "decoy_seq":decoy_seq_score}, f)
    end = timeit.default_timer()
    print("Running time: %s Seconds"%(end-start))















