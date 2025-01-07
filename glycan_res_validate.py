import glob
import os
import pandas as pd
def results_validation(input_folder):
    isomers_file_path = glob.glob(os.path.join(input_folder,  'test*_isomers.csv'), recursive=True)[0]
    result = []
    if not isomers_file_path:
        print(f"No isomers files found in '{input_folder}'.")
    else:
        isomers_df = pd.read_csv(isomers_file_path)
        for index,glycan in enumerate(isomers_df["Cpd-AddOnMass"]):
            glycan_num = glycan.split("-")[0].split("_")
            h_num = int(glycan_num[0])
            n_num = int(glycan_num[1])
            f_num = int(glycan_num[2])
            a_num = int(glycan_num[3])
            ff_num = f_num
            while(ff_num>=2):
                ff_num = ff_num-2
                aa_num = a_num+1
                new_glycan_com = f'{h_num}_{n_num}_{ff_num}_{aa_num}-0'
                matching_rows = isomers_df[isomers_df["Cpd-AddOnMass"].str.contains(new_glycan_com, na=False)]
                reference_RT = isomers_df.iloc[index]["RT"]
                for _, row in matching_rows.iterrows():
                    if row["RT"] == reference_RT:
                        result.append(glycan)
    return result





if __name__ == "__main__":
    input_folder= "ExampleDataset/Results/stagger_Nglycan_ExampleData/"
    glycan_fake = results_validation(input_folder)
    print(glycan_fake)