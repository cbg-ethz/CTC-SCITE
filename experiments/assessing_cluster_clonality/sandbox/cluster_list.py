#!/usr/bin python3

from pathlib import Path
import pandas as pd

input_folder = Path.home() / "Documents" / "projects" / "CTC_backup" / "input_folder"

files = input_folder.glob("splittingSummary*")



data_frames = [pd.read_csv(file) for file in files]

for data_frame, filename in zip(data_frames, input_folder.glob("splittingSummary*")):
    sample_name = "_".join(filename.stem.split(".")[0].split("_")[1:])
    data_frame["Sample Name"] = [sample_name] * data_frame.shape[0]

data_frames_valid = [data_frame for data_frame in data_frames if data_frame.shape[1] == 6 ]
full_data_frame = pd.concat(data_frames_valid)

full_data_frame.to_csv(input_folder / "splittingSummary_full.tsv", sep='\t')
