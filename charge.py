import click, pickle, inquirer



from waffles.utils.utils import print_colored
import waffles.input.raw_hdf5_reader as reader
import re

import os, io, click, subprocess, stat, math, shlex
from array import array
from tqdm import tqdm
import numpy as np
from XRootD import client
from typing import List, Optional
from hdf5libs import HDF5RawDataFile
import tempfile

import uproot # import daqdataformats
from daqdataformats import FragmentType
#from rawdatautils.unpack.daphne import *
from rawdatautils.unpack.utils  import *
import detdataformats
import fddetdataformats

from rawdatautils.unpack.dataclasses import *
import rawdatautils.unpack.wibeth

from multiprocessing import Pool, current_process, cpu_count


def write_permission(directory_path: str) -> bool:
    """This function gets the path to a directory, and checks
    whether the running process has write permissions in such
    directory.

    Parameters
    ----------
    directory_path: str
        Path to an existing directory. If the directory does
        not exist, an exception is raised.

    Returns
    ----------
    bool
        True if the the running process has write permissions
        in the specified directory. False if else."""

    try:
        with tempfile.TemporaryFile(
            dir=directory_path
        ):
            pass

        return True

    except FileNotFoundError:
        raise we.NonExistentDirectory(
            we.GenerateExceptionMessage(
                1,
                'write_permission()',
                f"The specified directory ({directory_path}) does not exist."
            )
        )

    except (OSError, PermissionError):
        return False



def main():
    debug=True
    run= "28210"
    if run is None: 
        q = [ inquirer.Text("run", message="Please provide the run(s) number(s) to be analysed, separated by commas:)")]
        run_list = inquirer.prompt(q)["run"].split(",")
    else:
        run_list = run.split(",")
    
    for r in run_list:
        rucio_filepath = f"/eos/experiment/neutplatform/protodune/experiments/ProtoDUNE-II/PDS_Commissioning/waffles/1_rucio_paths/{str(r).zfill(6)}.txt"
        rucio_filepath="028210.txt"
        if debug: 
            print_colored(f"Processing {str(r).zfill(6)}...", color="DEBUG")

    
    filepaths = reader.get_filepaths_from_rucio(rucio_filepath)
    #print(filepaths)
    for my_file in filepaths[18:]:
        print_colored(my_file,color="DEBUG")
        temporal_copy_directory = '/tmp'

        if "/eos" not in my_file and "/nfs" not in my_file and "/afs" not in my_file:
            print("Using XROOTD")

            if write_permission(temporal_copy_directory):

                subprocess.call(
                    shlex.split(f"xrdcp {my_file} {temporal_copy_directory}"),
                    shell=False
                )
                fUsedXRootD = True

                my_file = os.path.join(
                    temporal_copy_directory,
                    my_file.split('/')[-1]
                )

            else:
                raise Exception(
                    GenerateExceptionMessage(
                        1,
                        'WaveformSet_from_hdf5_file()',
                        f"Attempting to temporarily copy {my_file} into "
                        f"{temporal_copy_directory}, but the current process "
                        f"has no write permissions there. Please specify a "
                        "valid directory."
                    )
                )

        else:
            fUsedXRootD = False

        h5_file = HDF5RawDataFile(my_file)
        records = h5_file.get_all_record_ids()

        wvfm_index = 0
        det="HD_TPC"
        
        # Create a root file
        file_output = "/afs/cern.ch/work/g/gbotogos/data/charge/"+my_file.split('/')[-1].removesuffix(".hdf5")+"_charge"
        print(file_output)
        
        crate_geo_list = []
        slot_geo_list = []
        stream_geo_list = []
        timestamp_list = []
        channel_list = []
        waveform_list = []
        offline_channel_list=[]
        apa_list=[]
        plane_list=[]
        index_save=0
        index_name=0

        ch_map = detchannelmaps.make_map("PD2HDChannelMap")

        print(len(tqdm(records)))
        for i, r in enumerate(tqdm(records)):
            geo_ids = list(h5_file.get_geo_ids_for_subdetector(
                r, detdataformats.DetID.string_to_subdetector(det)))
            
            for gid in geo_ids:
                try:
                    frag = h5_file.get_frag(r, gid)
                    trig = h5_file.get_trh(r)
                    
                except:
                    None

                #wf = fddetdataformats.WIBFrame(frag.get_data())
                #wh = wf.get_wib_header()
            
                crate_from_geo = 0xffff & (gid >> 16)
                slot_from_geo = 0xffff & (gid >> 32)
                stream_from_geo = 0xffff & (gid >> 48)
                
                adcs = WIBEthUnpacker.unpacker.np_array_adc(frag).T
                #for a in adcs:
                #    print(a)
                #print((adcs[1]))

                #input("Pressione Enter para continuar...")

                timestamps=WIBEthUnpacker.unpacker.np_array_timestamp(frag)[0] #cada deltat sao 32 unidades de timestamps ( eh muito ?)
                n_frame = WIBEthUnpacker.unpacker.get_n_frames(frag)
                
                #offline_ch_num_dict = [ch_map.get_offline_channel_from_crate_slot_fiber_chan(wh.crate_no, wh.slot_no, wh.fiber_no, c) for c in channels]
                offline_ch_num_dict = [ ch_map.get_offline_channel_from_crate_slot_stream_chan(crate_from_geo, slot_from_geo, stream_from_geo, c) for c in range(WIBEthUnpacker.N_CHANNELS_PER_FRAME) ]
                offline_ch_num_dict = offline_ch_num_dict
                planes =[ch_map.get_plane_from_offline_channel(c) for c in offline_ch_num_dict]
                APAs = [ch_map.get_tpc_element_from_offline_channel(c) for c in offline_ch_num_dict]

                for j in range(WIBEthUnpacker.N_CHANNELS_PER_FRAME):
                    crate_geo_list.append(crate_from_geo)
                    slot_geo_list.append(slot_from_geo)
                    stream_geo_list.append(stream_from_geo)
                    timestamp_list.append(timestamps)
                    #channel_list.append(channels[j])
                    waveform_list.append(adcs[j])
                    #print(len(adcs[j]))
                    apa_list.append(APAs[j])
                    plane_list.append(planes[j])
                    offline_channel_list.append(offline_ch_num_dict[j])      
                    
                if(index_save==5):
                        
                    # Converter listas para arrays NumPy
                    crate_geo = np.array(crate_geo_list, dtype=np.float32)
                    crate_geo_list.clear()
                    slot_geo = np.array(slot_geo_list, dtype=np.float32)
                    slot_geo_list.clear()
                    stream_geo = np.array(stream_geo_list, dtype=np.float32)
                    stream_geo_list.clear()
                    timestamp = np.array(timestamp_list, dtype=np.uint64)
                    timestamp_list.clear()
                    #channel = np.array(channel_list, dtype=np.float32)
                    #channel_list.clear()
                    offline_channel = np.array(offline_channel_list, dtype=np.float32)
                    offline_channel_list.clear()
                    apa = np.array(apa_list, dtype=str)
                    apa_list.clear()
                    plane = np.array(plane_list, dtype=np.float32)
                    plane_list.clear()
                    waveform_array = [np.array(w, dtype=np.float32) for w in waveform_list]  # np.array(waveform_list, dtype=np.float32)  # 2D array 
                    waveform_list.clear()

                    # Criar e salvar a TTree no arquivo ROOT
                    with uproot.recreate(file_output+f"_{index_name}.root") as f:
                        f["my_tree"] = {
                            "crate_geo": crate_geo,
                            "slot_geo": slot_geo,
                            "stream_geo": stream_geo,
                            "timestamp": timestamp,
                            #"channel": channel,
                            "wf": waveform_array,
                            "apa": apa,
                            "off_channel": offline_channel,
                            "plane": plane
                        }
                    index_name=index_name+1
                    index_save=0
                index_save=index_save+1

if __name__ == "__main__":
    main()
