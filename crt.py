import click, pickle, inquirer

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

import struct
import uproot
import awkward as ak

def save_triggers_with_uproot(triggers, filename):
    # Construir listas para cada campo
    modules = []
    time50MHz_list = []
    raw_backend_time_list = []
    timestamp_list = []
    nhits_list = []
    channels_list = []
    adcs_list = []

    for trig in triggers:
        modules.append(trig['module'])
        time50MHz_list.append(trig['time50MHz'])
        raw_backend_time_list.append(trig['raw_backend_time'])
        timestamp_list.append(trig['timestamp'])

        nhits_list.append(len(trig['hits']))
        channels_list.append([h['offline_channel'] for h in trig['hits']])
        adcs_list.append([h['adc'] for h in trig['hits']])

    # Criar uma árvore (TTree) com arrays (Awkward Array para listas variáveis)
    array_data = {
        "module": np.array(modules, dtype=np.int32),
        "time50MHz": np.array(time50MHz_list, dtype=np.uint64),
        "raw_backend_time": np.array(raw_backend_time_list, dtype=np.uint64),
        "timestamp": np.array(timestamp_list, dtype=np.uint64),
        "nhits": np.array(nhits_list, dtype=np.int32),
        "channels": ak.Array(channels_list),
        "adcs": ak.Array(adcs_list),
    }

    with uproot.recreate(filename) as f:
        f["trigger_tree"] = array_data


def map_offline_channel(module_num, channel):
  
    if module_num in [14,15,8,9,10,11,4,5,30,31,24,25,26,27,20,21]:
        if channel < 32:
            offline_channel = (31 - channel) * 2
        else:
            offline_channel = (63 - channel) * 2 + 1
    else:
        if channel < 32:
            offline_channel = channel * 2
        else:
            offline_channel = (channel - 32) * 2 + 1

    if offline_channel % 2 == 0:
        offline_channel += 1
    else:
        offline_channel -= 1

    return offline_channel


def write_permission(directory_path: str) -> bool:
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

    rucio_filepath = [
        "root://eospublic.cern.ch:1094//eos/experiment/neutplatform/protodune/dune/hd-protodune/13/63/np04hd_raw_run028612_0000_dataflow2_datawriter_0_20240812T221909.hdf5"]
   
    for my_file in rucio_filepath:

        temporal_copy_directory = '.'

        try:
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

        except:
            fUsedXRootD = False

        print("oi")
        h5_file = HDF5RawDataFile(my_file)
        records = h5_file.get_all_record_ids()

        wvfm_index = 0
        det="HD_CRT"
        
        # Create a root file
        file_output = "./" + my_file.split('/')[-1].removesuffix(".hdf5")+"_crt.root"
        print(file_output)

        fChannelMap = [
            24, 25, 30, 18, 15, 7, 13, 12, 11, 10, 6, 14,
            19, 31, 26, 27, 22, 23, 29, 17, 8, 0, 3, 2,
            5, 4, 1, 9, 16, 28, 20, 21
        ] #32 channels

        triggers = []

        for i, r in enumerate(tqdm(records)):
            geo_ids = list(h5_file.get_geo_ids_for_subdetector(
                r, detdataformats.DetID.string_to_subdetector(det)))
            for gid in geo_ids:
                try:
                    frag = h5_file.get_frag(r, gid)
                    trig = h5_file.get_trh(r)
                except:
                    continue  # melhor que None aqui, senão continua com valores quebrados

                det_link = 0xffff & (gid >> 48)
                det_slot = 0xffff & (gid >> 32)
                det_crate = 0xffff & (gid >> 16)
                det_id = 0xffff & gid
                subdet = detdataformats.DetID.Subdetector(det_id)
                det_name = detdataformats.DetID.subdetector_to_string(subdet)

                #print(f"link: { det_link}")
                #print(f"slot: { det_slot}")
                #print(f"crate: { det_crate}")
                #print(f"id: { det_id}")

                fragment_bytes = frag.get_data_bytes()

                start_index = fragment_bytes.find(b'M')  # busca o byte 0x4D
                if start_index == -1:
                    print("Magic Byte'M' not found. Jumping fragment.")
                    continue

                if len(fragment_bytes) < start_index + 16:
                    print("Very little fragment. SKIPPING IT.")
                    continue

                header_format = "<IIQ"
                fields = struct.unpack_from(header_format, fragment_bytes)

                timestamp = fields[2]

                header_format = "<BBHIQ"
                header_size = 16
                hit_size = 4

                #print(fragment_bytes[0:16].hex())
                #print(list(fragment_bytes[0:16])) 
                
                header = struct.unpack_from(header_format, fragment_bytes, start_index)
                magic, nhit, module_num, fifty_mhz_time, raw_back_time = header

                """ print( magic, nhit, module_num, fifty_mhz_time, raw_back_time)"""
                #input("pause") 

                if magic != ord('M'):
                    print(f"Magic header bytenot found: expecting: 'M', found: {magic}")
                    continue

                hits = []
                offset_hits = start_index + header_size
                for i_hit in range(nhit):
                    hit_bytes = fragment_bytes[offset_hits + i_hit * hit_size: offset_hits + (i_hit + 1) * hit_size]
                    if len(hit_bytes) < hit_size:
                        print("Incomplete hit bytes, jumping.")
                        continue

                    hit_magic, channel, adc = struct.unpack("<BBh", hit_bytes)

                    #print(channel)
                    #input()

                    #print(hit_magic)
                    if hit_magic != ord('H'):
                        print(f"WRONG: Hit {i_hit} magic byte : {hit_magic}")
                        continue

                    offline_channel = map_offline_channel(module_num, channel)
                    hits.append({'offline_channel': offline_channel, 'adc': adc})

                try:
                    offline_module_num = fChannelMap[module_num]
                except IndexError:
                    print(f"Module {module_num} outside channel map.")
                    continue

                # Agora montamos o CRT::Trigger no estilo do C++
                trigger = {
                    'module': offline_module_num,
                    'time50MHz': fifty_mhz_time,
                    'raw_backend_time': raw_back_time,
                    'timestamp': timestamp,
                    'hits': hits
                }

                triggers.append(trigger)
                
        save_triggers_with_uproot(triggers, file_output)
        print(f"Arquivo ROOT salvo com uproot: {file_output}")
        # Exemplo: imprimir todos os triggers deste evento
        #for trig in triggers:
            #print("Trigger montado:", trig)

            
if __name__ == "__main__":
    main()
