'''
This script is for ST data SME modified - modified to adjust to error in Read10X from stlearn and scanpy
'''
#!/usr/bin/env python
# conding:utf-8 -*-
import argparse
import stlearn as st
import scanpy as sc
import numpy as np
from numpy import random,mat
from pathlib import Path
import pandas as pd
from scipy import io,sparse
import os

from PIL import Image
from anndata import AnnData
import json
import matplotlib.pyplot as plt
from matplotlib.image import imread
import logging as log

def ME_normalize_new(inDir,outDir,sample):
    print (sample, "start SME normalize")

    #read data
    data=read_10X_custom(path = inDir)
    data.var_names_make_unique()
    data.layers['raw_count']=data.X
    #tile data
    TILE_PATH=Path(os.path.join(outDir,'{0}_tile'.format(sample)))
    TILE_PATH.mkdir(parents=True,exist_ok=True)
    
    #tile morphology
    st.pp.tiling(data,TILE_PATH,crop_size=40)
    st.pp.extract_feature(data)

    ###process data
    st.pp.normalize_total(data)
    st.pp.log1p(data)

    #gene pca dimention reduction
    st.em.run_pca(data,n_comps=50,random_state=0)
    
    #stSME to normalise log transformed data
    st.spatial.SME.SME_normalize(data, use_data="raw",weights =  "weights_matrix_gd_md")

    #convert SME_norm data to sparesmatrix
    raw_SME_normalized = mat(data.obsm['raw_SME_normalized'])
    raw_SME_normalizedA = sparse.csr_matrix(raw_SME_normalized)
    print ("matrix convert ok!")
    
    io.mmwrite(os.path.join(outDir,'{0}_raw_SME_normalizeA.mtx'.format(sample)),raw_SME_normalizedA)
    print("Morphology adjusted ok!")

    return raw_SME_normalizedA
   
def read_10X_custom(path):
    # establish parameters
    path = Path(path)
    count_file = "filtered_feature_bc_matrix.h5"
    genome = None
    library_id = None
    load_images = True
    image_path = None
    quality = "lowres"

    # create scanpy AnnData object
    adata = sc.read_10x_h5(path / count_file, genome=genome)
    adata.uns["spatial"] = dict()

    # library_id and attrs
    from h5py import File
    with File(path / count_file, mode="r") as f:
        attrs = dict(f.attrs)
    if  library_id is None:
        library_id = attrs.pop("library_ids")[0]
    
    adata.uns["spatial"][library_id] = dict()

    # establish files
    tissue_positions_file = (
        path / "spatial/tissue_positions.csv"
        if (path / "spatial/tissue_positions.csv").exists()
        else path / "spatial/tissue_positions_list.csv"
    )   

    # image processing
    if load_images:
        files = dict(
            tissue_positions_file=tissue_positions_file,
            scalefactors_json_file=path / "spatial/scalefactors_json.json",
            #hires_image=path / "spatial/tissue_hires_image.png",
            lowres_image=path / "spatial/tissue_lowres_image.png",
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not f.exists():
                if any(x in str(f) for x in ["lowres_image"]):
                    log.warning(
                        f"You seem to be missing an image file.\n"
                        f"Could not find '{f}'."
                    )
                else:
                    raise OSError(f"Could not find line105 '{f}'")

        adata.uns["spatial"][library_id]["images"] = dict()
        for res in ["lowres"]:
            try:
                adata.uns["spatial"][library_id]["images"][res] = imread(
                    str(files[f"{res}_image"])
                )
            except Exception:
                print(str(files[f"{res}_image"]))
                raise OSError(f"Could not find line114'{res}_image'")

        # read json scalefactors
        adata.uns["spatial"][library_id]["scalefactors"] = json.loads(
            files["scalefactors_json_file"].read_bytes()
        )

        adata.uns["spatial"][library_id]["metadata"] = {
            k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
            for k in ("chemistry_description", "software_version")
            if k in attrs
        }

        # read coordinates
        positions = pd.read_csv(files["tissue_positions_file"], header=None)
        positions.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
        ]
        positions.index = positions["barcode"]

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm["spatial"] = (
            adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]]
            .to_numpy()
            .astype(int)
        )
        adata.obs.drop(
            columns=["barcode", "pxl_row_in_fullres", "pxl_col_in_fullres"],
            inplace=True,
        )

        # put image path in uns
        if image_path is not None:
            # get an absolute path
            image_path = str(Path(image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
                image_path
            )

    adata.var_names_make_unique()

    if library_id is None:
        library_id = list(adata.uns["spatial"].keys())[0]

    if quality == "fulres":
        image_coor = adata.obsm["spatial"]
        img = plt.imread(image_path, 0)
        adata.uns["spatial"][library_id]["images"]["fulres"] = img
    else:
        scale = adata.uns["spatial"][library_id]["scalefactors"][
            "tissue_" + quality + "_scalef"
        ]
        image_coor = adata.obsm["spatial"] * scale

    adata.obs["imagecol"] = image_coor[:, 0]
    adata.obs["imagerow"] = image_coor[:, 1]
    adata.uns["spatial"][library_id]["use_quality"] = quality

    adata.obs["array_row"] = adata.obs["array_row"].astype(int)
    adata.obs["array_col"] = adata.obs["array_col"].astype(int)
    adata.obsm["spatial"] = adata.obsm["spatial"].astype("int64")

    return adata