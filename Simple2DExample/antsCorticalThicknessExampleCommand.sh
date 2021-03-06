#! /bin/sh

DATA_DIR=./
TEMPLATE_DIR=${DATA_DIR}/template/
OUT_DIR=${DATA_DIR}/output/

# set ANTSPATH
export ANTSPATH=/Users/ntustison/Pkg/ANTs/bin/bin/

bash ${ANTSPATH}antsCorticalThickness.sh -d 2 \
  -a ${DATA_DIR}KKI2009-08-MPRAGE_slice165.nii.gz \
  -e ${TEMPLATE_DIR}template_slice80.nii.gz \
  -m ${TEMPLATE_DIR}template_cerebrum_mask_slice80.nii.gz \
  -p ${TEMPLATE_DIR}prior%d_slice80.nii.gz \
  -f ${TEMPLATE_DIR}template_extraction_mask_slice80.nii.gz \
  -w 0.25 \
  -n 3 \
  -k 0 \
  -b Aristotle[1] \
  -l 1[1.0,0.5] \
  -l 2[1.0,0.5] \
  -l 3[1.0,0.5] \
  -l 4[5.0,0.25] \
  -o ${OUT_DIR}antsCorticalThicknessExample_
