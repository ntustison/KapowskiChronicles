#! /bin/sh

DATA_DIR=./
TEMPLATE_DIR=${DATA_DIR}/template/
OUT_DIR=${DATA_DIR}/output2/

# set ANTSPATH
export ANTSPATH=/Users/ntustison/Pkg/ANTs/bin/bin/

bash ${ANTSPATH}antsCorticalThickness.sh -d 2 \
  -a ${DATA_DIR}/IXI002-Guys-0828-T1_slice90.nii.gz \
  -a ${DATA_DIR}/IXI002-Guys-0828-T2_slice90.nii.gz \
  -e ${TEMPLATE_DIR}template_slice80.nii.gz \
  -m ${TEMPLATE_DIR}template_cerebrum_mask_slice80.nii.gz \
  -p ${TEMPLATE_DIR}prior%d_slice80.nii.gz \
  -f ${TEMPLATE_DIR}template_extraction_mask_slice80.nii.gz \
  -w 0.25 \
  -n 3 \
  -k 0 \
  -o ${OUT_DIR}antsCorticalThicknessMultivariateExample_
