Usage:

In python script directory:

`$ python make_bias_dark_opticflat_biasdarksub_science.py path_to_file/filelist.txt type`

In target directory:

`$ python path_to_script/make_bias_dark_opticflat_biasdarksub_science.py filelist.txt type`

type: bias, dark, opticflat, biasdarksub, science, biasdarksub2science

This python script defines functions to make bias.fits, dark300.fits, dark600.fits, optic flats, biasdarksub_xxx.fits and science_xxx.fits. For unnecessary re-running this script, it will find and use existing bias/dark/optic/twilight flat... from environment variable $_90PRIME_BIAS_DARK_FLATS_DIR


Note: "science" here means corrected frame, not the final masked science frame for stacking
