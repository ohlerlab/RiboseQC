# sed -n '/host/,/test/p' recipes/riboseqc/meta.yaml | grep -e '-' | sed 's/.*\ -//' | xargs conda create -n buildtest
# cp /fast/users/dharnet_m/work/Applications/conda_RiboseQC/{build*,meta*} recipes/riboseqc/
# mkdir tmp;cd tmp
# cat ../recipes/riboseqc/meta.yaml | grep url: | sed 's/url: //' | xargs wget